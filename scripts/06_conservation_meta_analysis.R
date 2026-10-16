# ======================================================================
# Step 6: Cross-Species Conservation Meta-Analysis
# 
# Description: Fisher combined probability test for organelle conservation across species
# 
# Inputs: results/enrichment/gocc_enrichment_classified.csv, results/deg/all_degs_with_orthologs.csv
# Outputs: results/meta/organelle_conservation_fisher.csv, organelle_conservation_summary.csv, conserved_organelles.csv
# 
# Language: R
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- R code block 1 ---
suppressPackageStartupMessages({ library(data.table) })

# Load DEGs with orthologs
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gene_id)]

# Load organelle classification
enr_exploded <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_classified.csv')

# For each organelle category, collect the DEG gene lists per species
# and run Fisher's combined p-value on the shared orthologous gene set

# First, get the gene IDs from each enriched GOCC term
# The geneID column in the enrichment results has the genes
enr <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')

# For each organelle category, get all DEG genes that are in enriched GOCC terms
# Map these to human orthologs, then test conservation

# Build organelle -> genes mapping per species
organelle_genes <- list()
for (i in 1:nrow(enr_exploded)) {
  row <- enr_exploded[i]
  org <- row$organelle
  sp <- row$species
  dir <- row$direction
  go_id <- row$ID
  
  # Get genes from the original enrichment result
  go_row <- enr[ID == go_id & species == sp & direction == dir]
  if (nrow(go_row) > 0) {
    genes <- strsplit(go_row$geneID, "/")[[1]]
    key <- paste(org, sp, dir, sep = "||")
    if (is.null(organelle_genes[[key]])) organelle_genes[[key]] <- c()
    organelle_genes[[key]] <- unique(c(organelle_genes[[key]], genes))
  }
}

# Convert to data.table
org_gene_dt <- rbindlist(lapply(names(organelle_genes), function(key) {
  parts <- strsplit(key, "||", fixed = TRUE)[[1]]
  data.table(organelle = parts[1], species = parts[2], direction = parts[3],
             gene_id = organelle_genes[[key]])
}))

cat("=== Organelle gene counts ===\n")
print(org_gene_dt[, .N, by = .(organelle, species, direction)][order(organelle, species)])

# Map genes to human orthologs
# For human, gene_id is already human ensembl
# For other species, use the ortholog mapping
ortho_long <- fread('/mnt/shared-workspace/shared/orthodb/unified_orthology_matrix.csv')

# Build species_gene -> human_ensembl mapping
ortho_map <- ortho_long[, .(species, species_gene_id, human_ensembl)]
ortho_map[, human_ensembl := sub("\\.\\d+$", "", human_ensembl)]

# Merge organelle genes with orthologs
org_gene_dt[, gene_id_base := sub("\\.\\d+$", "", gene_id)]
org_gene_human <- merge(org_gene_dt, ortho_map, 
                        by.x = c("species", "gene_id_base"),
                        by.y = c("species", "species_gene_id"),
                        all.x = TRUE)

# For human, map to itself
human_mask <- org_gene_human$species == "Homo sapiens"
org_gene_human[human_mask, human_ensembl := gene_id_base]

cat("\n=== Organelle genes mapped to human orthologs ===\n")
org_gene_human_mapped <- org_gene_human[!is.na(human_ensembl)]
print(org_gene_human_mapped[, .N, by = .(organelle, species, direction)])

# --- Conservation test, as specified in METHODS.md ---------------------------
# For each ortholog-organelle pair:
#   * one p-value per species: the smallest p among that species' rows
#     (datasets and paralogs mapping to the human ortholog), Sidak-corrected
#     for the number of rows, 1 - (1 - p_min)^k;
#   * Fisher's combined test over those per-species p-values (>= 2 species);
#   * BH-adjusted across all ortholog-organelle tests;
#   * "conserved" if Fisher FDR < 0.05, DEG in >= 3 species, and >= 80% of the
#     DEG species share one direction. A species whose DEG rows disagree in
#     direction counts as a DEG species agreeing with neither.
# An earlier version pooled every row as an independent Fisher input, counted
# species that merely had a p-value toward ">= 3 species", called a gene "down"
# if every DEG call was down (even a single one), dropped p = 1 rows, and
# adjusted within each organelle.

deg_pvals <- deg[, .(species, gene_id_base = sub("\\.\\d+$", "", gene_id),
                     pvalue, is_deg = as.character(is_deg) %in% c("TRUE", "True", "true"),
                     deg_direction, human_ortholog)]
deg_pvals[species == "Homo sapiens", human_ortholog := gene_id_base]
deg_pvals <- deg_pvals[!is.na(pvalue) & !is.na(human_ortholog) & human_ortholog != ""]
deg_pvals[, pvalue := pmin(pmax(pvalue, 1e-300), 1)]

per_species <- deg_pvals[, {
  dirs <- unique(deg_direction[is_deg])
  list(k = .N,
       p_species = -expm1(.N * log1p(-min(pvalue))),   # Sidak: 1 - (1 - p_min)^k
       species_direction = if (length(dirs) == 0) "ns" else if (length(dirs) == 1) dirs else "mixed")
}, by = .(human_ortholog, species)]
per_species[, p_species := pmax(p_species, 1e-300)]

candidates <- unique(org_gene_human_mapped[, .(organelle, human_ensembl)])
tests <- merge(candidates, per_species, by.x = "human_ensembl", by.y = "human_ortholog",
               allow.cartesian = TRUE)

conservation_dt <- tests[, {
  n_tested <- .N
  fisher_stat <- -2 * sum(log(p_species))
  n_up <- sum(species_direction == "up")
  n_down <- sum(species_direction == "down")
  n_mixed <- sum(species_direction == "mixed")
  n_deg <- n_up + n_down + n_mixed
  agreement <- if (n_deg > 0) max(n_up, n_down) / n_deg else NA_real_
  list(n_species_tested = n_tested,
       n_deg_species = n_deg,
       n_up_species = n_up,
       n_down_species = n_down,
       n_mixed_species = n_mixed,
       direction_agreement = agreement,
       species_list = paste(sort(species), collapse = ","),
       fisher_stat = fisher_stat,
       fisher_p = pchisq(fisher_stat, df = 2 * n_tested, lower.tail = FALSE))
}, by = .(organelle, human_ensembl)][n_species_tested >= 2]

conservation_dt[, fisher_padj := p.adjust(fisher_p, method = "BH")]
conservation_dt[, direction_conserved := fifelse(
  !is.na(direction_agreement) & direction_agreement >= 0.8,
  fifelse(n_up_species >= n_down_species, "up", "down"), "mixed")]
conservation_dt[n_deg_species == 0, direction_conserved := "ns"]
conservation_dt[, conserved := fisher_padj < 0.05 & n_deg_species >= 3 & direction_conserved %in% c("up", "down")]
setorder(conservation_dt, fisher_p)

fwrite(conservation_dt, '/mnt/shared-workspace/shared/results/meta/organelle_conservation_fisher.csv')

cat("\n=== Conservation meta-analysis results ===\n")
cat("Ortholog-organelle tests:", nrow(conservation_dt), "\n")
cat("Fisher FDR < 0.05:", sum(conservation_dt$fisher_padj < 0.05), "\n")
cat("Conserved (FDR < 0.05, DEG in >= 3 species, >= 80% same direction):",
    sum(conservation_dt$conserved), "\n")

org_summary <- conservation_dt[, .(
  n_orthologs = .N,
  n_fisher_sig = sum(fisher_padj < 0.05),
  n_deg_3plus = sum(n_deg_species >= 3),
  n_conserved = sum(conserved),
  n_up = sum(conserved & direction_conserved == "up"),
  n_down = sum(conserved & direction_conserved == "down"),
  min_padj = min(fisher_padj)
), by = organelle][order(-n_conserved, min_padj)]
print(org_summary)
fwrite(org_summary, '/mnt/shared-workspace/shared/results/meta/organelle_conservation_summary.csv')

cat("\n=== Conserved mitochondrial genes ===\n")
print(conservation_dt[organelle == "Mitochondrion" & conserved,
                      .(human_ensembl, n_species_tested, n_deg_species, n_up_species, n_down_species,
                        n_mixed_species, fisher_padj, direction_conserved)])

fwrite(conservation_dt[conserved == TRUE],
       '/mnt/shared-workspace/shared/results/meta/organelle_conservation_final.csv')
