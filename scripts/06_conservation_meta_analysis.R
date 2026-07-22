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

# For each organelle, find human orthologs that appear as DEG-associated in multiple species
# and run Fisher's combined p-value
# We need p-values per gene per species - use the DEG p-values

# Get per-gene p-values from DEG data
deg_pvals <- deg[, .(species, gene_id, gene_id_base = sub("\\.\\d+$", "", gene_id), 
                      pvalue, padj, is_deg, deg_direction, human_ortholog)]

# For human, set human_ortholog to gene_id_base
deg_pvals[species == "Homo sapiens", human_ortholog := gene_id_base]

# For each organelle, collect p-values per species for the shared orthologous gene set
conservation_results <- list()

for (org in unique(org_gene_human_mapped$organelle)) {
  org_genes <- org_gene_human_mapped[organelle == org]
  
  # Get unique human orthologs for this organelle
  human_orthos <- unique(org_genes$human_ensembl)
  
  # For each human ortholog, collect p-values across species
  for (ho in human_orthos) {
    # Get p-values for this human ortholog in each species
    sp_pvals <- deg_pvals[human_ortholog == ho & !is.na(pvalue) & pvalue > 0 & pvalue < 1, 
                          .(species, pvalue, deg_direction)]
    
    if (nrow(sp_pvals) >= 2) {
      # Fisher's combined p-value: -2 * sum(log(p))
      fisher_stat <- -2 * sum(log(sp_pvals$pvalue))
      fisher_p <- pchisq(fisher_stat, df = 2 * nrow(sp_pvals), lower.tail = FALSE)
      
      # Check direction consistency
      dirs <- unique(sp_pvals$deg_direction)
      dirs <- dirs[dirs != "ns"]
      direction_conserved <- if (length(dirs) == 1) dirs[1] else "mixed"
      
      conservation_results[[length(conservation_results) + 1]] <- data.table(
        organelle = org,
        human_ensembl = ho,
        n_species = nrow(sp_pvals),
        species_list = paste(sp_pvals$species, collapse = ","),
        fisher_stat = fisher_stat,
        fisher_p = fisher_p,
        fisher_padj = NA_real_,
        direction_conserved = direction_conserved
      )
    }
  }
}

conservation_dt <- rbindlist(conservation_results)
# BH correction within each organelle
conservation_dt[, fisher_padj := p.adjust(fisher_p, method = "BH"), by = organelle]

fwrite(conservation_dt, '/mnt/shared-workspace/shared/results/meta/organelle_conservation_fisher.csv')

cat("\n=== Conservation meta-analysis results ===\n")
cat("Total ortholog-organelle tests:", nrow(conservation_dt), "\n")
cat("Significant (Fisher padj < 0.05):", sum(conservation_dt$fisher_padj < 0.05), "\n")

cat("\n=== Conserved organelle responses (Fisher padj < 0.05, >=3 species, same direction) ===\n")
conserved_sig <- conservation_dt[fisher_padj < 0.05 & n_species >= 3 & direction_conserved != "mixed"]
cat("Count:", nrow(conserved_sig), "\n")
if (nrow(conserved_sig) > 0) print(conserved_sig[order(fisher_padj)])

cat("\n=== Summary by organelle (>=3 species, any direction) ===\n")
org_summary <- conservation_dt[n_species >= 3, .(
  n_orthologs = .N,
  n_sig = sum(fisher_padj < 0.05),
  n_conserved = sum(fisher_padj < 0.05 & direction_conserved != "mixed"),
  n_up = sum(fisher_padj < 0.05 & direction_conserved == "up"),
  n_down = sum(fisher_padj < 0.05 & direction_conserved == "down"),
  min_padj = min(fisher_padj)
), by = organelle][order(-n_sig)]
print(org_summary)
fwrite(org_summary, '/mnt/shared-workspace/shared/results/meta/organelle_conservation_summary.csv')

# --- R code block 2 ---
suppressPackageStartupMessages({ library(data.table) })

conservation_dt <- fread('/mnt/shared-workspace/shared/results/meta/organelle_conservation_fisher.csv')

# Fix: count UNIQUE species (not total entries)
conservation_dt[, unique_species := sapply(strsplit(species_list, ","), function(x) length(unique(x)))]

# Recompute summary with unique species count
org_summary_clean <- conservation_dt[unique_species >= 3, .(
  n_orthologs = .N,
  n_sig = sum(fisher_padj < 0.05),
  n_conserved = sum(fisher_padj < 0.05 & direction_conserved != "mixed"),
  n_up = sum(fisher_padj < 0.05 & direction_conserved == "up"),
  n_down = sum(fisher_padj < 0.05 & direction_conserved == "down"),
  min_padj = min(fisher_padj)
), by = organelle][order(-n_conserved)]

cat("=== Clean conservation summary (>=3 unique species) ===\n")
print(org_summary_clean)
fwrite(org_summary_clean, '/mnt/shared-workspace/shared/results/meta/organelle_conservation_summary.csv')

# Key finding: conserved downregulated mitochondrial genes
mito_down <- conservation_dt[organelle == "Mitochondrion" & fisher_padj < 0.05 & direction_conserved == "down" & unique_species >= 3]
cat("\n=== Conserved downregulated mitochondrial genes (>=3 species) ===\n")
cat("Count:", nrow(mito_down), "\n")
if (nrow(mito_down) > 0) print(mito_down[order(fisher_padj), .(human_ensembl, unique_species, fisher_padj, direction_conserved)][1:20])

# Top conserved responses across all organelles
cat("\n=== Top 20 conserved orthologs (>=4 unique species, same direction) ===\n")
top_conserved <- conservation_dt[unique_species >= 4 & fisher_padj < 0.05 & direction_conserved != "mixed"][order(fisher_padj)]
print(top_conserved[1:20, .(organelle, human_ensembl, unique_species, fisher_padj, direction_conserved)])

# Save the final conservation table
fwrite(conservation_dt[, .(organelle, human_ensembl, unique_species, species_list, fisher_stat, fisher_p, fisher_padj, direction_conserved)],
      '/mnt/shared-workspace/shared/results/meta/organelle_conservation_final.csv')

