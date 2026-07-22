# ======================================================================
# Step 3: Orthology Mapping (OrthoDB v12 + babelgene)
# 
# Description: Build human-anchored orthology matrix using OrthoDB v12 and babelgene
# 
# Inputs: OrthoDB v12 bulk files, babelgene R package
# Outputs: results/orthology/orthology_matrix_wide.csv, unified_orthology_matrix.csv
# 
# Language: R
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- R code block 1 ---
suppressPackageStartupMessages(library(babelgene))

# Test babelgene ortholog mapping
# Map a human gene to mouse ortholog
hs_to_mm <- orthologs(genes = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"),
                      species = "Homo sapiens",
                      to_species = "Mus musculus")
print(head(hs_to_mm))

# Check what species babelgene supports
print("\nSupported species:")
print(head(babelgene::species))

# Map human to all 5 other species
for (to_sp in c("Mus musculus", "Drosophila melanogaster", "Caenorhabditis elegans", 
                "Saccharomyces cerevisiae", "Arabidopsis thaliana")) {
  ortho <- orthologs(genes = c("ENSG00000000003", "ENSG00000000419", "ENSG00000141510"),
                     species = "Homo sapiens", to_species = to_sp)
  cat(sprintf("\nHuman -> %s: %d orthologs found\n", to_sp, nrow(ortho)))
  if (nrow(ortho) > 0) print(head(ortho[,c("ortholog_gene","symbol","species","ortholog_species","evidence","support")]))
}

# --- R code block 2 ---
suppressPackageStartupMessages(library(babelgene))
# Check function signature
args(orthologs)
# Try correct syntax
hs_to_mm <- orthologs(genes = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"),
                      species = "Homo sapiens",
                      human = FALSE)
print(head(hs_to_mm))
cat("\nColumns:", paste(colnames(hs_to_mm), collapse=", "), "\n")

# --- R code block 3 ---
suppressPackageStartupMessages(library(babelgene))
# Check function args
print(args(orthologs))

# Try with species names from the list
hs_to_mm <- orthologs(genes = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"),
                      species = "human")
print(head(hs_to_mm))
cat("\nColumns:", paste(colnames(hs_to_mm), collapse=", "), "\n")
cat("Unique species:", unique(hs_to_mm$species), "\n")

# --- R code block 4 ---
suppressPackageStartupMessages(library(babelgene))
# "Homo sapiens" is not in the list but "human" should work... let me check the species list
sp <- babelgene::species
print(head(sp))
print(class(sp))
# Try "Mus musculus"
hs_to_mm <- orthologs(genes = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"),
                      species = "Mus musculus")
print(head(hs_to_mm))

# --- R code block 5 ---
suppressPackageStartupMessages(library(babelgene))
library(data.table)

# Get all human genes from our DEG data
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv')
human_genes <- unique(deg[species == "Homo sapiens", gene_id])
cat("Human genes in DEG data:", length(human_genes), "\n")

# Map human genes to all 5 other species (babelgene supports these)
target_species <- c("Mus musculus", "Drosophila melanogaster", "Caenorhabditis elegans", "Saccharomyces cerevisiae")

all_orthologs <- list()
for (to_sp in target_species) {
  cat(sprintf("Mapping human -> %s...\n", to_sp))
  ortho <- orthologs(genes = human_genes, species = to_sp, human = TRUE, min_support = 1, top = TRUE)
  ortho$target_species <- to_sp
  all_orthologs[[to_sp]] <- ortho
  cat(sprintf("  Found %d ortholog pairs\n", nrow(ortho)))
}

# Combine
ortho_all <- rbindlist(all_orthologs)
cat("\nTotal ortholog pairs:", nrow(ortho_all), "\n")

# Save
fwrite(ortho_all, '/mnt/shared-workspace/shared/orthodb/babelgene_orthologs.csv')
cat("Saved to babelgene_orthologs.csv\n")

# Summary
cat("\n=== Ortholog mapping summary ===\n")
print(ortho_all[, .N, by = target_species])
cat("\nHuman genes with orthologs in each species:\n")
for (sp in target_species) {
  n <- ortho_all[target_species == sp, uniqueN(human_ensembl)]
  cat(sprintf("  %s: %d human genes have orthologs\n", sp, n))
}

# --- R code block 6 ---
suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

# Use Ensembl Plants to get Arabidopsis orthologs
# Try Ensembl Plants
tryCatch({
  ensembl_plants <- useEnsembl(biomart = "plants_mart", dataset = "athaliana_eg_gene")
  cat("Connected to Ensembl Plants\n")
  
  # Get Arabidopsis genes and their human orthologs
  ath_genes <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv')
  ath_genes <- unique(ath_genes[species == "Arabidopsis thaliana", gene_id])
  cat("Arabidopsis genes in DEG data:", length(ath_genes), "\n")
  
  # Get orthologs to human
  ath_ortho <- getBM(
    attributes = c("ensembl_gene_id", "name_1006", "ensembl_gene_id_hsapiens"),
    filters = "ensembl_gene_id",
    values = ath_genes,
    mart = ensembl_plants
  )
  cat("Arabidopsis orthologs found:", nrow(ath_ortho), "\n")
  print(head(ath_ortho))
  
  fwrite(ath_ortho, '/mnt/shared-workspace/shared/orthodb/arabidopsis_orthologs.csv')
}, error = function(e) {
  cat("Ensembl Plants error:", e$message, "\n")
  cat("Will use OrthoDB API instead\n")
})

# --- R code block 7 ---
suppressPackageStartupMessages({
  library(data.table)
  library(babelgene)
})

# Load all DEGs
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv')

# === 1. Babelgene orthologs (human -> mouse, fly, worm, yeast) ===
babel_ortho <- fread('/mnt/shared-workspace/shared/orthodb/babelgene_orthologs.csv')
cat("Babelgene orthologs:", nrow(babel_ortho), "\n")

# === 2. Arabidopsis orthologs (via OrthoDB OGs) ===
ath_ortho <- fread('/mnt/shared-workspace/shared/orthodb/arabidopsis_to_human_orthologs.csv')
# Split pipe-separated human orthologs and take first
ath_ortho[, human_ensembl := sapply(strsplit(human_ortholog, "\\|"), function(x) x[1])]
ath_ortho <- ath_ortho[human_ensembl != "" & !is.na(human_ensembl)]
cat("Arabidopsis orthologs:", nrow(ath_ortho), "\n")

# === Build unified orthology matrix ===
# Format: human_ensembl | species | species_gene_id | og_id (where available)

# Babelgene: human_ensembl -> target species gene
babel_unified <- babel_ortho[, .(human_ensembl, target_species, 
                                  species_gene_id = ensembl,
                                  source = "babelgene")]
setnames(babel_unified, "target_species", "species")

# Arabidopsis: human_ensembl -> Arabidopsis gene
ath_unified <- ath_ortho[, .(human_ensembl, 
                              species = "Arabidopsis thaliana",
                              species_gene_id = gene_id,
                              source = "orthodb")]
ath_unified[, human_ensembl := sub("\\.\\d+$", "", human_ensembl)]  # strip version

# Combine
ortho_matrix <- rbind(babel_unified, ath_unified)
cat("\n=== Unified orthology matrix ===\n")
cat("Total ortholog pairs:", nrow(ortho_matrix), "\n")
print(ortho_matrix[, .N, by = .(species, source)])

# Get unique human genes that have orthologs in multiple species
ortho_wide <- dcast(ortho_matrix, human_ensembl ~ species, 
                    value.var = "species_gene_id", 
                    fun.aggregate = function(x) x[1])
cat("\nOrthology matrix (wide):", nrow(ortho_wide), "human genes\n")
cat("Columns:", paste(colnames(ortho_wide), collapse=", "), "\n")

# Count how many species each human gene has orthologs in
ortho_wide[, n_species := rowSums(!is.na(.SD)), .SDcols = c("Mus musculus","Drosophila melanogaster","Caenorhabditis elegans","Saccharomyces cerevisiae","Arabidopsis thaliana")]
cat("\nDistribution of ortholog count per human gene:\n")
print(table(ortho_wide$n_species))

# Save
fwrite(ortho_matrix, '/mnt/shared-workspace/shared/orthodb/unified_orthology_matrix.csv')
fwrite(ortho_wide, '/mnt/shared-workspace/shared/orthodb/orthology_matrix_wide.csv')
cat("\nSaved unified orthology matrix\n")

