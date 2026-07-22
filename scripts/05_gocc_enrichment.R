# ======================================================================
# Step 5: GOCC Subcellular Enrichment + Organelle Classification
# 
# Description: Run clusterProfiler enrichGO (CC ontology) per species; classify terms into organelle categories
# 
# Inputs: results/deg/all_degs_with_orthologs.csv, species-specific org.db packages
# Outputs: results/enrichment/gocc_enrichment_all.csv, gocc_enrichment_classified.csv, organelle_enrichment_summary.csv
# 
# Language: R
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- R code block 1 ---
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(data.table)
  library(GO.db)
})

# Load DEGs with orthologs
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gene_id)]

# Species-specific org.db packages and key types
species_config <- list(
  "Homo sapiens" = list(orgdb = "org.Hs.eg.db", keytype = "ENSEMBL", id_strip_version = TRUE),
  "Mus musculus" = list(orgdb = "org.Mm.eg.db", keytype = "ENSEMBL", id_strip_version = TRUE),
  "Drosophila melanogaster" = list(orgdb = "org.Dm.eg.db", keytype = "ENSEMBL", id_strip_version = TRUE),
  "Caenorhabditis elegans" = list(orgdb = "org.Ce.eg.db", keytype = "ENSEMBL", id_strip_version = TRUE),
  "Saccharomyces cerevisiae" = list(orgdb = "org.Sc.sgd.db", keytype = "ORF", id_strip_version = FALSE),
  "Arabidopsis thaliana" = list(orgdb = "org.At.tair.db", keytype = "TAIR", id_strip_version = FALSE)
)

# Run GOCC enrichment per species (separately for up and down DEGs)
all_enrichment <- list()

for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  cat(sprintf("\n=== %s ===\n", sp))
  
  # Get DEGs for this species
  sp_degs <- deg[species == sp & is_deg == TRUE]
  if (nrow(sp_degs) == 0) {
    cat("  No DEGs, skipping\n")
    next
  }
  
  # Get all genes (universe) for this species
  sp_all <- deg[species == sp]
  
  # Prepare gene IDs
  prepare_ids <- function(ids, strip_version) {
    ids <- unique(na.omit(ids))
    if (strip_version) ids <- sub("\\.\\d+$", "", ids)
    return(ids)
  }
  
  deg_ids_up <- prepare_ids(sp_degs[deg_direction == "up", gene_id], cfg$id_strip_version)
  deg_ids_down <- prepare_ids(sp_degs[deg_direction == "down", gene_id], cfg$id_strip_version)
  universe_ids <- prepare_ids(sp_all$gene_id, cfg$id_strip_version)
  
  cat(sprintf("  Up DEGs: %d, Down DEGs: %d, Universe: %d\n", length(deg_ids_up), length(deg_ids_down), length(universe_ids)))
  
  # Load org.db
  orgdb <- getNamespace(cfg$orgdb)
  
  for (direction in c("up", "down")) {
    gene_list <- if (direction == "up") deg_ids_up else deg_ids_down
    if (length(gene_list) < 5) {
      cat(sprintf("  %s: too few DEGs (%d), skipping\n", direction, length(gene_list)))
      next
    }
    
    tryCatch({
      ego <- enrichGO(gene = gene_list,
                      universe = universe_ids,
                      OrgDb = orgdb,
                      keyType = cfg$keytype,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = FALSE)
      
      if (nrow(as.data.frame(ego)) > 0) {
        res <- as.data.frame(ego)
        res$species <- sp
        res$direction <- direction
        res$osd_datasets <- paste(unique(sp_degs$osd), collapse=",")
        all_enrichment[[paste(sp, direction, sep="_")]] <- res
        cat(sprintf("  %s: %d enriched GOCC terms\n", direction, nrow(res)))
      } else {
        cat(sprintf("  %s: no enriched terms\n", direction))
      }
    }, error = function(e) {
      cat(sprintf("  %s: ERROR %s\n", direction, e$message))
    })
  }
}

# Combine all enrichment results
if (length(all_enrichment) > 0) {
  enrichment_all <- rbindlist(all_enrichment, fill = TRUE)
  fwrite(enrichment_all, '/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')
  cat(sprintf("\n=== Total enriched GOCC terms: %d ===\n", nrow(enrichment_all)))
  cat("\nEnrichment counts per species/direction:\n")
  print(enrichment_all[, .N, by = .(species, direction)])
} else {
  cat("\nNo enrichment results!\n")
}

# --- R code block 2 ---
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(data.table)
  library(GO.db)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Dm.eg.db)
  library(org.Ce.eg.db)
  library(org.Sc.sgd.db)
  library(org.At.tair.db)
})

# Load DEGs
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gene_id)]

# Species config: use the actual OrgDb objects
species_config <- list(
  "Homo sapiens" = list(orgdb = org.Hs.eg.db, keytype = "ENSEMBL", id_strip_version = TRUE),
  "Mus musculus" = list(orgdb = org.Mm.eg.db, keytype = "ENSEMBL", id_strip_version = TRUE),
  "Drosophila melanogaster" = list(orgdb = org.Dm.eg.db, keytype = "ENSEMBL", id_strip_version = TRUE),
  "Caenorhabditis elegans" = list(orgdb = org.Ce.eg.db, keytype = "ENSEMBL", id_strip_version = TRUE),
  "Saccharomyces cerevisiae" = list(orgdb = org.Sc.sgd.db, keytype = "ORF", id_strip_version = FALSE),
  "Arabidopsis thaliana" = list(orgdb = org.At.tair.db, keytype = "TAIR", id_strip_version = FALSE)
)

all_enrichment <- list()

for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  cat(sprintf("\n=== %s ===\n", sp))
  
  sp_degs <- deg[species == sp & is_deg == TRUE]
  if (nrow(sp_degs) == 0) { cat("  No DEGs\n"); next }
  
  sp_all <- deg[species == sp]
  
  prepare_ids <- function(ids, strip_version) {
    ids <- unique(na.omit(ids))
    if (strip_version) ids <- sub("\\.\\d+$", "", ids)
    return(ids)
  }
  
  deg_ids_up <- prepare_ids(sp_degs[deg_direction == "up", gene_id], cfg$id_strip_version)
  deg_ids_down <- prepare_ids(sp_degs[deg_direction == "down", gene_id], cfg$id_strip_version)
  universe_ids <- prepare_ids(sp_all$gene_id, cfg$id_strip_version)
  
  cat(sprintf("  Up: %d, Down: %d, Universe: %d\n", length(deg_ids_up), length(deg_ids_down), length(universe_ids)))
  
  for (direction in c("up", "down")) {
    gene_list <- if (direction == "up") deg_ids_up else deg_ids_down
    if (length(gene_list) < 5) { cat(sprintf("  %s: too few (%d)\n", direction, length(gene_list))); next }
    
    tryCatch({
      ego <- enrichGO(gene = gene_list, universe = universe_ids,
                      OrgDb = cfg$orgdb, keyType = cfg$keytype,
                      ont = "CC", pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE)
      if (nrow(as.data.frame(ego)) > 0) {
        res <- as.data.frame(ego)
        res$species <- sp; res$direction <- direction
        all_enrichment[[paste(sp, direction, sep="_")]] <- res
        cat(sprintf("  %s: %d terms\n", direction, nrow(res)))
      } else { cat(sprintf("  %s: none\n", direction)) }
    }, error = function(e) { cat(sprintf("  %s: ERROR %s\n", direction, e$message)) })
  }
}

if (length(all_enrichment) > 0) {
  enrichment_all <- rbindlist(all_enrichment, fill = TRUE)
  fwrite(enrichment_all, '/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')
  cat(sprintf("\n=== Total: %d enriched GOCC terms ===\n", nrow(enrichment_all)))
  print(enrichment_all[, .N, by = .(species, direction)])
}

# --- R code block 3 ---
suppressPackageStartupMessages({
  library(data.table)
  library(GO.db)
})

# Load enrichment results
enr <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')

# Define broad organelle categories with their GOCC term patterns
# Based on SubcellulaRVis categories and GO hierarchy
organelle_categories <- list(
  "Nucleus" = c("nucleus", "nuclear", "chromatin", "nucleolus", "nucleoplasm", "nuclear envelope", "nuclear pore"),
  "Mitochondrion" = c("mitochondr", "mitochondrial", "mitochondria", "oxidative phosphorylation", "respiratory chain", "mitochondrial matrix", "mitochondrial membrane", "mitochondrial inner membrane", "mitochondrial outer membrane"),
  "Endoplasmic Reticulum" = c("endoplasmic reticulum", "endoplasmic", "ER membrane", "rough ER", "smooth ER", "ER lumen"),
  "Golgi apparatus" = c("golgi", "golgi apparatus", "golgi membrane", "golgi lumen", "trans-golgi", "cis-golgi"),
  "Peroxisome" = c("peroxisome", "peroxisomal"),
  "Lysosome/Vacuole" = c("lysosome", "lysosomal", "vacuole", "vacuolar", "lytic vacuole"),
  "Cytoskeleton" = c("cytoskeleton", "microtubule", "actin", "intermediate filament", "microfilament", "tubulin", "spindle"),
  "Ribosome" = c("ribosome", "ribosomal", "ribonucleoprotein"),
  "Proteasome" = c("proteasome", "proteasomal"),
  "Plasma membrane" = c("plasma membrane", "cell surface", "cell junction", "focal adhesion"),
  "Extracellular" = c("extracellular", "extracellular space", "extracellular matrix", "extracellular region", "secreted"),
  "Chloroplast/Plastid" = c("chloroplast", "plastid", "thylakoid", "photosystem"),
  "Endosome" = c("endosome", "endosomal", "early endosome", "late endosome"),
  "Cilium/Flagellum" = c("cilium", "ciliary", "flagellum", "flagellar", "axoneme"),
  "Mitochondrial ribosome" = c("mitochondrial ribosome", "mitoribosome")
)

# Classify each GOCC term into organelle categories
classify_term <- function(term_name, term_id) {
  text <- tolower(paste(term_name, term_id))
  categories <- c()
  for (cat_name in names(organelle_categories)) {
    patterns <- organelle_categories[[cat_name]]
    if (any(sapply(patterns, function(p) grepl(p, text, fixed = TRUE)))) {
      categories <- c(categories, cat_name)
    }
  }
  if (length(categories) == 0) categories <- "Other"
  return(paste(categories, collapse = "|"))
}

enr[, organelle := mapply(classify_term, Description, ID)]
enr[, neg_log10_padj := -log10(p.adjust)]

# Explode pipe-separated organelle categories
enr_exploded <- enr[, .(organelle = unlist(strsplit(organelle, "\\|"))), 
                     by = .(ID, Description, species, direction, p.adjust, qvalue, neg_log10_padj, Count, GeneRatio, BgRatio, osd_datasets)]

# Aggregate by organelle category per species/direction
organelle_summary <- enr_exploded[, .(
  n_terms = .N,
  min_padj = min(p.adjust),
  combined_neg_log10_padj = sum(neg_log10_padj),
  n_genes = sum(as.integer(Count)),
  go_terms = paste(unique(Description), collapse = "; ")
), by = .(organelle, species, direction)]

# Also create a signed score: positive for up, negative for down
organelle_summary[, signed_score := ifelse(direction == "up", combined_neg_log10_padj, -combined_neg_log10_padj)]

# Save
fwrite(organelle_summary, '/mnt/shared-workspace/shared/results/enrichment/organelle_enrichment_summary.csv')
fwrite(enr_exploded, '/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_classified.csv')

cat("=== Organelle enrichment summary ===\n")
cat("Total organelle-species-direction entries:", nrow(organelle_summary), "\n\n")
print(organelle_summary[, .(organelle, species, direction, n_terms, min_padj, n_genes)][order(organelle, species, direction)])

cat("\n\n=== Organelles enriched in multiple species ===\n")
multi_species_organelle <- organelle_summary[, .(
  n_species = uniqueN(species),
  species_list = paste(unique(species), collapse = ", "),
  n_directions_up = sum(direction == "up"),
  n_directions_down = sum(direction == "down"),
  min_padj = min(min_padj)
), by = organelle][order(-n_species)]
print(multi_species_organelle)

# --- R code block 4 ---
suppressPackageStartupMessages({ library(data.table); library(GO.db) })

enr <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')
cat("Columns:", paste(colnames(enr), collapse=", "), "\n")

organelle_categories <- list(
  "Nucleus" = c("nucleus", "nuclear", "chromatin", "nucleolus", "nucleoplasm", "nuclear envelope", "nuclear pore"),
  "Mitochondrion" = c("mitochondr", "mitochondrial", "oxidative phosphorylation", "respiratory chain", "mitochondrial matrix", "mitochondrial membrane"),
  "Endoplasmic Reticulum" = c("endoplasmic reticulum", "endoplasmic", "er membrane", "er lumen"),
  "Golgi apparatus" = c("golgi"),
  "Peroxisome" = c("peroxisome", "peroxisomal"),
  "Lysosome/Vacuole" = c("lysosome", "lysosomal", "vacuole", "vacuolar"),
  "Cytoskeleton" = c("cytoskeleton", "microtubule", "actin", "intermediate filament", "microfilament", "tubulin", "spindle"),
  "Ribosome" = c("ribosome", "ribosomal", "ribonucleoprotein"),
  "Proteasome" = c("proteasome", "proteasomal"),
  "Plasma membrane" = c("plasma membrane", "cell surface", "cell junction", "focal adhesion"),
  "Extracellular" = c("extracellular"),
  "Chloroplast/Plastid" = c("chloroplast", "plastid", "thylakoid", "photosystem"),
  "Endosome" = c("endosome", "endosomal"),
  "Cilium/Flagellum" = c("cilium", "ciliary", "flagellum", "flagellar", "axoneme")
)

classify_term <- function(term_name) {
  text <- tolower(term_name)
  categories <- c()
  for (cat_name in names(organelle_categories)) {
    if (any(sapply(organelle_categories[[cat_name]], function(p) grepl(p, text, fixed = TRUE)))) {
      categories <- c(categories, cat_name)
    }
  }
  if (length(categories) == 0) categories <- "Other"
  return(paste(categories, collapse = "|"))
}

enr[, organelle := sapply(Description, classify_term)]
enr[, neg_log10_padj := -log10(p.adjust)]

# Explode
enr_exploded <- enr[, .(organelle = unlist(strsplit(organelle, "\\|"))), 
                     by = .(ID, Description, species, direction, p.adjust, qvalue, neg_log10_padj, Count)]

organelle_summary <- enr_exploded[, .(
  n_terms = .N,
  min_padj = min(p.adjust),
  combined_neg_log10_padj = sum(neg_log10_padj),
  n_genes = sum(as.integer(Count)),
  go_terms = paste(unique(Description), collapse = "; ")
), by = .(organelle, species, direction)]

organelle_summary[, signed_score := ifelse(direction == "up", combined_neg_log10_padj, -combined_neg_log10_padj)]

fwrite(organelle_summary, '/mnt/shared-workspace/shared/results/enrichment/organelle_enrichment_summary.csv')
fwrite(enr_exploded, '/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_classified.csv')

cat("\n=== Organelles enriched in multiple species ===\n")
multi <- organelle_summary[, .(
  n_species = uniqueN(species),
  species_list = paste(unique(species), collapse = ", "),
  n_up = sum(direction == "up"),
  n_down = sum(direction == "down"),
  min_padj = min(min_padj)
), by = organelle][order(-n_species)]
print(multi)

cat("\n=== Conserved organelle responses (>=3 species, same direction) ===\n")
conserved <- organelle_summary[, {
  n_sp = uniqueN(species)
  dir_consistent = (all(direction == "up") || all(direction == "down"))
  list(n_species = n_sp, direction_conserved = dir_consistent, 
       direction = if(all(direction=="up")) "up" else if(all(direction=="down")) "down" else "mixed",
       species = paste(unique(species), collapse=", "))
}, by = organelle][n_species >= 3][order(-n_species)]
print(conserved)
fwrite(conserved, '/mnt/shared-workspace/shared/results/meta/conserved_organelles.csv')

# --- R code block 5 ---
suppressPackageStartupMessages({ library(data.table) })

enr_exploded <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_classified.csv')

# Check "Other" terms
other_terms <- enr_exploded[organelle == "Other", .(Description, species, direction, p.adjust)]
cat("=== 'Other' GOCC terms ===\n")
print(other_terms[order(Description)])

# Refine: add more patterns
cat("\n=== All unique GOCC terms ===\n")
print(enr_exploded[, .(Description, species, direction)][order(Description)])

# --- R code block 6 ---
suppressPackageStartupMessages({ library(data.table) })

enr <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')

# Refined organelle categories with more patterns
organelle_categories <- list(
  "Nucleus" = c("nucleus", "nuclear", "chromatin", "nucleolus", "nucleoplasm", "nuclear envelope", "nuclear pore",
                "cenp", "synaptonemal", "cmg complex", "replication fork", "chromosome", "kinetochore", "condensin",
                "cohesin", "nucleolar", "nucleic acid", "dna repair", "replication complex"),
  "Mitochondrion" = c("mitochondr", "oxidative phosphorylation", "respiratory chain", "mitochondrial matrix", 
                       "mitochondrial membrane", "mitochondrial ribosome", "mitoribosome", "mrc", "complex i", 
                       "complex ii", "complex iii", "complex iv", "complex v", "atp synthase"),
  "Endoplasmic Reticulum" = c("endoplasmic reticulum", "endoplasmic", "er membrane", "er lumen", "er-golgi",
                               "rough endoplasmic", "smooth endoplasmic"),
  "Golgi apparatus" = c("golgi", "trans-golgi", "cis-golgi"),
  "Peroxisome" = c("peroxisome", "peroxisomal"),
  "Lysosome/Vacuole" = c("lysosome", "lysosomal", "vacuole", "vacuolar"),
  "Cytoskeleton" = c("cytoskeleton", "microtubule", "actin", "intermediate filament", "microfilament", 
                     "tubulin", "spindle", "a band", "i band", "sarcomere", "myofibril", "muscle fiber",
                     "stress fiber", "filamentous", "cortical cytoskeleton"),
  "Ribosome" = c("ribosome", "ribosomal", "ribonucleoprotein", "preribosome", "polysome", "translation"),
  "Proteasome" = c("proteasome", "proteasomal", "ubiquitin ligase complex", "ring ubiquitin", "cul2-ring",
                   "cullin-ring", "e3 ubiquitin", "proteasome complex"),
  "Plasma membrane" = c("plasma membrane", "cell surface", "cell junction", "focal adhesion", "cell-substrate"),
  "Extracellular" = c("extracellular", "secreted", "lipoprotein particle"),
  "Chloroplast/Plastid" = c("chloroplast", "plastid", "thylakoid", "photosystem", "plastoglobule", "chloroplast envelope"),
  "Endosome" = c("endosome", "endosomal"),
  "Cilium/Flagellum" = c("cilium", "ciliary", "flagellum", "flagellar", "axoneme", "basal body"),
  "Membrane transport" = c("transmembrane transporter", "transporter complex", "membrane transport", 
                            "ion channel", "solute carrier", "porin"),
  "Cell wall/Extracellular matrix" = c("cell wall", "casparian", "extracellular matrix", "collagen", "basement membrane"),
  "Lipid particle" = c("lipid particle", "lipid droplet", "oil body")
)

classify_term <- function(term_name) {
  text <- tolower(term_name)
  categories <- c()
  for (cat_name in names(organelle_categories)) {
    if (any(sapply(organelle_categories[[cat_name]], function(p) grepl(p, text, fixed = TRUE)))) {
      categories <- c(categories, cat_name)
    }
  }
  if (length(categories) == 0) categories <- "Other"
  return(paste(categories, collapse = "|"))
}

enr[, organelle := sapply(Description, classify_term)]
enr[, neg_log10_padj := -log10(p.adjust)]

enr_exploded <- enr[, .(organelle = unlist(strsplit(organelle, "\\|"))), 
                     by = .(ID, Description, species, direction, p.adjust, qvalue, neg_log10_padj, Count)]

organelle_summary <- enr_exploded[, .(
  n_terms = .N,
  min_padj = min(p.adjust),
  combined_neg_log10_padj = sum(neg_log10_padj),
  n_genes = sum(as.integer(Count)),
  go_terms = paste(unique(Description), collapse = "; ")
), by = .(organelle, species, direction)]

organelle_summary[, signed_score := ifelse(direction == "up", combined_neg_log10_padj, -combined_neg_log10_padj)]

fwrite(organelle_summary, '/mnt/shared-workspace/shared/results/enrichment/organelle_enrichment_summary.csv')
fwrite(enr_exploded, '/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_classified.csv')

cat("=== Organelles enriched in >=2 species ===\n")
multi <- organelle_summary[, .(
  n_species = uniqueN(species),
  species_list = paste(unique(species), collapse = ", "),
  n_up = sum(direction == "up"),
  n_down = sum(direction == "down"),
  min_padj = min(min_padj)
), by = organelle][n_species >= 2][order(-n_species)]
print(multi)

cat("\n=== Conserved organelle responses (>=3 species) ===\n")
conserved <- organelle_summary[, {
  n_sp = uniqueN(species)
  dirs = unique(direction)
  list(n_species = n_sp, 
       direction = if(length(dirs)==1) dirs else "mixed",
       species = paste(unique(species), collapse=", "),
       min_padj = min(min_padj))
}, by = organelle][n_species >= 3][order(-n_species)]
print(conserved)
fwrite(conserved, '/mnt/shared-workspace/shared/results/meta/conserved_organelles.csv')

cat("\n=== Remaining 'Other' terms ===\n")
print(enr_exploded[organelle == "Other", .(Description, species, direction)][order(Description)])

# --- R code block 7 ---
suppressPackageStartupMessages({ library(data.table) })

enr <- fread('/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_all.csv')

# Final refined categories
organelle_categories <- list(
  "Nucleus" = c("nucleus", "nuclear", "chromatin", "nucleolus", "nucleoplasm", "nuclear envelope", "nuclear pore",
                "cenp", "synaptonemal", "cmg complex", "replication fork", "chromosome", "kinetochore", "condensin",
                "cohesin", "nucleolar", "dna repair", "replication complex", "dna polymerase", "mcm complex",
                "replisome", "protein-dna complex", "nucleosome", "chromocenter", "chromosomal region",
                "anaphase-promoting", "cyclin-dependent", "xy body", "site of dna damage", "dna replication",
                "lateral element", "cell division site", "cleavage furrow", "midbody", "spindle",
                "germ plasm", "pole plasm", "p granule", "nucleoid"),
  "Mitochondrion" = c("mitochondr", "oxidative phosphorylation", "respiratory chain", "mitochondrial matrix", 
                       "mitochondrial membrane", "mitochondrial ribosome", "mitoribosome", "atp synthase",
                       "nadh dehydrogenase", "cytochrome complex", "oxidoreductase complex", "proton-transporting",
                       "alpha-ketoacid dehydrogenase"),
  "Endoplasmic Reticulum" = c("endoplasmic reticulum", "endoplasmic", "er membrane", "er lumen", "er body"),
  "Golgi apparatus" = c("golgi", "trans-golgi", "cis-golgi", "secretory vesicle"),
  "Peroxisome" = c("peroxisome", "peroxisomal", "microbody"),
  "Lysosome/Vacuole" = c("lysosome", "lysosomal", "vacuole", "vacuolar"),
  "Cytoskeleton" = c("cytoskeleton", "microtubule", "actin", "intermediate filament", "microfilament", 
                     "tubulin", "a band", "i band", "m band", "z disc", "sarcomere", "myofibril", "muscle fiber",
                     "stress fiber", "myofilament", "myosin", "kinesin", "microvillus", "brush border",
                     "supramolecular fiber", "supramolecular polymer", "keratin filament", "pseudopodium",
                     "cleavage furrow"),
  "Ribosome" = c("ribosome", "ribosomal", "ribonucleoprotein", "preribosome", "polysome", "u5 snrnp",
                 "risc complex", "spliceosomal"),
  "Proteasome" = c("proteasome", "proteasomal", "ubiquitin ligase", "ring ubiquitin", "cul2-ring",
                   "cullin-ring", "e3 ubiquitin", "proteasome complex", "endopeptidase clp"),
  "Plasma membrane" = c("plasma membrane", "cell surface", "cell junction", "focal adhesion", "cell-substrate",
                        "adherens junction", "desmosome", "cell projection membrane", "receptor complex",
                        "apical part", "basal part", "intercellular bridge"),
  "Extracellular" = c("extracellular", "secreted", "lipoprotein particle", "laminin", "cornified envelope",
                      "external encapsulating", "apoplast"),
  "Chloroplast/Plastid" = c("chloroplast", "plastid", "thylakoid", "photosystem", "plastoglobule", 
                             "chloroplast envelope", "photosynthetic membrane", "stromule", "nad(p)h dehydrogenase"),
  "Endosome" = c("endosome", "endosomal"),
  "Cilium/Flagellum" = c("cilium", "ciliary", "flagellum", "flagellar", "axoneme", "basal body"),
  "Membrane transport" = c("transmembrane transporter", "transporter complex", "membrane transport", 
                            "ion channel", "solute carrier", "porin", "extrinsic component of membrane"),
  "Cell wall/Extracellular matrix" = c("cell wall", "casparian", "extracellular matrix", "collagen", 
                                        "basement membrane", "supramolecular complex"),
  "Lipid particle" = c("lipid particle", "lipid droplet", "oil body", "protein-lipid complex")
)

classify_term <- function(term_name) {
  text <- tolower(term_name)
  categories <- c()
  for (cat_name in names(organelle_categories)) {
    if (any(sapply(organelle_categories[[cat_name]], function(p) grepl(p, text, fixed = TRUE)))) {
      categories <- c(categories, cat_name)
    }
  }
  if (length(categories) == 0) categories <- "Other"
  return(paste(categories, collapse = "|"))
}

enr[, organelle := sapply(Description, classify_term)]
enr[, neg_log10_padj := -log10(p.adjust)]

enr_exploded <- enr[, .(organelle = unlist(strsplit(organelle, "\\|"))), 
                     by = .(ID, Description, species, direction, p.adjust, qvalue, neg_log10_padj, Count)]

organelle_summary <- enr_exploded[, .(
  n_terms = .N,
  min_padj = min(p.adjust),
  combined_neg_log10_padj = sum(neg_log10_padj),
  n_genes = sum(as.integer(Count)),
  go_terms = paste(unique(Description), collapse = "; ")
), by = .(organelle, species, direction)]

organelle_summary[, signed_score := ifelse(direction == "up", combined_neg_log10_padj, -combined_neg_log10_padj)]

fwrite(organelle_summary, '/mnt/shared-workspace/shared/results/enrichment/organelle_enrichment_summary.csv')
fwrite(enr_exploded, '/mnt/shared-workspace/shared/results/enrichment/gocc_enrichment_classified.csv')

cat("=== Organelles enriched in >=2 species ===\n")
multi <- organelle_summary[, .(
  n_species = uniqueN(species),
  species_list = paste(unique(species), collapse=", "),
  n_up = sum(direction == "up"),
  n_down = sum(direction == "down"),
  min_padj = min(min_padj)
), by = organelle][n_species >= 2][order(-n_species)]
print(multi)

cat("\n=== Remaining 'Other' count ===\n")
cat(enr_exploded[organelle == "Other", .N], "terms\n")
print(enr_exploded[organelle == "Other", .(Description, species, direction)])

