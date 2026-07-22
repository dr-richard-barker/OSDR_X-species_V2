# ======================================================================
# Step 7: KEGG Pathway and Cofactor Mapping
# 
# Description: Fetch KEGG pathway-gene and cofactor-gene links via REST API; build log2FC tables
# 
# Inputs: KEGG REST API (rest.kegg.jp), results/deg/all_degs_with_orthologs.csv
# Outputs: results/meta/cofactor_gene_mapping.csv, results/meta/pathway_gene_log2fc.csv, results/kegg/*.csv
# 
# Language: R
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- R code block 1 ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggkegg)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
})

# First, identify top conserved pathways to visualize
# Based on our organelle enrichment results, the key conserved pathways are:
# 1. Oxidative phosphorylation (mitochondrion - conserved downregulation)
# 2. Proteasome (conserved)
# 3. Ribosome (conserved)
# 4. TCA cycle (mitochondrion)

# Load DEGs with orthologs to get per-species log2FC for pathway genes
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gene_id)]

# Get human gene symbols for pathway mapping
library(org.Hs.eg.db)
human_ens <- unique(deg[species == "Homo sapiens", gene_id])
human_sym <- AnnotationDbi::select(org.Hs.eg.db, keys = sub("\\.\\d+$", "", human_ens), 
                                    columns = c("SYMBOL","ENTREZID"), keytype = "ENSEMBL")
setnames(human_sym, "ENSEMBL", "human_ensembl")
deg[species == "Homo sapiens", human_ensembl := sub("\\.\\d+$", "", gene_id)]
deg <- merge(deg, human_sym[, .(human_ensembl, SYMBOL, ENTREZID)], by = "human_ensembl", all.x = TRUE)

# For each species, get the mean log2FC per human ortholog (across datasets)
ortho_fc <- deg[!is.na(human_ensembl) & !is.na(log2fc), 
                .(mean_log2fc = mean(log2fc, na.rm = TRUE),
                  n_datasets = uniqueN(osd),
                  n_deg = sum(is_deg),
                  direction = ifelse(mean(log2fc, na.rm = TRUE) > 0, "up", "down")),
                by = .(species, human_ensembl, SYMBOL, ENTREZID)]

# Get KEGG orthology mapping for key pathway genes
# Oxidative phosphorylation: hsa00190
# Proteasome: hsa03050
# Ribosome: hsa03010
# TCA cycle: hsa00020

# Map human ENTREZID to KEGG
library(KEGGREST)
# Get KEGG gene mapping for human
hsa_kegg <- keggLink("hsa", "pathway")
hsa_kegg_df <- data.table(
  kegg_gene = names(hsa_kegg),
  pathway = hsa_kegg
)
hsa_kegg_df[, entrez := sub("hsa:", "", kegg_gene)]
hsa_kegg_df[, pathway_id := sub("path:hsa", "", pathway)]

cat("KEGG gene-pathway mappings:", nrow(hsa_kegg_df), "\n")
cat("Pathways available:\n")
print(hsa_kegg_df[, .N, by = pathway_id][order(-N)][1:20])

# Key pathways to visualize
key_pathways <- list(
  "Oxidative phosphorylation" = "00190",
  "Proteasome" = "03050",
  "Ribosome" = "03010",
  "TCA cycle" = "00020",
  "Glutathione metabolism" = "00480"
)

# For each pathway, get the genes and their per-species log2FC
pathway_gene_fc <- list()
for (pname in names(key_pathways)) {
  pid <- key_pathways[[pname]]
  pathway_genes <- hsa_kegg_df[pathway_id == pid, entrez]
  pathway_fc <- ortho_fc[ENTREZID %in% pathway_genes]
  pathway_fc[, pathway := pname]
  pathway_fc[, pathway_id := pid]
  pathway_gene_fc[[pname]] <- pathway_fc
  cat(sprintf("\n%s (hsa%s): %d genes, %d with ortholog FC data\n", 
              pname, pid, length(pathway_genes), nrow(pathway_fc)))
  cat("  Species coverage:", paste(unique(pathway_fc$species), collapse=", "), "\n")
}

pathway_fc_all <- rbindlist(pathway_gene_fc, fill = TRUE)
fwrite(pathway_fc_all, '/mnt/shared-workspace/shared/results/meta/pathway_gene_log2fc.csv')
cat("\nSaved pathway gene log2FC table\n")

# --- R code block 2 ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggkegg)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(KEGGREST)
  library(org.Hs.eg.db)
})

# Load DEGs with orthologs
deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gid_id := gene_id)]

# Get human gene symbols
human_ens <- unique(deg[species == "Homo sapiens", sub("\\.\\d+$", "", gene_id)])
human_sym <- AnnotationDbi::select(org.Hs.eg.db, keys = human_ens, 
                                    columns = c("SYMBOL","ENTREZID"), keytype = "ENSEMBL")
setnames(human_sym, "ENSEMBL", "human_ensembl")
deg[species == "Homo sapiens", human_ensembl := sub("\\.\\d+$", "", gene_id)]
deg <- merge(deg, human_sym[, .(human_ensembl, SYMBOL, ENTREZID)], by = "human_ensembl", all.x = TRUE)

# Per-species mean log2FC per human ortholog
ortho_fc <- deg[!is.na(human_ensembl) & !is.na(log2fc), 
                .(mean_log2fc = mean(log2fc, na.rm = TRUE),
                  n_datasets = uniqueN(osd),
                  n_deg = sum(is_deg)),
                by = .(species, human_ensembl, SYMBOL, ENTREZID)]

# Get KEGG gene-pathway mapping
hsa_kegg <- keggLink("hsa", "pathway")
hsa_kegg_df <- data.table(kegg_gene = names(hsa_kegg), pathway = hsa_kegg)
hsa_kegg_df[, entrez := sub("hsa:", "", kegg_gene)]
hsa_kegg_df[, pathway_id := sub("path:hsa", "", pathway)]

cat("KEGG mappings:", nrow(hsa_kegg_df), "\n")

# Key pathways
key_pathways <- list(
  "Oxidative phosphorylation" = "00190",
  "Proteasome" = "03050",
  "Ribosome" = "03010",
  "TCA cycle" = "00020"
)

pathway_gene_fc <- list()
for (pname in names(key_pathways)) {
  pid <- key_pathways[[pname]]
  pathway_genes <- hsa_kegg_df[pathway_id == pid, entrez]
  pathway_fc <- ortho_fc[ENTREZID %in% pathway_genes]
  pathway_fc[, pathway := pname]
  pathway_fc[, pathway_id := pid]
  pathway_gene_fc[[pname]] <- pathway_fc
  cat(sprintf("%s (hsa%s): %d KEGG genes, %d with FC data\n", pname, pid, length(pathway_genes), nrow(pathway_fc)))
}

pathway_fc_all <- rbindlist(pathway_gene_fc, fill = TRUE)
fwrite(pathway_fc_all, '/mnt/shared-workspace/shared/results/meta/pathway_gene_log2fc.csv')

# Now build cofactor mapping for pathway enzymes
# Key cofactors: NAD+, FAD, CoA, Fe-S, PLP, TPP, biotin, lipoate, heme, Zn2+
# Map via KEGG Compound/Reaction
cofactors <- data.table(
  cofactor = c("NAD+", "NADH", "FAD", "FADH2", "CoA", "Fe-S cluster", "PLP (pyridoxal phosphate)",
               "TPP (thiamine pyrophosphate)", "Biotin", "Lipoamide", "Heme", "Zn2+", "FMN", "NADP+", "NADPH"),
  kegg_compound = c("C00003", "C00004", "C00016", "C01352", "C00010", "C00824", "C00018",
                    "C00068", "C00120", "C00248", "C00032", "C00038", "C00061", "C00006", "C00005"),
  category = c("Electron transfer", "Electron transfer", "Electron transfer", "Electron transfer",
               "Acyl transfer", "Electron transfer", "Amino transfer", "Decarboxylation",
               "Carboxylation", "Acyl transfer", "Electron transfer", "Structural/catalytic",
               "Electron transfer", "Electron transfer", "Electron transfer")
)

# Get KEGG compound-gene (enzyme) relationships for cofactors
cofactor_genes <- list()
for (i in 1:nrow(cofactors)) {
  comp <- cofactors$kegg_compound[i]
  cf <- cofactors$cofactor[i]
  tryCatch({
    # Get reactions involving this compound
    reactions <- keggLink("reaction", paste0("cpd:", comp))
    if (length(reactions) > 0) {
      # Get enzymes for these reactions
      enzymes <- keggLink("enzyme", "reaction")
      # This is complex - let's use a simpler approach
      # Get genes that use this compound
      comp_genes <- keggLink("hsa", paste0("cpd:", comp))
      if (length(comp_genes) > 0) {
        gene_entrez <- sub("hsa:", "", comp_genes)
        cofactor_genes[[cf]] <- data.table(cofactor = cf, kegg_compound = comp, 
                                            entrez = gene_entrez, category = cofactors$category[i])
      }
    }
  }, error = function(e) {})
  cat(sprintf("  %s: %d genes\n", cf, length(cofactor_genes[[cf]])))
}

cofactor_gene_dt <- rbindlist(cofactor_genes, fill = TRUE)
fwrite(cofactor_gene_dt, '/mnt/shared-workspace/shared/results/meta/cofactor_gene_mapping.csv')
cat("\nTotal cofactor-gene mappings:", nrow(cofactor_gene_dt), "\n")

# Merge cofactor mapping with pathway gene FC
cofactor_pathway_fc <- merge(cofactor_gene_dt, pathway_fc_all, by = "entrez", all.x = TRUE, allow.cartesian = TRUE)
fwrite(cofactor_pathway_fc, '/mnt/shared-workspace/shared/results/meta/cofactor_pathway_log2fc.csv')
cat("Cofactor-pathway-FC rows:", nrow(cofactor_pathway_fc), "\n")

# --- R code block 3 ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggkegg)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(KEGGREST)
  library(org.Hs.eg.db)
})

deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gene_id)]

# Get human gene symbols
human_ens <- unique(deg[species == "Homo sapiens", sub("\\.\\d+$", "", gene_id)])
human_sym <- AnnotationDbi::select(org.Hs.eg.db, keys = human_ens, 
                                    columns = c("SYMBOL","ENTREZID"), keytype = "ENSEMBL")
setnames(human_sym, "ENSEMBL", "human_ensembl")
deg[species == "Homo sapiens", human_ensembl := sub("\\.\\d+$", "", gene_id)]
deg <- merge(deg, human_sym[, .(human_ensembl, SYMBOL, ENTREZID)], by = "human_ensembl", all.x = TRUE)

# Per-species mean log2FC per human ortholog
ortho_fc <- deg[!is.na(human_ensembl) & !is.na(log2fc), 
                .(mean_log2fc = mean(log2fc, na.rm = TRUE),
                  n_datasets = uniqueN(osd),
                  n_deg = sum(is_deg)),
                by = .(species, human_ensembl, SYMBOL, ENTREZID)]

# KEGG gene-pathway mapping
hsa_kegg <- keggLink("hsa", "pathway")
hsa_kegg_df <- data.table(kegg_gene = names(hsa_kegg), pathway = hsa_kegg)
hsa_kegg_df[, entrez := sub("hsa:", "", kegg_gene)]
hsa_kegg_df[, pathway_id := sub("path:hsa", "", pathway)]

key_pathways <- list(
  "Oxidative phosphorylation" = "00190",
  "Proteasome" = "03050",
  "Ribosome" = "03010",
  "TCA cycle" = "00020"
)

pathway_gene_fc <- list()
for (pname in names(key_pathways)) {
  pid <- key_pathways[[pname]]
  pathway_genes <- hsa_kegg_df[pathway_id == pid, entrez]
  pathway_fc <- ortho_fc[ENTREZID %in% pathway_genes]
  pathway_fc[, pathway := pname]
  pathway_fc[, pathway_id := pid]
  pathway_gene_fc[[pname]] <- pathway_fc
  cat(sprintf("%s (hsa%s): %d KEGG genes, %d with FC data\n", pname, pid, length(pathway_genes), nrow(pathway_fc)))
}

pathway_fc_all <- rbindlist(pathway_gene_fc, fill = TRUE)
fwrite(pathway_fc_all, '/mnt/shared-workspace/shared/results/meta/pathway_gene_log2fc.csv')

# Cofactor mapping
cofactors <- data.table(
  cofactor = c("NAD+", "NADH", "FAD", "FADH2", "CoA", "Fe-S cluster", "PLP",
               "TPP", "Biotin", "Lipoamide", "Heme", "Zn2+", "FMN", "NADP+", "NADPH"),
  kegg_compound = c("C00003", "C00004", "C00016", "C01352", "C00010", "C00824", "C00018",
                    "C00068", "C00120", "C00248", "C00032", "C00038", "C00061", "C00006", "C00005"),
  category = c("Electron transfer", "Electron transfer", "Electron transfer", "Electron transfer",
               "Acyl transfer", "Electron transfer", "Amino transfer", "Decarboxylation",
               "Carboxylation", "Acyl transfer", "Electron transfer", "Structural/catalytic",
               "Electron transfer", "Electron transfer", "Electron transfer")
)

cofactor_genes <- list()
for (i in 1:nrow(cofactors)) {
  comp <- cofactors$kegg_compound[i]
  cf <- cofactors$cofactor[i]
  tryCatch({
    comp_genes <- keggLink("hsa", paste0("cpd:", comp))
    if (length(comp_genes) > 0) {
      gene_entrez <- sub("hsa:", "", comp_genes)
      cofactor_genes[[cf]] <- data.table(cofactor = cf, kegg_compound = comp, 
                                          entrez = gene_entrez, category = cofactors$category[i])
    }
  }, error = function(e) {})
  cat(sprintf("  %s: %d genes\n", cf, ifelse(is.null(cofactor_genes[[cf]]), 0, nrow(cofactor_genes[[cf]]))))
}

cofactor_gene_dt <- rbindlist(cofactor_genes, fill = TRUE)
fwrite(cofactor_gene_dt, '/mnt/shared-workspace/shared/results/meta/cofactor_gene_mapping.csv')
cat("\nTotal cofactor-gene mappings:", nrow(cofactor_gene_dt), "\n")

# Merge cofactor with pathway FC
cofactor_pathway_fc <- merge(cofactor_gene_dt, pathway_fc_all, by = "entrez", all.x = TRUE, allow.cartesian = TRUE)
fwrite(cofactor_pathway_fc, '/mnt/shared-workspace/shared/results/meta/cofactor_pathway_log2fc.csv')
cat("Cofactor-pathway-FC rows:", nrow(cofactor_pathway_fc), "\n")

# --- R code block 4 ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggkegg)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(KEGGREST)
  library(org.Hs.eg.db)
})

deg <- fread('/mnt/shared-workspace/shared/results/deg/all_degs_with_orthologs.csv')
deg[, gene_id := as.character(gene_id)]

# Get human gene symbols
human_ens <- unique(deg[species == "Homo sapiens", sub("\\.\\d+$", "", gene_id)])
human_sym <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = human_ens, 
                                    columns = c("SYMBOL","ENTREZID"), keytype = "ENSEMBL"))
setnames(human_sym, "ENSEMBL", "human_ensembl")
deg[species == "Homo sapiens", human_ensembl := sub("\\.\\d+$", "", gene_id)]
deg <- merge(deg, human_sym[, .(human_ensembl, SYMBOL, ENTREZID)], by = "human_ensembl", all.x = TRUE)

# Per-species mean log2FC per human ortholog
ortho_fc <- deg[!is.na(human_ensembl) & !is.na(log2fc), 
                .(mean_log2fc = mean(log2fc, na.rm = TRUE),
                  n_datasets = uniqueN(osd),
                  n_deg = sum(is_deg)),
                by = .(species, human_ensembl, SYMBOL, ENTREZID)]

# KEGG gene-pathway mapping
hsa_kegg <- keggLink("hsa", "pathway")
hsa_kegg_df <- data.table(kegg_gene = names(hsa_kegg), pathway = hsa_kegg)
hsa_kegg_df[, entrez := sub("hsa:", "", kegg_gene)]
hsa_kegg_df[, pathway_id := sub("path:hsa", "", pathway)]

key_pathways <- list(
  "Oxidative phosphorylation" = "00190",
  "Proteasome" = "03050",
  "Ribosome" = "03010",
  "TCA cycle" = "00020"
)

pathway_gene_fc <- list()
for (pname in names(key_pathways)) {
  pid <- key_pathways[[pname]]
  pathway_genes <- hsa_kegg_df[pathway_id == pid, entrez]
  pathway_fc <- ortho_fc[ENTREZID %in% pathway_genes]
  pathway_fc[, pathway := pname]
  pathway_fc[, pathway_id := pid]
  pathway_gene_fc[[pname]] <- pathway_fc
  cat(sprintf("%s (hsa%s): %d KEGG genes, %d with FC data\n", pname, pid, length(pathway_genes), nrow(pathway_fc)))
}

pathway_fc_all <- rbindlist(pathway_gene_fc, fill = TRUE)
fwrite(pathway_fc_all, '/mnt/shared-workspace/shared/results/meta/pathway_gene_log2fc.csv')

# Cofactor mapping
cofactors <- data.table(
  cofactor = c("NAD+", "NADH", "FAD", "FADH2", "CoA", "Fe-S cluster", "PLP",
               "TPP", "Biotin", "Lipoamide", "Heme", "Zn2+", "FMN", "NADP+", "NADPH"),
  kegg_compound = c("C00003", "C00004", "C00016", "C01352", "C00010", "C00824", "C00018",
                    "C00068", "C00120", "C00248", "C00032", "C00038", "C00061", "C00006", "C00005"),
  category = c("Electron transfer", "Electron transfer", "Electron transfer", "Electron transfer",
               "Acyl transfer", "Electron transfer", "Amino transfer", "Decarboxylation",
               "Carboxylation", "Acyl transfer", "Electron transfer", "Structural/catalytic",
               "Electron transfer", "Electron transfer", "Electron transfer")
)

cofactor_genes <- list()
for (i in 1:nrow(cofactors)) {
  comp <- cofactors$kegg_compound[i]
  cf <- cofactors$cofactor[i]
  tryCatch({
    comp_genes <- keggLink("hsa", paste0("cpd:", comp))
    if (length(comp_genes) > 0) {
      gene_entrez <- sub("hsa:", "", comp_genes)
      cofactor_genes[[cf]] <- data.table(cofactor = cf, kegg_compound = comp, 
                                          entrez = gene_entrez, category = cofactors$category[i])
    }
  }, error = function(e) {})
  n <- ifelse(is.null(cofactor_genes[[cf]]), 0, nrow(cofactor_genes[[cf]]))
  cat(sprintf("  %s: %d genes\n", cf, n))
}

cofactor_gene_dt <- rbindlist(cofactor_genes, fill = TRUE)
fwrite(cofactor_gene_dt, '/mnt/shared-workspace/shared/results/meta/cofactor_gene_mapping.csv')
cat("\nTotal cofactor-gene mappings:", nrow(cofactor_gene_dt), "\n")

# Merge cofactor with pathway FC
cofactor_pathway_fc <- merge(cofactor_gene_dt, pathway_fc_all, by = "entrez", all.x = TRUE, allow.cartesian = TRUE)
fwrite(cofactor_pathway_fc, '/mnt/shared-workspace/shared/results/meta/cofactor_pathway_log2fc.csv')
cat("Cofactor-pathway-FC rows:", nrow(cofactor_pathway_fc), "\n")
cat("Done\n")

# --- R code block 5 ---

suppressMessages({
  library(org.Hs.eg.db); library(data.table)
})
keys <- keys(org.Hs.eg.db, keytype="ENSEMBL")
mp <- AnnotationDbi::select(org.Hs.eg.db, keys=keys, columns=c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
setDT(mp)
mp <- mp[!is.na(ENTREZID)]
fwrite(mp, "/workspace/kegg_cache/human_ensembl_entrez_symbol.csv")
cat("Rows:", nrow(mp), "
")
print(head(mp))


