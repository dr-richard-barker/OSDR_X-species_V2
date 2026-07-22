# ======================================================================
# Step 2: DEG Extraction and Differential Expression
# 
# Description: Extract DEGs from GeneLab RCP pre-computed tables; run DESeq2 on OSD-96 and limma on OSD-35
# 
# Inputs: Processed DE tables and count matrices from Step 1
# Outputs: results/deg/all_degs_combined.csv, per-dataset DEG tables
# 
# Language: R
# Extracted from analysis notebook (worker-0.ipynb)
# See METHODS.md for full parameter details.
# ======================================================================

# --- R code block 1 ---
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

# === OSD-96: DESeq2 on Drosophila RSEM counts ===
counts96 <- read.csv('/mnt/shared-workspace/shared/osdr_data/OSD-96_counts.csv', row.names=1, check.names=FALSE)
meta96 <- read.csv('/mnt/shared-workspace/shared/osdr_data/OSD-96_metadata.csv')

# Ensure column order matches
counts96 <- counts96[, meta96$sample]
# Round counts to integers for DESeq2
counts96 <- round(as.matrix(counts96))
mode(counts96) <- "integer"

# Filter low-count genes
keep <- rowSums(counts96 >= 10) >= 3
counts96 <- counts96[keep,]
cat("OSD-96: genes after filtering:", nrow(counts96), "\n")

# DESeq2 with condition + timepoint + generation as covariates
meta96$condition <- factor(meta96$condition, levels=c("Ground_Control","Space_Flight"))
meta96$timepoint <- factor(meta96$timepoint)
meta96$generation <- factor(meta96$generation)

dds <- DESeqDataSetFromMatrix(counts96, meta96, design=~ generation + timepoint + condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Space_Flight","Ground_Control"), alpha=0.05)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$is_deg <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
res_df$deg_direction <- ifelse(res_df$is_deg & res_df$log2FoldChange > 0, "up",
                        ifelse(res_df$is_deg & res_df$log2FoldChange < 0, "down", "ns"))
res_df$osd <- "OSD-96"
res_df$species <- "Drosophila melanogaster"

n_deg <- sum(res_df$is_deg, na.rm=TRUE)
n_up <- sum(res_df$deg_direction=="up", na.rm=TRUE)
n_down <- sum(res_df$deg_direction=="down", na.rm=TRUE)
cat(sprintf("OSD-96 (Drosophila): %d genes, %d DEGs (%d up, %d down)\n", nrow(res_df), n_deg, n_up, n_down))

out96 <- res_df[,c("gene_id","log2FoldChange","pvalue","padj","is_deg","deg_direction","osd","species")]
colnames(out96)[2] <- "log2fc"
out96$symbol <- NA
out96 <- out96[,c("gene_id","symbol","log2fc","pvalue","padj","is_deg","deg_direction","osd","species")]
write.csv(out96, '/mnt/shared-workspace/shared/results/deg/OSD-96_deg.csv', row.names=FALSE)

# === OSD-35: limma on C. elegans microarray normalized expression ===
expr35 <- read.csv('/mnt/shared-workspace/shared/osdr_data/processed/OSD-35_counts_or_expr_GLDS-35_array_normalized_expression_probeset_GLmicroarray.csv', check.names=FALSE)
cat("\nOSD-35 columns:", paste(colnames(expr35)[1:15], collapse=", "), "\n")

# Sample columns are after annotation columns
annot_cols <- c("ENSEMBL","SYMBOL","GENENAME","REFSEQ","ENTREZID","STRING_id","GOSLIM_IDS","ProbesetID","count_ENSEMBL_mappings")
sample_cols <- setdiff(colnames(expr35), annot_cols)
cat("OSD-35 sample columns:", paste(sample_cols, collapse=", "), "\n")

# Build metadata from sample names
# Cele_N2_wo_GC, Cele_N2_wo_FLT, Cele_N2_wo_HG, Cele_N2_wo_VIV, Cele_N2_wo_suG
meta35 <- data.frame(
  sample = sample_cols,
  condition = ifelse(grepl("FLT", sample_cols), "Space_Flight",
              ifelse(grepl("GC", sample_cols) | grepl("suG", sample_cols), "Ground_Control",
              ifelse(grepl("VIV", sample_cols), "Vivarium_Control",
              ifelse(grepl("HG", sample_cols), "Hypergravity", "unknown"))))
)
print(meta35)

# Filter to FLT vs GC only
flt_gc_samples <- sample_cols[meta35$condition %in% c("Space_Flight","Ground_Control")]
expr_matrix <- as.matrix(expr35[, flt_gc_samples])
rownames(expr_matrix) <- expr35$ENSEMBL
mode(expr_matrix) <- "numeric"

# Remove rows with NA gene IDs
valid <- !is.na(rownames(expr_matrix)) & rownames(expr_matrix) != ""
expr_matrix <- expr_matrix[valid,]
# Take first if duplicated
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)),]
cat("OSD-35: genes after filtering:", nrow(expr_matrix), "\n")

# limma
condition35 <- factor(meta35$condition[meta35$sample %in% flt_gc_samples], levels=c("Ground_Control","Space_Flight"))
design <- model.matrix(~ condition35)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
res35 <- topTable(fit, coef="condition35Space_Flight", number=Inf, sort.by="none")
res35$gene_id <- rownames(res35)
res35$is_deg <- res35$adj.P.Val < 0.05 & abs(res35$logFC) > 1
res35$deg_direction <- ifelse(res35$is_deg & res35$logFC > 0, "up",
                       ifelse(res35$is_deg & res35$logFC < 0, "down", "ns"))
res35$osd <- "OSD-35"
res35$species <- "Caenorhabditis elegans"

n_deg35 <- sum(res35$is_deg)
cat(sprintf("OSD-35 (C. elegans): %d genes, %d DEGs (%d up, %d down)\n", nrow(res35), n_deg35,
            sum(res35$deg_direction=="up"), sum(res35$deg_direction=="down")))

# Add symbol mapping
sym_map <- setNames(expr35$SYMBOL[match(rownames(res35), expr35$ENSEMBL)], rownames(res35))
res35$symbol <- sym_map[rownames(res35)]

out35 <- res35[,c("gene_id","symbol","logFC","P.Value","adj.P.Val","is_deg","deg_direction","osd","species")]
colnames(out35)[3:5] <- c("log2fc","pvalue","padj")
write.csv(out35, '/mnt/shared-workspace/shared/results/deg/OSD-35_deg.csv', row.names=FALSE)

cat("\nDone with OSD-96 and OSD-35\n")

