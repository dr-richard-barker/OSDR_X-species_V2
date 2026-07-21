#!/usr/bin/env Rscript
# Figures 1, 2, 3, 4, 7 for the cross-species spaceflight transcriptomics manuscript
suppressMessages({
  library(ggplot2); library(dplyr); library(data.table); library(viridis)
  library(patchwork); library(ggrepel); library(UpSetR)
  library(ComplexHeatmap); library(circlize); library(RColorBrewer)
})

FIG_DIR <- "/mnt/shared-workspace/shared/results/figures"
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)

# Global theme
theme_set(theme_classic(base_family="Liberation Sans", base_size=9) +
          theme(panel.grid.minor=element_blank()))

# Species colors (viridis-based, colorblind-friendly, fixed globally)
SPECIES_COLORS <- c(
  "Homo sapiens"            = "#440154FF",
  "Mus musculus"            = "#3B528BFF",
  "Drosophila melanogaster" = "#21918CFF",
  "Caenorhabditis elegans"  = "#5EC962FF",
  "Saccharomyces cerevisiae"= "#FDE725FF",
  "Arabidopsis thaliana"    = "#9B1946FF"
)

# ============================================================
# Fig 1: Study overview - dataset counts per species + DEG counts
# ============================================================
cat("Building Fig 1: Study overview...\n")
degs <- fread("/mnt/shared-workspace/shared/results/deg/all_degs_combined.csv")
contrast_sel <- fread("/mnt/shared-workspace/shared/osdr_data/contrast_selection.csv")

# Per-species summary
species_summary <- degs[is_deg==TRUE, .(n_deg=.N, n_up=sum(deg_direction=="up"),
                       n_down=sum(deg_direction=="down")), by=species]
dataset_summary <- contrast_sel[, .(n_datasets=.N), by=species]
fig1_data <- merge(species_summary, dataset_summary, by="species", all=TRUE)
fig1_data[is.na(n_datasets), n_datasets := 0]
setorder(fig1_data, species)

p1a <- ggplot(fig1_data, aes(x=species, y=n_datasets, fill=species)) +
  geom_col(width=0.7) +
  geom_text(aes(label=n_datasets), vjust=-0.3, size=2.5, family="Liberation Sans") +
  scale_fill_manual(values=SPECIES_COLORS) +
  labs(title="A. Datasets per species", x=NULL, y="Number of datasets") +
  theme(axis.text.x=element_text(angle=30, hjust=1, face="italic"),
        legend.position="none")

# DEG counts - up vs down
fig1_deg_long <- melt(fig1_data[, .(species, n_up, n_down)],
                      id.vars="species", variable.name="direction", value.name="count")
fig1_deg_long[, direction := ifelse(direction=="n_up", "Upregulated", "Downregulated")]
p1b <- ggplot(fig1_deg_long, aes(x=species, y=count, fill=species, alpha=direction)) +
  geom_col(position="dodge", width=0.7) +
  scale_fill_manual(values=SPECIES_COLORS) +
  scale_alpha_manual(values=c(Upregulated=1, Downregulated=0.5)) +
  labs(title="B. DEGs per species", x=NULL, y="Number of DEGs", alpha=NULL) +
  theme(axis.text.x=element_text(angle=30, hjust=1, face="italic"),
        legend.position="right")

fig1 <- p1a / p1b + plot_layout(guides="collect")
ggsave(file.path(FIG_DIR, "fig1_study_overview.svg"), fig1, width=7, height=6, units="in")
ggsave(file.path(FIG_DIR, "fig1_study_overview.png"), fig1, width=7, height=6, units="in", dpi=300)
cat("  Saved fig1_study_overview\n")

# ============================================================
# Fig 2: Volcano plots (one panel per species, representative dataset)
# ============================================================
cat("Building Fig 2: Volcano plots...\n")
# Pick representative dataset per species (most DEGs)
rep_dataset <- degs[is_deg==TRUE, .(n_deg=.N), by=.(species, osd)][order(-n_deg)][, .SD[1], by=species]
cat("  Representative datasets:\n"); print(rep_dataset)

volcano_data <- merge(degs, rep_dataset[, .(species, osd)], by=c("species","osd"))
volcano_data[, neg_log10_p := -log10(pvalue)]
volcano_data[is_deg==TRUE, sig := ifelse(log2fc>0, "Up", "Down")]
volcano_data[is_deg==FALSE, sig := "NS"]

p2 <- ggplot(volcano_data, aes(x=log2fc, y=neg_log10_p, color=sig)) +
  geom_point(size=0.3, alpha=0.5) +
  scale_color_manual(values=c(Up="#B2182B", Down="#2166AC", NS="grey80")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50", linewidth=0.3) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50", linewidth=0.3) +
  facet_wrap(~species, scales="free", ncol=3, labeller=labeller(species=label_both)) +
  labs(x="log2 fold change (spaceflight vs ground)", y="-log10(p-value)",
       title="Volcano plots (representative dataset per species)") +
  theme(strip.background=element_rect(fill="grey95", color=NA),
        strip.text=element_text(face="italic", size=7),
        legend.position="bottom",
        legend.title=element_blank()) +
  guides(color=guide_legend(override.aes=list(size=2)))

ggsave(file.path(FIG_DIR, "fig2_volcano_plots.svg"), p2, width=9, height=6, units="in")
ggsave(file.path(FIG_DIR, "fig2_volcano_plots.png"), p2, width=9, height=6, units="in", dpi=300)
cat("  Saved fig2_volcano_plots\n")

# ============================================================
# Fig 3: Orthology UpSet plot
# ============================================================
cat("Building Fig 3: Orthology UpSet...\n")
ortho_wide <- fread("/mnt/shared-workspace/shared/orthodb/orthology_matrix_wide.csv")
# Human is anchor (column 1 = human_ensembl, always present)
# Non-human species columns
nonhuman_cols <- c("Mus musculus","Drosophila melanogaster","Caenorhabditis elegans",
                   "Saccharomyces cerevisiae","Arabidopsis thaliana")
all_cols <- c("Homo sapiens", nonhuman_cols)

# Build binary matrix (human always 1 since it's the anchor)
ortho_bin <- as.data.frame(ortho_wide[, ..nonhuman_cols])
ortho_bin[is.na(ortho_bin)] <- ""
ortho_bin[ortho_bin != ""] <- 1
ortho_bin[ortho_bin == ""] <- 0
for (col in nonhuman_cols) ortho_bin[[col]] <- as.integer(ortho_bin[[col]])
ortho_bin[["Homo sapiens"]] <- 1L  # anchor species

svg(file.path(FIG_DIR, "fig3_orthology_upset.svg"), width=8, height=5, family="Liberation Sans")
print(
upset(ortho_bin, sets=all_cols, order.by="freq", nsets=6,
      sets.bar.color=unname(SPECIES_COLORS[all_cols]),
      main.bar.color="grey40", matrix.color="grey20",
      point.size=2.5, line.size=0.8,
      text.scale=c(1.2, 1, 1, 0.9, 1.1, 0.9),
      mb.ratio=c(0.65, 0.35),
      mainbar.y.label="Ortholog group size", sets.x.label="Genes with ortholog")
)
dev.off()
png(file.path(FIG_DIR, "fig3_orthology_upset.png"), width=8, height=5, units="in", res=300, family="Liberation Sans")
print(
upset(ortho_bin, sets=all_cols, order.by="freq", nsets=6,
      sets.bar.color=unname(SPECIES_COLORS[all_cols]),
      main.bar.color="grey40", matrix.color="grey20",
      point.size=2.5, line.size=0.8,
      text.scale=c(1.2, 1, 1, 0.9, 1.1, 0.9),
      mb.ratio=c(0.65, 0.35),
      mainbar.y.label="Ortholog group size", sets.x.label="Genes with ortholog")
)
dev.off()
cat("  Saved fig3_orthology_upset\n")

# ============================================================
# Fig 4: Organelle enrichment heatmap (signed -log10 padj)
# ============================================================
cat("Building Fig 4: Organelle enrichment heatmap...\n")
enrich <- fread("/mnt/shared-workspace/shared/results/enrichment/organelle_enrichment_summary.csv")
# Pivot to organelle x species matrix using signed score
fig4_mat <- dcast(enrich, organelle ~ species, value.var="signed_score",
                  fun.aggregate=sum, fill=0)
mat <- as.matrix(fig4_mat[, -1])
rownames(mat) <- fig4_mat$organelle
# Order organelles by total abs score
ord <- order(rowSums(abs(mat)), decreasing=TRUE)
mat <- mat[ord, ]

col_fun <- colorRamp2(c(-50, 0, 50), c("#2166AC", "white", "#B2182B"))

svg(file.path(FIG_DIR, "fig4_organelle_enrichment_heatmap.svg"), width=7, height=6, family="Liberation Sans")
h <- Heatmap(mat, name="Signed\n-log10(padj)", col=col_fun,
             cluster_rows=FALSE, cluster_columns=FALSE,
             row_names_side="left", row_names_gp=gpar(fontsize=8, fontfamily="Liberation Sans"),
             column_names_gp=gpar(fontsize=8, fontfamily="Liberation Sans", fontface="italic"),
             column_title="Subcellular organelle enrichment across species",
             column_title_gp=gpar(fontsize=10, fontfamily="Liberation Sans", fontface="bold"),
             cell_fun=function(j,i,x,y,width,height,fill) {
               if (abs(mat[i,j]) > 5) {
                 grid.text(sprintf("%.0f", mat[i,j]), x, y, gp=gpar(fontsize=6, fontfamily="Liberation Sans",
                           col=ifelse(mat[i,j]>0, "white", "white")))
               }
             })
draw(h)
dev.off()

png(file.path(FIG_DIR, "fig4_organelle_enrichment_heatmap.png"), width=7, height=6, units="in", res=300, family="Liberation Sans")
draw(h)
dev.off()
cat("  Saved fig4_organelle_enrichment_heatmap\n")

# ============================================================
# Fig 7: Cofactor summary - cofactor x species mean FC for affected genes
# ============================================================
cat("Building Fig 7: Cofactor summary...\n")
cof_summary <- fread("/mnt/shared-workspace/shared/results/kegg/cofactor_pathway_summary.csv")

# Heatmap: cofactor x species, value = mean FC of genes using that cofactor
fig7_mat <- dcast(cof_summary, cofactor_name ~ species, value.var="mean_fc",
                  fun.aggregate=mean, fill=0, na.rm=TRUE)
mat7 <- as.matrix(fig7_mat[, -1])
rownames(mat7) <- fig7_mat$cofactor_name
# Order cofactors
cof_order <- c("NAD+","NADH","NADP+","NADPH","FAD","FADH2","FMN","CoA",
               "TPP (thiamine PP)","Lipoamide","Biotin","Heme")
mat7 <- mat7[intersect(cof_order, rownames(mat7)), ]
# Order species
sp_order <- c("Homo sapiens","Mus musculus","Drosophila melanogaster",
              "Caenorhabditis elegans","Saccharomyces cerevisiae","Arabidopsis thaliana")
mat7 <- mat7[, intersect(sp_order, colnames(mat7))]

col_fun7 <- colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B"))

svg(file.path(FIG_DIR, "fig7_cofactor_summary.svg"), width=7, height=4.5, family="Liberation Sans")
h7 <- Heatmap(mat7, name="Mean log2FC\n(cofactor-using genes)", col=col_fun7,
              cluster_rows=FALSE, cluster_columns=FALSE,
              row_names_side="left", row_names_gp=gpar(fontsize=8, fontfamily="Liberation Sans"),
              column_names_gp=gpar(fontsize=8, fontfamily="Liberation Sans", fontface="italic"),
              column_title="Cofactor-dependent enzyme expression across species",
              column_title_gp=gpar(fontsize=10, fontfamily="Liberation Sans", fontface="bold"),
              cell_fun=function(j,i,x,y,width,height,fill) {
                grid.text(sprintf("%.2f", mat7[i,j]), x, y, gp=gpar(fontsize=6, fontfamily="Liberation Sans",
                          col=ifelse(abs(mat7[i,j])>0.25, "white", "black")))
              })
draw(h7)
dev.off()

png(file.path(FIG_DIR, "fig7_cofactor_summary.png"), width=7, height=4.5, units="in", res=300, family="Liberation Sans")
draw(h7)
dev.off()
cat("  Saved fig7_cofactor_summary\n")

cat("\nAll figures 1-4, 7 done.\n")
