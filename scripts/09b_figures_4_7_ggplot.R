#!/usr/bin/env Rscript
# ======================================================================
# Step 9b: Figures 4, 7 (Organelle Enrichment Heatmap, Cofactor Summary)
#
# Description: Generate Fig 4 (organelle x species enrichment heatmap) and
#   Fig 7 (cofactor x species fold-change heatmap) using ggplot2 geom_tile
#   for proper SVG text rendering
#
# Inputs: results/enrichment/organelle_enrichment_summary.csv,
#         results/meta/cofactor_gene_mapping.csv, results/meta/pathway_gene_log2fc.csv
# Outputs: figures/fig4_organelle_enrichment_heatmap.{svg,png},
#          figures/fig7_cofactor_summary.{svg,png}
#
# Language: R
# See METHODS.md for full parameter details.
# ======================================================================
suppressMessages({
  library(ggplot2); library(dplyr); library(data.table); library(viridis)
})

FIG_DIR <- "/mnt/shared-workspace/shared/results/figures"

theme_set(theme_classic(base_family="Liberation Sans", base_size=10) +
          theme(panel.grid=element_blank()))

# ============================================================
# Fig 4: Organelle enrichment heatmap (ggplot tile)
# ============================================================
cat("Building Fig 4 (ggplot)...\n")
enrich <- fread("/mnt/shared-workspace/shared/results/enrichment/organelle_enrichment_summary.csv")
# Sum signed scores per organelle x species
fig4_data <- enrich[, .(signed_score=sum(signed_score, na.rm=TRUE)), by=.(organelle, species)]
# Order organelles by total abs score
org_order <- fig4_data[, .(tot=sum(abs(signed_score))), by=organelle][order(-tot), organelle]
fig4_data[, organelle := factor(organelle, levels=rev(org_order))]
sp_order <- c("Homo sapiens","Mus musculus","Drosophila melanogaster",
              "Caenorhabditis elegans","Saccharomyces cerevisiae","Arabidopsis thaliana")
fig4_data[, species := factor(species, levels=sp_order)]

p4 <- ggplot(fig4_data, aes(x=species, y=organelle, fill=signed_score)) +
  geom_tile(color="white", linewidth=0.3) +
  geom_text(aes(label=ifelse(abs(signed_score)>5, sprintf("%.0f", signed_score), "")),
            size=2.5, family="Liberation Sans",
            color=ifelse(abs(fig4_data$signed_score)>30, "white", "black")) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                       name="Signed\n-log10(padj)\n(down | up)") +
  labs(title="Subcellular organelle enrichment across species",
       x=NULL, y=NULL,
       caption="Signed enrichment score: negative = downregulated enrichment, positive = upregulated enrichment") +
  theme(axis.text.x=element_text(angle=30, hjust=1, face="italic"),
        axis.text.y=element_text(size=8),
        plot.title=element_text(face="bold"),
        legend.position="right",
        plot.caption=element_text(size=7, color="grey40"))

ggsave(file.path(FIG_DIR, "fig4_organelle_enrichment_heatmap.svg"), p4, width=7, height=6, units="in")
ggsave(file.path(FIG_DIR, "fig4_organelle_enrichment_heatmap.png"), p4, width=7, height=6, units="in", dpi=300)
cat("  Saved fig4_organelle_enrichment_heatmap\n")

# ============================================================
# Fig 7: Cofactor summary heatmap (ggplot tile)
# ============================================================
cat("Building Fig 7 (ggplot)...\n")
cof_summary <- fread("/mnt/shared-workspace/shared/results/kegg/cofactor_pathway_summary.csv")
fig7_data <- cof_summary[, .(mean_fc=mean(mean_fc, na.rm=TRUE)), by=.(cofactor_name, species)]

cof_order <- c("NAD+","NADH","NADP+","NADPH","FAD","FADH2","FMN","CoA",
               "TPP (thiamine PP)","Lipoamide","Biotin","Heme")
fig7_data[, cofactor_name := factor(cofactor_name, levels=rev(cof_order))]
fig7_data[, species := factor(species, levels=sp_order)]

p7 <- ggplot(fig7_data, aes(x=species, y=cofactor_name, fill=mean_fc)) +
  geom_tile(color="white", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.2f", mean_fc)), size=2.5, family="Liberation Sans",
            color=ifelse(abs(fig7_data$mean_fc)>0.25, "white", "black")) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                       limits=c(-0.5,0.5), oob=scales::squish,
                       name="Mean log2FC\n(cofactor-using\ngenes)") +
  labs(title="Cofactor-dependent enzyme expression across species",
       x=NULL, y="Cofactor",
       caption="Mean log2FC of genes using each cofactor, averaged across pathways (OxPhos, Proteasome, Ribosome, TCA)") +
  theme(axis.text.x=element_text(angle=30, hjust=1, face="italic"),
        axis.text.y=element_text(size=8),
        plot.title=element_text(face="bold"),
        legend.position="right",
        plot.caption=element_text(size=7, color="grey40"))

ggsave(file.path(FIG_DIR, "fig7_cofactor_summary.svg"), p7, width=7, height=4.5, units="in")
ggsave(file.path(FIG_DIR, "fig7_cofactor_summary.png"), p7, width=7, height=4.5, units="in", dpi=300)
cat("  Saved fig7_cofactor_summary\n")

cat("\nDone.\n")
