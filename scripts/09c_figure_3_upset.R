#!/usr/bin/env Rscript
# ======================================================================
# Step 9c: Figure 3 (Orthology UpSet Plot)
#
# Description: Generate orthology UpSet plot showing shared ortholog groups
#   across species, using a custom ggplot2 implementation for clean SVG text
#
# Inputs: results/orthology/orthology_matrix_wide.csv
# Outputs: figures/fig3_orthology_upset.{svg,png}
#
# Language: R
# See METHODS.md for full parameter details.
# ======================================================================
suppressMessages({
  library(ggplot2); library(dplyr); library(data.table); library(viridis)
  library(patchwork)
})

FIG_DIR <- "/mnt/shared-workspace/shared/results/figures"

SPECIES_COLORS <- c(
  "Homo sapiens"            = "#440154FF",
  "Mus musculus"            = "#3B528BFF",
  "Drosophila melanogaster" = "#21918CFF",
  "Caenorhabditis elegans"  = "#5EC962FF",
  "Saccharomyces cerevisiae"= "#FDE725FF",
  "Arabidopsis thaliana"    = "#9B1946FF"
)

cat("Building Fig 3: Orthology UpSet (ggplot)...\n")
ortho_wide <- fread("/mnt/shared-workspace/shared/orthodb/orthology_matrix_wide.csv")
nonhuman_cols <- c("Mus musculus","Drosophila melanogaster","Caenorhabditis elegans",
                   "Saccharomyces cerevisiae","Arabidopsis thaliana")
all_cols <- c("Homo sapiens", nonhuman_cols)

# Build binary matrix
ortho_bin <- as.data.frame(ortho_wide[, ..nonhuman_cols])
ortho_bin[is.na(ortho_bin)] <- ""
ortho_bin[ortho_bin != ""] <- 1
ortho_bin[ortho_bin == ""] <- 0
for (col in nonhuman_cols) ortho_bin[[col]] <- as.integer(ortho_bin[[col]])
ortho_bin[["Homo sapiens"]] <- 1L

# Compute intersection sizes
n_total <- nrow(ortho_bin)
# All non-empty combinations
combo_counts <- ortho_bin[, all_cols] |>
  (\(x) as.data.table(x))() |>
  (\(dt) dt[, .N, by=all_cols])()
setorder(combo_counts, -N)
combo_counts <- head(combo_counts, 15)  # top 15 combos

# Build combination labels
combo_labels <- apply(combo_counts[, ..all_cols], 1, function(row) {
  present <- all_cols[row == 1]
  paste(present, collapse="\n")
})
combo_counts[, combo_label := combo_labels]
combo_counts[, combo_idx := seq_len(.N)]

# Melt for matrix display
combo_long <- melt(combo_counts[, ..all_cols], measure.vars=all_cols,
                   variable.name="species", value.name="present")
combo_long[, combo_idx := rep(seq_len(nrow(combo_counts)), each=length(all_cols))]
combo_long[, species := factor(species, levels=rev(all_cols))]

# Top panel: bar chart of combination sizes
p_top <- ggplot(combo_counts, aes(x=combo_idx, y=N)) +
  geom_col(fill="grey40", width=0.7) +
  geom_text(aes(label=N), vjust=-0.3, size=2.5, family="Liberation Sans") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
  labs(y="Ortholog group size") +
  theme_classic(base_family="Liberation Sans", base_size=9) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), plot.margin=margin(5,5,0,5))

# Matrix panel: dots for present species
p_matrix <- ggplot(combo_long, aes(x=combo_idx, y=species)) +
  geom_point(aes(color=ifelse(present==1, "yes", "no")), size=3) +
  scale_color_manual(values=c(yes="grey20", no="white"), guide="none") +
  # connect lines within each combo
  geom_line(data=combo_long[present==1], aes(group=combo_idx), color="grey20", linewidth=0.5) +
  theme_classic(base_family="Liberation Sans", base_size=9) +
  labs(x="Ortholog combination", y=NULL) +
  theme(axis.text.y=element_text(face="italic", size=8, color=unname(SPECIES_COLORS[rev(all_cols)])),
        axis.text.x=element_blank(), axis.ticks=element_blank(),
        panel.grid=element_blank(), plot.margin=margin(0,5,5,5))

# Side panel: set sizes (genes per species)
set_sizes <- sapply(all_cols, function(col) sum(ortho_bin[[col]]))
set_dt <- data.table(species=all_cols, n_genes=set_sizes)
set_dt[, species := factor(species, levels=rev(all_cols))]

p_side <- ggplot(set_dt, aes(x=n_genes, y=species, fill=species)) +
  geom_col(width=0.7) +
  geom_text(aes(label=n_genes), hjust=-0.1, size=2.5, family="Liberation Sans") +
  scale_fill_manual(values=SPECIES_COLORS, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0, 0.15))) +
  labs(x="Genes with ortholog") +
  theme_classic(base_family="Liberation Sans", base_size=9) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank(), plot.margin=margin(0,5,5,0))

# Combine: top + matrix, with side panel
p_combined <- (p_top + plot_spacer() + plot_layout(widths=c(4, 1))) /
  (p_matrix + p_side + plot_layout(widths=c(4, 1))) +
  plot_layout(heights=c(1, 2)) +
  plot_annotation(title="Ortholog group sharing across species (human-anchored, OrthoDB v12)",
                  theme=theme(plot.title=element_text(size=11, face="bold", family="Liberation Sans")))

ggsave(file.path(FIG_DIR, "fig3_orthology_upset.svg"), p_combined, width=8, height=5, units="in")
ggsave(file.path(FIG_DIR, "fig3_orthology_upset.png"), p_combined, width=8, height=5, units="in", dpi=300)
cat("  Saved fig3_orthology_upset\n")

# Print key stats
cat("\nKey orthology stats:\n")
cat(sprintf("  Total human genes with orthologs: %d\n", n_total))
cat(sprintf("  Genes with ortholog in all 6 species: %d\n", sum(rowSums(ortho_bin[, all_cols])==6)))
cat(sprintf("  Genes with ortholog in >=4 species: %d\n", sum(rowSums(ortho_bin[, all_cols])>=4)))
cat(sprintf("  Genes with ortholog in >=3 species: %d\n", sum(rowSums(ortho_bin[, all_cols])>=3)))
cat("\nSet sizes:\n"); print(set_dt)
cat("\nTop combinations:\n"); print(combo_counts[, .(combo_label, N)])
