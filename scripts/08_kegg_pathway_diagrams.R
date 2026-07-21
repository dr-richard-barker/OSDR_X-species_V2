#!/usr/bin/env Rscript
# ggkegg pathway diagrams with per-species log2FC overlay + cofactor annotations
suppressMessages({
  library(ggkegg); library(ggraph); library(igraph); library(tidygraph)
  library(ggplot2); library(dplyr); library(data.table); library(viridis)
  library(ggnewscale); library(patchwork); library(stringr)
})

options(datatable.optimize=1)  # avoid GForce issues

# ---- Config ----
FIG_DIR <- "/mnt/shared-workspace/shared/results/figures"
KGML_DIR <- "/workspace/kegg_cache/kgml"
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(KGML_DIR, showWarnings=FALSE, recursive=TRUE)

FC_SCALE <- scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
  midpoint=0, limits=c(-2,2), oob=scales::squish, name="mean log2FC\n(spaceflight vs ground)")

# ---- Load data ----
fc_long <- fread("/mnt/shared-workspace/shared/results/meta/pathway_gene_log2fc.csv")
fc_long[, entrez := as.character(entrez)]
fc_long[, n_deg := as.integer(n_deg)]
fc_long[, min_padj := as.numeric(min_padj)]
cof_gene <- fread("/mnt/shared-workspace/shared/results/meta/cofactor_gene_mapping.csv")
cof_gene[, entrez := as.character(entrez)]

# ---- Helper: parse hsa node name into entrez list ----
parse_hsa_names <- function(name_str) {
  if (is.na(name_str)) return(character(0))
  ids <- str_extract_all(name_str, "hsa:\\d+")[[1]]
  gsub("hsa:", "", ids)
}

# ---- Helper: safe process_line ----
safe_process_line <- function(g) {
  nd <- g %>% activate(nodes) %>% as_tibble()
  has_line <- any(nd$type == "line", na.rm=TRUE)
  if (has_line) {
    tryCatch(process_line(g), error=function(e) g)
  } else {
    g
  }
}

# ---- Helper: build per-pathway multi-species plot ----
plot_pathway <- function(pid, pname, fc_data, species_order) {
  cat(sprintf("Building %s (%s)...\n", pid, pname))
  g <- pathway(pid, directory=KGML_DIR, use_cache=TRUE)
  g <- safe_process_line(g)

  nd <- g %>% activate(nodes) %>% as_tibble()
  cat(sprintf("  Nodes: %d (types: %s)\n", nrow(nd),
              paste(names(table(nd$type)), collapse=",")))

  # For each gene node, parse entrez IDs and compute per-species mean FC
  node_fc_list <- lapply(seq_len(nrow(nd)), function(i) {
    nm <- nd$name[i]
    if (is.na(nm) || nd$type[i] != "gene") return(NULL)
    ids <- parse_hsa_names(nm)
    if (length(ids)==0) return(NULL)
    sub <- fc_data[entrez %in% ids]
    if (nrow(sub)==0) return(NULL)
    spec_fc <- sub[, .(mean_fc=mean(mean_log2fc, na.rm=TRUE),
                       n_deg=sum(n_deg, na.rm=TRUE),
                       min_padj=min(min_padj, na.rm=TRUE)), by=species]
    spec_fc[, node_idx := i]
    return(spec_fc)
  })
  node_fc <- rbindlist(node_fc_list[!sapply(node_fc_list, is.null)])
  if (nrow(node_fc)==0) {
    cat(sprintf("  No FC data for any gene node in %s\n", pid))
    return(NULL)
  }
  cat(sprintf("  Genes with FC: %d nodes across %d species\n",
              node_fc[, uniqueN(node_idx)], node_fc[, uniqueN(species)]))

  plots <- lapply(species_order, function(sp) {
    sp_fc <- node_fc[species==sp]
    if (nrow(sp_fc)==0) return(NULL)
    g2 <- g
    # Build FC vector indexed by node position
    fc_vec <- rep(NA_real_, vcount(g2))
    deg_vec <- rep(0L, vcount(g2))
    for (j in seq_len(nrow(sp_fc))) {
      idx <- sp_fc$node_idx[j]
      fc_vec[idx] <- sp_fc$mean_fc[j]
      deg_vec[idx] <- as.integer(sp_fc$n_deg[j])
    }
    g2 <- g2 %>% activate(nodes) %>%
      mutate(mean_fc = fc_vec, n_deg_node = deg_vec)

    p <- tryCatch({
      ggraph(g2, x=x, y=y) +
        geom_node_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                           fill=mean_fc, filter=type=="gene" & !is.na(mean_fc)),
                       linewidth=0.1) +
        geom_node_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                           filter=type=="gene" & is.na(mean_fc)),
                       fill="grey90", linewidth=0.1) +
        geom_edge_link0(width=0.2, color="grey70") +
        geom_node_point(aes(x=x, y=y, filter=type!="gene" & type!="line"),
                        size=1, color="grey50") +
        FC_SCALE +
        ggtitle(sp) +
        theme_void(base_family="Liberation Sans") +
        theme(plot.title=element_text(size=7, hjust=0.5, face="italic"),
              legend.position="none",
              plot.margin=margin(2,2,2,2))
    }, error=function(e) {
      cat(sprintf("    plot error for %s: %s\n", sp, conditionMessage(e)))
      NULL
    })
    return(p)
  })
  plots <- plots[!sapply(plots, is.null)]
  if (length(plots)==0) return(NULL)
  cat(sprintf("  Panels: %d species\n", length(plots)))

  combined <- wrap_plots(plots, ncol=3) +
    plot_annotation(title=pname,
      subtitle=sprintf("KEGG %s | node fill = mean log2FC (spaceflight vs ground)", pid),
      theme=theme(plot.title=element_text(size=11, face="bold", family="Liberation Sans"),
                  plot.subtitle=element_text(size=7, family="Liberation Sans")))

  svg_path <- file.path(FIG_DIR, sprintf("fig6_pathway_%s.svg", pid))
  png_path <- file.path(FIG_DIR, sprintf("fig6_pathway_%s.png", pid))
  ggsave(svg_path, combined, width=14, height=10, units="in", limitsize=FALSE)
  ggsave(png_path, combined, width=14, height=10, units="in", dpi=200, limitsize=FALSE)
  cat(sprintf("  Saved: %s\n", basename(svg_path)))
  return(combined)
}

# ---- Run for all pathways ----
species_order <- c("Homo sapiens","Mus musculus","Drosophila melanogaster",
                   "Caenorhabditis elegans","Saccharomyces cerevisiae","Arabidopsis thaliana")

pathways <- list(
  "hsa00190" = "Oxidative phosphorylation",
  "hsa03050" = "Proteasome",
  "hsa03010" = "Ribosome",
  "hsa00020" = "TCA cycle"
)

results <- list()
for (pid in names(pathways)) {
  results[[pid]] <- plot_pathway(pid, pathways[[pid]], fc_long, species_order)
}
cat("\nDone.\n")
