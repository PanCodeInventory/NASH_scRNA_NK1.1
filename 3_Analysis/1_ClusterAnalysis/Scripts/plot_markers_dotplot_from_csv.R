#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--input"), type = "character", help = "Path to markers_all_clusters.csv"),
  make_option(c("--outdir_plots"), type = "character", default = "3_Analysis/1_ClusterAnalysis/plots/markers"),
  make_option(c("--outdir_tables"), type = "character", default = "3_Analysis/1_ClusterAnalysis/results/tables"),
  make_option(c("--topn"), type = "integer", default = 10),
  make_option(c("--pval"), type = "double", default = 0.05, help = "p_val_adj threshold"),
  make_option(c("--logfc"), type = "double", default = 0.25, help = "avg_log2FC threshold (>)"),
  make_option(c("--dpct"), type = "double", default = 0.10, help = "(pct.1 - pct.2) threshold (>)")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$input), file.exists(opt$input))

message(sprintf("[INFO] Reading: %s", opt$input))
df <- suppressMessages(readr::read_csv(opt$input, show_col_types = FALSE))

# Basic checks
req_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
missing <- setdiff(req_cols, colnames(df))
if (length(missing) > 0) {
  stop(sprintf("[ERROR] Missing required columns: %s", paste(missing, collapse=", ")))
}

# Ensure types
num_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
df <- df %>% mutate(across(all_of(num_cols), as.numeric))

# Cluster as character then factor (numeric-aware ordering if possible)
cl_try <- suppressWarnings(as.numeric(df$cluster))
if (all(!is.na(cl_try))) {
  cluster_levels_all <- sort(unique(cl_try))
  df <- df %>% mutate(cluster = factor(as.character(cluster), levels = as.character(cluster_levels_all)))
} else {
  cluster_levels_all <- sort(unique(as.character(df$cluster)))
  df <- df %>% mutate(cluster = factor(as.character(cluster), levels = cluster_levels_all))
}

# Compute delta percentage and filter by thresholds
sel <- df %>%
  mutate(dpct = pct.1 - pct.2) %>%
  filter(p_val_adj < opt$pval, avg_log2FC > opt$logfc, dpct > opt$dpct)

# Rank within cluster and take top N by avg_log2FC desc, then p_val_adj asc, then pct.1 desc
sel <- sel %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, desc(pct.1), .by_group = TRUE) %>%
  slice_head(n = opt$topn) %>%
  ungroup()

# If nothing selected, stop gracefully
if (nrow(sel) == 0) {
  stop("[ERROR] No markers passed the thresholds; try relaxing thresholds or decreasing dpct.")
}

# Create output dirs
if (!dir.exists(opt$outdir_plots)) dir.create(opt$outdir_plots, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(opt$outdir_tables)) dir.create(opt$outdir_tables, recursive = TRUE, showWarnings = FALSE)

# Export selected table
out_csv <- file.path(opt$outdir_tables, sprintf("dotplot_markers_top%d.csv", opt$topn))
readr::write_csv(sel, out_csv)
message(sprintf("[INFO] Selected markers table written: %s (%d rows)", out_csv, nrow(sel)))

# Prepare plotting data: unique gene label per (gene, cluster) for y-axis to avoid duplicates
plot_df <- sel %>%
  mutate(gene_label = sprintf("%s (C%s)", gene, as.character(cluster)))

# Recompute factor levels for cluster limited to selected clusters for a tidy x-axis
cluster_levels_sel <- levels(droplevels(plot_df$cluster))

# Order gene labels by cluster then avg_log2FC desc
plot_df <- plot_df %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(gene_label = forcats::fct_inorder(gene_label)) %>%
  ungroup()

y_count <- length(unique(plot_df$gene_label))
height <- max(6, min(30, 0.28 * y_count + 3))

p <- ggplot(plot_df, aes(x = factor(cluster, levels = cluster_levels_sel), y = gene_label)) +
  geom_point(aes(size = pct.1, color = avg_log2FC), alpha = 0.9) +
  scale_size_continuous(range = c(1.5, 8), breaks = scales::pretty_breaks(4)) +
  scale_color_gradient(low = "#6BAED6", high = "#DE2D26") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = sprintf("Top %d markers per cluster (p_adj<%.2g, log2FC>%.2f, dpct>%.2f)", opt$topn, opt$pval, opt$logfc, opt$dpct),
    x = "Cluster",
    y = "Gene (cluster)",
    size = "pct.1",
    color = "avg_log2FC"
  )

out_png <- file.path(opt$outdir_plots, sprintf("dotplot_markers_top%d.png", opt$topn))
message(sprintf("[INFO] Saving plot: %s (height=%.1f)", out_png, height))

ggsave(out_png, plot = p, width = 12, height = height, dpi = 300, limitsize = FALSE)
message("[INFO] All done.")
