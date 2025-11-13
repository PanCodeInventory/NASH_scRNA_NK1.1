#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- ifelse(length(args) >= 1, args[1], "3_Analysis/2_TimepointAnalysis/2_SpecificGeneVariation/scripts/config.yaml")

message("[INFO] Using config: ", config_path)
cfg <- yaml::read_yaml(config_path)

# Resolve outputs
outputs <- cfg$outputs
plots_dir <- outputs$plots_dir
tables_dir <- outputs$tables_dir
log_dir <- outputs$log_dir
base_dir <- outputs$base_dir

safe_dir_create <- function(p) {
  if (!dir.exists(p)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}
safe_dir_create(plots_dir)
safe_dir_create(tables_dir)
safe_dir_create(log_dir)
safe_dir_create(base_dir)

log_file <- file.path(log_dir, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
con <- file(log_file, open = "wt")
writeLog <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", file = con)
}
writeLog("Start Lgals1 timepoint variation analysis")
writeLog("RDS path:", cfg$rds_path)

# Load object
obj <- readRDS(cfg$rds_path)
writeLog("Seurat object loaded. Cells:", ncol(obj), " Genes:", nrow(obj))

# Validate gene and timepoint column
gene <- cfg$gene
tp_col <- cfg$timepoint_col
threshold <- as.numeric(cfg$expr_threshold)

if (!gene %in% rownames(obj)) {
  writeLog("ERROR: gene not found:", gene)
  close(con)
  stop(paste("Gene not found:", gene))
}

if (!tp_col %in% colnames(obj@meta.data)) {
  writeLog("ERROR: timepoint column not found:", tp_col)
  close(con)
  stop(paste("Timepoint column not found:", tp_col))
}

# Order timepoints if provided
if (!is.null(cfg$timepoint_order)) {
  tp_levels <- cfg$timepoint_order
  obj@meta.data[[tp_col]] <- factor(obj@meta.data[[tp_col]], levels = tp_levels)
}

# Prepare data
Idents(obj) <- obj@meta.data[[tp_col]]
fetch_df <- FetchData(obj, vars = c(gene, tp_col))
colnames(fetch_df) <- c("expr", "timepoint")

# Summary statistics by timepoint
summary_df <- fetch_df %>%
  group_by(timepoint) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    pct_expr_gt_threshold = mean(expr > threshold, na.rm = TRUE) * 100,
    n_cells = n(),
    .groups = "drop"
  )

summary_path <- file.path(tables_dir, paste0(gene, "_byTimepoint_summary.csv"))
writeLog("Writing summary table:", summary_path)
readr::write_csv(summary_df, summary_path)

# Plot functions
save_plot <- function(p, base_filename) {
  png_path <- file.path(plots_dir, paste0(base_filename, ".png"))
  svg_path <- file.path(plots_dir, paste0(base_filename, ".svg"))
  ggsave(png_path, plot = p, width = 8, height = 5, dpi = 300)
  ggsave(svg_path, plot = p, width = 8, height = 5)
  writeLog("Saved:", png_path)
  writeLog("Saved:", svg_path)
}

# Violin plot - all cells
writeLog("Generating violin plot (all cells)")
p_all <- VlnPlot(obj, features = gene, group.by = tp_col, pt.size = 0) +
  ggtitle(paste0(gene, " expression by ", tp_col, " (all cells)")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  ) + ylab(paste0(gene, " expression")) +
  geom_point(
    data = fetch_df,
    aes(x = timepoint, y = expr),
    color = "black", size = 0.2, alpha = 0.2,
    position = position_jitter(width = 0.2),
    inherit.aes = FALSE
  )

save_plot(p_all, paste0("Violin_", gene, "_byTimepoint_allcells"))

# Violin plot - expr > threshold
writeLog("Generating violin plot (expr > ", threshold, ")")
keep_cells <- rownames(fetch_df)[fetch_df$expr > threshold]
fetch_df_pos <- fetch_df[fetch_df$expr > threshold, ]
obj_pos <- obj[, keep_cells]
if (ncol(obj_pos) > 0) {
  Idents(obj_pos) <- obj_pos@meta.data[[tp_col]]
  p_pos <- VlnPlot(obj_pos, features = gene, group.by = tp_col, pt.size = 0) +
    ggtitle(paste0(gene, " expression by ", tp_col, " (expr > ", threshold, ")")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    ) + ylab(paste0(gene, " expression")) +
    geom_point(
      data = fetch_df_pos,
      aes(x = timepoint, y = expr),
      color = "black", size = 0.2, alpha = 0.2,
      position = position_jitter(width = 0.2),
      inherit.aes = FALSE
    )
  save_plot(p_pos, paste0("Violin_", gene, "_byTimepoint_expr_gt", threshold))
} else {
  writeLog("WARNING: No cells with expression > threshold. Skipping second plot.")
}

writeLog("Analysis completed successfully.")
close(con)
