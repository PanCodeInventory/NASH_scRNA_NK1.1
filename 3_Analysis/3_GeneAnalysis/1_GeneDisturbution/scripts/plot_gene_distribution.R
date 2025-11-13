#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(optparse)
  library(patchwork)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("--rds"), type = "character", help = "Input Seurat RDS path", metavar = "file"),
  make_option(c("--geneset_dir"), type = "character", default = "3_Analysis/4_GeneDisturbution",
              help = "Directory containing Geneset_*.txt [default %default]"),
  make_option(c("--outdir"), type = "character", default = "3_Analysis/4_GeneDisturbution/results",
              help = "Output directory for results [default %default]"),
  make_option(c("--group-var"), type = "character", default = "seurat_clusters",
              help = "Metadata column to group cells by (e.g., seurat_clusters) [default %default]"),
  make_option(c("--split-var"), type = "character", default = "auto",
              help = "Metadata column to split by (e.g., timepoint). 'auto' to autodetect; '' to disable [default %default]"),
  make_option(c("--max-featureplots"), type = "integer", default = 12,
              help = "Max number of features for FeaturePlot per geneset [default %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "DPI for saved figures [default %default]"),
  make_option(c("--mode"), type = "character", default = "umap_only",
              help = "Plot mode: 'umap_only' for one-gene UMAP only; 'all' to include DotPlot/VlnPlot [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$rds), file.exists(opt$rds))

# ---------- Helpers ----------

msg <- function(...) cat(sprintf("[INFO] %s\n", sprintf(...)))
warn <- function(...) cat(sprintf("[WARN] %s\n", sprintf(...)))

ensure_dir <- function(d) { dir.create(d, showWarnings = FALSE, recursive = TRUE); d }

read_geneset <- function(path) {
  if (!file.exists(path)) return(character(0))
  txt <- read_file(path)
  if (is.na(txt) || txt == "") return(character(0))
  # split by comma, newline, semicolon, tab, space
  genes <- unlist(strsplit(txt, "[\n,;\t ]+", perl = TRUE))
  genes <- unique(genes[genes != "" & !is.na(genes)])
  genes
}

present_cols <- function(md, candidates) candidates[candidates %in% colnames(md)]

autodetect_split_var <- function(md) {
  cands <- c("timepoint", "Timepoint", "time", "group", "condition", "orig.ident", "dataset", "sample")
  cands <- present_cols(md, cands)
  if (length(cands) == 0) return(NA_character_)
  cands[[1]]
}

ensure_umap <- function(obj, dims = 1:10) {
  has_umap <- tryCatch({Embeddings(obj, reduction = "umap"); TRUE}, error = function(e) FALSE)
  if (has_umap) return(obj)
  msg("UMAP not found; attempting to compute using PCA dims %s", paste(dims, collapse=","))
  has_pca <- tryCatch({Embeddings(obj, reduction = "pca"); TRUE}, error = function(e) FALSE)
  if (!has_pca) {
    DefaultAssay(obj) <- DefaultAssay(obj)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(dims), verbose = FALSE)
  }
  obj <- RunUMAP(obj, dims = dims, reduction = "pca", verbose = FALSE)
  obj
}

save_plot <- function(p, path, width = 10, height = 6, dpi = 300) {
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi, limitsize = FALSE)
  msg("Saved: %s", path)
}

subset_features <- function(obj, features) {
  feats <- intersect(features, rownames(obj))
  missing <- setdiff(features, feats)
  if (length(missing) > 0) warn("Missing features (%d): %s", length(missing), paste(missing, collapse=", "))
  feats
}

# ---------- Main ----------

msg("Loading RDS: %s", opt$rds)
obj <- readRDS(opt$rds)
md <- obj@meta.data

plots_dir <- file.path(opt$outdir, "plots")
ensure_dir(plots_dir)

# Choose expression assay for plotting to maximize gene availability
assays_avail <- Assays(obj)
expr_assay <- if ("SCT" %in% assays_avail) "SCT" else if ("RNA" %in% assays_avail) "RNA" else DefaultAssay(obj)
DefaultAssay(obj) <- expr_assay
msg("Using expression assay for plotting: %s", expr_assay)

# Resolve grouping/splitting variables
if (!(opt$`group-var` %in% colnames(md))) {
  stop(sprintf("group-var '%s' not found in meta.data. Available: %s", opt$`group-var`, paste(colnames(md), collapse=", ")))
}

split_var <- opt$`split-var`
if (identical(split_var, "auto")) {
  sv <- autodetect_split_var(md)
  if (is.na(sv)) {
    msg("No split-var detected; proceeding without split.")
    split_var <- ""
  } else {
    split_var <- sv
    msg("Auto-detected split-var: %s", split_var)
  }
}

# Read gene sets
alias_map <- c("Cd49a"="Itga1","Cd49b"="Itga2","Cd127"="Il7r","Cd200r"="Cd200r1","Hobit"="Zfp683")
map_alias <- function(v) {
  out <- ifelse(v %in% names(alias_map), alias_map[v], v)
  as.character(out)
}
chemokines <- read_geneset(file.path(opt$geneset_dir, "Geneset_Chemokines.txt"))
markers_raw <- read_geneset(file.path(opt$geneset_dir, "Geneset_Marker.txt"))
markers <- unique(map_alias(markers_raw))
if (length(markers_raw) > 0 && !identical(sort(markers_raw), sort(markers))) {
  msg("Marker aliases applied: %s -> %s",
      paste(markers_raw, collapse=", "),
      paste(markers, collapse=", "))
}
specific   <- read_geneset(file.path(opt$geneset_dir, "Geneset_Specific.txt"))

msg("Geneset sizes: Chemokines=%d, Marker=%d, Specific=%d", length(chemokines), length(markers), length(specific))

# Ensure UMAP for FeaturePlot
obj <- ensure_umap(obj)

# Set identities to group-var for grouping
Idents(obj) <- factor(md[[opt$`group-var`]])

# Optionally create combined grouping for split visualization via DotPlot
if (!identical(split_var, "") && split_var %in% colnames(md)) {
  md$tp_cluster <- interaction(md[[split_var]], md[[opt$`group-var`]], drop = TRUE, sep = "|")
  obj@meta.data <- md
}

make_dotplot <- function(obj, feats, group_var, title_suffix) {
  if (length(feats) == 0) return(NULL)
  feats <- subset_features(obj, feats)
  if (length(feats) == 0) return(NULL)
  DotPlot(obj, features = feats, group.by = group_var) + 
    RotatedAxis() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = sprintf("DotPlot - %s", title_suffix), x = "Genes", y = group_var)
}

make_vlnplot <- function(obj, feats, group_var, title_suffix) {
  if (length(feats) == 0) return(NULL)
  feats <- subset_features(obj, feats)
  if (length(feats) == 0) return(NULL)
  p <- VlnPlot(obj, features = feats, group.by = group_var, pt.size = 0.1, combine = TRUE)
  p + theme_bw() + labs(title = sprintf("VlnPlot - %s", title_suffix))
}

make_featureplot_grid <- function(obj, feats, split_var, nmax = 12, title_prefix = "FeaturePlot") {
  if (length(feats) == 0) return(NULL)
  feats <- subset_features(obj, feats)
  if (length(feats) == 0) return(NULL)
  feats <- head(feats, nmax)
  plots <- lapply(feats, function(g) {
    if (!identical(split_var, "") && split_var %in% colnames(obj@meta.data)) {
      FeaturePlot(obj, features = g, split.by = split_var) + ggtitle(g)
    } else {
      FeaturePlot(obj, features = g) + ggtitle(g)
    }
  })
  wrap_plots(plots, ncol = 3) + plot_annotation(title = sprintf("%s (%s)", title_prefix, paste(feats, collapse=", ")))
}

# Helper: sanitize gene name for filenames
sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# Helper: single-gene FeaturePlot (no splitting)
make_featureplot_single <- function(obj, gene) {
  FeaturePlot(obj, features = gene) + ggtitle(gene)
}

# ----- Generate plots for each geneset -----

genesets <- list(
  Chemokines = chemokines,
  Marker = markers,
  Specific = specific
)

for (nm in names(genesets)) {
  feats <- genesets[[nm]]
  if (length(feats) == 0) { warn("Geneset '%s' is empty; skipping", nm); next }

  # UMAP: one gene per image, no splitting by samples/timepoint
  feats_plot <- subset_features(obj, feats)
  if (length(feats_plot) > 0) {
    for (g in feats_plot) {
      p1 <- make_featureplot_single(obj, g)
      out <- file.path(plots_dir, sprintf("UMAP_%s_%s.png", nm, sanitize_filename(g)))
      save_plot(p1, out, dpi = opt$`dpi`, width = 7, height = 6)
    }
  }

  if (tolower(opt$mode) != "umap_only") {
    # Keep DotPlot/VlnPlot outputs
    p_dot <- make_dotplot(obj, feats, opt$`group-var`, sprintf("%s by %s", nm, opt$`group-var`))
    if (!is.null(p_dot)) save_plot(p_dot, file.path(plots_dir, sprintf("DotPlot_%s_by_%s.png", nm, opt$`group-var`)), dpi = opt$`dpi`, width = 12, height = 6)

    if (!identical(split_var, "") && split_var %in% colnames(md)) {
      p_dot2 <- make_dotplot(obj, feats, "tp_cluster", sprintf("%s by %s|%s", nm, split_var, opt$`group-var`))
      if (!is.null(p_dot2)) save_plot(p_dot2, file.path(plots_dir, sprintf("DotPlot_%s_by_%s_%s.png", nm, split_var, opt$`group-var`)), dpi = opt$`dpi`, width = 14, height = 7)
    }

    p_vln <- make_vlnplot(obj, feats, opt$`group-var`, sprintf("%s by %s", nm, opt$`group-var`))
    if (!is.null(p_vln)) save_plot(p_vln, file.path(plots_dir, sprintf("VlnPlot_%s_by_%s.png", nm, opt$`group-var`)), dpi = opt$`dpi`, width = 14, height = 8)
  }
}

msg("All done. Outputs in: %s", plots_dir)
