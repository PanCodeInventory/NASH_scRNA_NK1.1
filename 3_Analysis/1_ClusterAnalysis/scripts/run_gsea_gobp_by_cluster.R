#!/usr/bin/env Rscript

# -----------------------
# Helper: package installer
# -----------------------
install_and_load <- function(pkgs_cran = c(), pkgs_bioc = c(), auto_install = TRUE) {
  load_one <- function(pkg) {
    suppressPackageStartupMessages(require(pkg, character.only = TRUE))
  }
  if (auto_install) {
    # CRAN
    for (p in pkgs_cran) {
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
      }
    }
    # Bioconductor
    if (length(pkgs_bioc) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
      }
      for (p in pkgs_bioc) {
        if (!requireNamespace(p, quietly = TRUE)) {
          BiocManager::install(p, ask = FALSE, quiet = TRUE, update = FALSE)
        }
      }
    }
  }
  # load
  ok <- vapply(c(pkgs_cran, pkgs_bioc), load_one, logical(1))
  invisible(all(ok))
}

# -----------------------
# CLI options
# -----------------------
optparse_pkgs <- c("optparse")
install_and_load(pkgs_cran = optparse_pkgs, pkgs_bioc = c(), auto_install = TRUE)
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--input"), type = "character", default = "3_Analysis/1_ClusterAnalysis/results/tables/markers_all_clusters.csv", help = "Input markers CSV path"),
  make_option(c("--outdir"), type = "character", default = "3_Analysis/1_ClusterAnalysis/results/enrichment_gsea_gobp", help = "Output directory for GSEA GO BP results"),
  make_option(c("--exclude"), type = "character", default = "6", help = "Comma-separated cluster ids to exclude (e.g., '6,9')"),
  make_option(c("--install_missing_pkgs"), type = "character", default = "true", help = "true/false to auto-install dependencies"),
  make_option(c("--logs"), type = "character", default = "logs", help = "Logs subdir name under outdir")
)

opt <- parse_args(OptionParser(option_list = option_list))

to_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("1", "true", "t", "yes", "y")
}
auto_install <- to_bool(opt$install_missing_pkgs)

# -----------------------
# Load required packages
# -----------------------
cran_pkgs <- c("data.table", "dplyr", "stringr", "readr", "ggplot2", "tibble", "tidyr", "glue", "forcats", "pheatmap")
bioc_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE")
install_and_load(cran_pkgs, bioc_pkgs, auto_install = auto_install)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(glue)
  library(forcats)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(pheatmap)
})

# -----------------------
# IO and dirs
# -----------------------
input <- opt$input
outdir <- opt$outdir
logs_dir <- file.path(outdir, opt$logs)
summary_dir <- file.path(outdir, "summary")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d-%H%M")
log_file <- file.path(logs_dir, paste0("gsea_gobp_run_", stamp, ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("[INFO] Starting GSEA GO BP pipeline\n")
cat(glue("[INFO] Input: {input}\n"))
cat(glue("[INFO] Outdir: {outdir}\n"))
cat(glue("[INFO] Exclude clusters: {opt$exclude}\n"))

# -----------------------
# Read and preprocess markers
# -----------------------
stopifnot(file.exists(input))
markers <- suppressWarnings(data.table::fread(input)) %>% as_tibble()

required_cols <- c("gene", "cluster", "avg_log2FC")
missing_cols <- setdiff(required_cols, colnames(markers))
if (length(missing_cols) > 0) {
  stop(glue("Missing required columns in input: {paste(missing_cols, collapse=\", \")}"))
}

markers_flt <- markers %>%
  transmute(
    gene = as.character(.data$gene),
    cluster = as.character(.data$cluster),
    avg_log2FC = as.numeric(.data$avg_log2FC)
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(avg_log2FC))

exclude_vec <- unique(strsplit(opt$exclude, ",")[[1]])
exclude_vec <- trimws(exclude_vec)
exclude_vec <- exclude_vec[nzchar(exclude_vec)]
if (length(exclude_vec) > 0) {
  markers_flt <- markers_flt %>% filter(!(.data$cluster %in% exclude_vec))
}

if (nrow(markers_flt) == 0) {
  stop("No markers remain after filtering and exclusion.")
}

clusters <- sort(unique(markers_flt$cluster))
cat(glue("[INFO] Clusters to analyze: {paste(clusters, collapse=\", \")}") %>% paste0("\n") )

# -----------------------
# GSEA wrapper
# -----------------------
run_gseGO <- function(gene_list) {
  if (length(gene_list) < 10) return(NULL)
  tryCatch({
    clusterProfiler::gseGO(
      geneList = gene_list,
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      nPerm = 10000,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
  }, error = function(e) {
    message("gseGO error: ", conditionMessage(e))
    NULL
  })
}

# -----------------------
# Helpers for output
# -----------------------
save_df_csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) {
    readr::write_csv(tibble(), path)
  } else {
    readr::write_csv(df, path)
  }
}

save_plot_dual <- function(p, basepath, width = 9, height = 7) {
  if(is.null(p)) return()
  ggplot2::ggsave(filename = paste0(basepath, ".pdf"), plot = p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(filename = paste0(basepath, ".png"), plot = p, width = width, height = height, dpi = 300)
}

# -----------------------
# Per-cluster GSEA and outputs
# -----------------------
summary_gsea <- list()

for (cl in clusters) {
  cl_sanit <- gsub("[^A-Za-z0-9_-]+", "_", cl)
  cat(glue("[INFO] Processing cluster {cl} ...\n"))
  out_cluster_dir <- file.path(outdir, paste0("cluster_", cl_sanit))
  dir.create(out_cluster_dir, recursive = TRUE, showWarnings = FALSE)

  # Prepare ranked gene list for the current cluster
  gene_list_df <- markers_flt %>%
    filter(cluster == cl) %>%
    distinct(gene, .keep_all = TRUE) %>% 
    arrange(desc(avg_log2FC))
  
  gene_list <- gene_list_df$avg_log2FC
  names(gene_list) <- gene_list_df$gene
  
  cat(glue("[INFO] Cluster {cl}: {length(gene_list)} genes in ranked list.\n"))

  # Run GSEA
  gsea_res <- run_gseGO(gene_list)
  gsea_df <- if (!is.null(gsea_res)) as.data.frame(gsea_res) else tibble()
  
  gsea_full_path <- file.path(out_cluster_dir, "GSEA_GO_BP_full.csv")
  save_df_csv(gsea_df, gsea_full_path)

  if (nrow(gsea_df) > 0) {
    gsea_top <- gsea_df %>% arrange(p.adjust) %>% slice_head(n = 20)
    gsea_top_path <- file.path(out_cluster_dir, "GSEA_GO_BP_top20.csv")
    save_df_csv(gsea_top, gsea_top_path)

    top_pathways <- head(gsea_res@result$ID, 5)
    if(length(top_pathways) > 0) {
        p_gsea <- enrichplot::gseaplot2(gsea_res, geneSetID = top_pathways, title = glue("GSEA GO BP for Cluster {cl}"))
        save_plot_dual(p_gsea, file.path(out_cluster_dir, "GSEA_plot_top5"), width = 10, height = 8)
    }

    p_dot <- enrichplot::dotplot(gsea_res, showCategory=15, split=".sign") + facet_grid(~.sign) +
             labs(title = glue("GSEA GO BP Dotplot for Cluster {cl}"))
    save_plot_dual(p_dot, file.path(out_cluster_dir, "GSEA_dotplot_top15"), width = 10, height = 8)

    summary_gsea[[cl]] <- gsea_df %>%
      transmute(cluster = cl, ID, Description, setSize, enrichmentScore, NES, p.adjust)
  } else {
    cat(glue("[WARN] Cluster {cl} has no significant GSEA GO BP terms.\n"))
  }
}

# -----------------------
# Summary bubble plots / heatmaps
# -----------------------
sum_gsea_df <- if (length(summary_gsea) > 0) bind_rows(summary_gsea) else tibble()
if (nrow(sum_gsea_df) > 0) {
  readr::write_csv(sum_gsea_df, file.path(summary_dir, "GSEA_GO_BP_summary_all_terms.csv"))

  top_terms_per_cluster <- sum_gsea_df %>%
    group_by(cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  wide_df <- top_terms_per_cluster %>%
    distinct(Description, cluster, .keep_all=TRUE) %>%
    select(Description, cluster, NES) %>%
    pivot_wider(names_from = cluster, values_from = NES, values_fill = 0) %>%
    column_to_rownames("Description")

  if(nrow(wide_df) > 1 && ncol(wide_df) > 1) {
    p_heatmap <- pheatmap(wide_df, 
                          cluster_rows = TRUE, 
                          cluster_cols = TRUE,
                          scale = "row",
                          main = "Heatmap of GSEA NES (Top 5 terms per cluster)")
    
    pdf(file.path(summary_dir, "GSEA_summary_heatmap.pdf"), width = 12, height = 10)
    print(p_heatmap)
    dev.off()
    png(file.path(summary_dir, "GSEA_summary_heatmap.png"), width = 1200, height = 1000, res=100)
    print(p_heatmap)
    dev.off()
  }
} else {
  cat("[WARN] No GSEA summary terms available.\n")
}

cat("[INFO] Finished.\n")
sink(type = "message")
sink()
close(log_con)
