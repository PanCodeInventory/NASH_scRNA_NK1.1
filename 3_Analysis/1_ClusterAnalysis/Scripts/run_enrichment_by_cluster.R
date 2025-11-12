#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  # nothing yet
})

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
  make_option(c("--input"), type = "character", default = "3_Analysis/1_ClusterAnalysis/data/markers_all_clusters.csv", help = "Input markers CSV path"),
  make_option(c("--outdir"), type = "character", default = "3_Analysis/1_ClusterAnalysis/results/enrichment", help = "Output directory"),
  make_option(c("--exclude"), type = "character", default = "6", help = "Comma-separated cluster ids to exclude (e.g., '6,9')"),
  make_option(c("--top"), type = "integer", default = 20L, help = "Top N terms after filtering"),
  make_option(c("--go_ont"), type = "character", default = "BP", help = "GO ontology: BP/CC/MF"),
  make_option(c("--filter_terms"), type = "character", default = "ribosom,mitochond", help = "Comma-separated keywords to remove from display e.g. 'ribosom,mitochond'"),
  make_option(c("--install_missing_pkgs"), type = "character", default = "true", help = "true/false to auto-install dependencies"),
  make_option(c("--logs"), type = "character", default = "logs", help = "Logs subdir name under outdir")
)

opt <- parse_args(OptionParser(option_list = option_list))

# normalize boolean
to_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("1", "true", "t", "yes", "y")
}

auto_install <- to_bool(opt$install_missing_pkgs)

# -----------------------
# Load required packages
# -----------------------
cran_pkgs <- c("data.table", "dplyr", "stringr", "readr", "ggplot2", "tibble", "tidyr", "glue", "forcats")
bioc_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE", "KEGGREST")
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
  library(DOSE)
  library(KEGGREST)
})

# -----------------------
# IO and dirs
# -----------------------
input <- opt$input
outdir <- opt$outdir
logs_dir <- file.path(outdir, opt$logs)
summary_dir <- file.path(outdir, "summary")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d-%H%M")
log_file <- file.path(logs_dir, paste0("enrichment_run_", stamp, ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("[INFO] Starting enrichment pipeline\n")
cat(glue("[INFO] Input: {input}\n"))
cat(glue("[INFO] Outdir: {outdir}\n"))
cat(glue("[INFO] Exclude clusters: {opt$exclude}\n"))
cat(glue("[INFO] GO ont: {opt$go_ont}, Top: {opt$top}\n"))
cat(glue("[INFO] Filter terms: {opt$filter_terms}\n"))

# -----------------------
# Read and preprocess markers
# -----------------------
stopifnot(file.exists(input))

# fread auto-detects; ensure needed columns exist after reading
markers <- suppressWarnings(data.table::fread(input)) %>% as_tibble()

required_cols <- c("gene", "cluster", "p_val_adj", "avg_log2FC")
missing_cols <- setdiff(required_cols, colnames(markers))
if (length(missing_cols) > 0) {
  stop(glue("Missing required columns in input: {paste(missing_cols, collapse=\", \")}"))
}

markers <- markers %>%
  transmute(
    gene = as.character(.data$gene),
    cluster = as.character(.data$cluster),
    p_val_adj = as.numeric(.data$p_val_adj),
    avg_log2FC = as.numeric(.data$avg_log2FC)
  )

# filter by thresholds
markers_flt <- markers %>%
  filter(!is.na(gene), gene != "", !is.na(p_val_adj), !is.na(avg_log2FC)) %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0)

# exclude clusters
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

# split gene lists per cluster
cluster_gene_sets <- markers_flt %>%
  distinct(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

# -----------------------
# ID mapping SYMBOL -> ENTREZ
# -----------------------
cat("[INFO] Mapping SYMBOL to ENTREZID...\n")
bitr_safely <- function(symbols) {
  if (length(symbols) == 0) return(tibble(SYMBOL = character(), ENTREZID = character()))
  res <- tryCatch({
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>% as_tibble()
  }, error = function(e) {
    message("bitr error: ", conditionMessage(e))
    tibble(SYMBOL = character(), ENTREZID = character())
  })
  distinct(res, SYMBOL, ENTREZID) %>% filter(!is.na(ENTREZID), nzchar(ENTREZID))
}

mapped_list <- lapply(seq_len(nrow(cluster_gene_sets)), function(i) {
  cl <- cluster_gene_sets$cluster[i]
  syms <- unlist(cluster_gene_sets$genes[i])
  mapdf <- bitr_safely(syms)
  n_in <- length(unique(syms)); n_map <- nrow(mapdf)
  cat(glue("[MAP] cluster={cl} input_symbols={n_in} mapped_ENTREZ={n_map} rate={sprintf('%.1f%%', 100*n_map/max(1,n_in))}\n"))
  list(cluster = cl, map = mapdf)
})

# background universe: union of all mapped ENTREZ across clusters
bg_entrez <- mapped_list %>% lapply(function(x) x$map$ENTREZID) %>% unlist(use.names = FALSE) %>% unique()
cat(glue("[INFO] Background universe size (ENTREZ): {length(bg_entrez)}\n"))

# -----------------------
# Enrichment wrappers
# -----------------------
run_enrichGO_bp <- function(entrez) {
  if (length(entrez) < 1) return(NULL)
  tryCatch({
    clusterProfiler::enrichGO(
      gene = entrez,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = opt$go_ont,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      universe = bg_entrez,
      readable = FALSE
    )
  }, error = function(e) {
    message("enrichGO error: ", conditionMessage(e))
    NULL
  })
}

# KEGG: prefer enrichKEGG with ENTREZ/ncbi-geneid; fallback to enricher with KEGGREST TERM2GENE
prepare_kegg_term2gene <- function() {
  cat("[INFO] Preparing KEGG TERM2GENE via KEGGREST...\n")
  # Map gene -> pathway
  link <- KEGGREST::keggLink("pathway", "mmu") # returns named character: names are like "mmu:1234", values like "path:mmuXXXXX"
  gene_ids <- sub("^mmu:", "", names(link))
  pathway_ids <- sub("^path:", "", as.character(link))
  term2gene <- tibble(term = pathway_ids, gene = gene_ids)
  term2name <- KEGGREST::keggList("pathway", "mmu")
  term2name <- tibble(term = sub("^path:", "", names(term2name)), name = as.character(term2name))
  list(T2G = term2gene, T2N = term2name)
}

run_enrichKEGG_mmu <- function(entrez, t2g_bundle = NULL) {
  if (length(entrez) < 1) return(NULL)
  # First try enrichKEGG
  ek <- tryCatch({
    clusterProfiler::enrichKEGG(
      gene = entrez,
      organism = "mmu",
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = bg_entrez
    )
  }, error = function(e) {
    message("enrichKEGG error: ", conditionMessage(e))
    NULL
  })
  if (!is.null(ek)) return(ek)
  # Fallback: enricher with TERM2GENE
  if (is.null(t2g_bundle)) t2g_bundle <- prepare_kegg_term2gene()
  t2g <- t2g_bundle$T2G; t2n <- t2g_bundle$T2N
  # filter to background space to keep universe consistent
  t2g_use <- t2g %>% filter(gene %in% unique(c(bg_entrez, entrez)))
  ek2 <- tryCatch({
    clusterProfiler::enricher(
      gene = entrez,
      TERM2GENE = t2g_use,
      TERM2NAME = t2n,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      universe = bg_entrez
    )
  }, error = function(e) {
    message("enricher(KEGG) error: ", conditionMessage(e))
    NULL
  })
  ek2
}

# -----------------------
# Helpers for output
# -----------------------
# Convert GeneRatio like "3/45" to numeric
ratio_to_numeric <- function(x) {
  ifelse(grepl("/", x), {
    sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
  }, as.numeric(x))
}

save_df_csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) {
    readr::write_csv(tibble(), path)
  } else {
    readr::write_csv(df, path)
  }
}

save_plot_dual <- function(p, basepath, width = 9, height = 7) {
  ggplot2::ggsave(filename = paste0(basepath, ".pdf"), plot = p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(filename = paste0(basepath, ".png"), plot = p, width = width, height = height, dpi = 300)
}

wrap_terms <- function(x, width = 50) stringr::str_wrap(x, width = width)

# filtering pattern
filter_keywords <- unique(trimws(unlist(strsplit(opt$filter_terms, ","))))
filter_keywords <- filter_keywords[nzchar(filter_keywords)]
filter_regex <- if (length(filter_keywords) > 0) paste(filter_keywords, collapse = "|") else "^$NEVERMATCH$"

# containers for summaries
summary_go <- list()
summary_kegg <- list()

# prepare KEGG bundle once for fallback path
kegg_bundle <- NULL

# -----------------------
# Per-cluster enrichment and outputs
# -----------------------
for (cl in clusters) {
  cl_sanit <- gsub("[^A-Za-z0-9_-]+", "_", cl)
  cat(glue("[INFO] Processing cluster {cl} ...\n"))
  out_cluster_dir <- file.path(outdir, paste0("cluster_", cl_sanit))
  if (!dir.exists(out_cluster_dir)) dir.create(out_cluster_dir, recursive = TRUE, showWarnings = FALSE)

  entrez <- mapped_list[[which(clusters == cl)]]$map$ENTREZID %>% unique()

  # GO BP
  ego <- run_enrichGO_bp(entrez)
  ego_df <- if (!is.null(ego)) as.data.frame(ego) else tibble()
  go_full_path <- file.path(out_cluster_dir, "GO_BP_full.csv")
  save_df_csv(ego_df, go_full_path)

  ego_top <- tibble()
  if (nrow(ego_df) > 0) {
    ego_top <- ego_df %>%
      filter(!str_detect(Description, regex(filter_regex, ignore_case = TRUE))) %>%
      arrange(p.adjust, pvalue) %>%
      slice_head(n = opt$top)
    go_top_path <- file.path(out_cluster_dir, "GO_BP_top20_filtered.csv")
    save_df_csv(ego_top, go_top_path)

    # dotplot for filtered set (build from data.frame to avoid S4 method issues)
    if (nrow(ego_top) > 0) {
      ego_top2 <- ego_top %>%
        mutate(GeneRatioNum = ratio_to_numeric(GeneRatio)) %>%
        mutate(Description_wrapped = wrap_terms(Description, width = 50)) %>%
        arrange(p.adjust, pvalue)
      p_go <- ggplot(ego_top2, aes(x = GeneRatioNum, y = forcats::fct_reorder(Description_wrapped, GeneRatioNum))) +
        geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.8) +
        scale_color_viridis_c(name = "-log10(p.adjust)") +
        scale_size_continuous(name = "GeneRatio", range = c(1.5, 8)) +
        labs(x = "GeneRatio", y = "Term", title = glue("Cluster {cl} GO {opt$go_ont} (Top{opt$top} filtered)")) +
        theme_bw(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
      save_plot_dual(p_go, file.path(out_cluster_dir, "GO_BP_dotplot_top20_filtered"), width = 9, height = 7)
    }

    # add to summary container
    if (nrow(ego_top) > 0) {
      summary_go[[cl]] <- ego_top %>%
        transmute(cluster = cl, ID, Description, GeneRatio = ratio_to_numeric(GeneRatio), p.adjust)
    }
  } else {
    cat(glue("[WARN] Cluster {cl} has no significant GO terms.\n"))
    # still emit empty top file
    save_df_csv(ego_top, file.path(out_cluster_dir, "GO_BP_top20_filtered.csv"))
  }

  # KEGG
  eke <- run_enrichKEGG_mmu(entrez, t2g_bundle = kegg_bundle)
  if (is.null(eke) && is.null(kegg_bundle)) {
    # try prepare once more if not already
    kegg_bundle <- prepare_kegg_term2gene()
    eke <- run_enrichKEGG_mmu(entrez, t2g_bundle = kegg_bundle)
  }
  eke_df <- if (!is.null(eke)) as.data.frame(eke) else tibble()
  kegg_full_path <- file.path(out_cluster_dir, "KEGG_full.csv")
  save_df_csv(eke_df, kegg_full_path)

  eke_top <- tibble()
  if (nrow(eke_df) > 0) {
    eke_top <- eke_df %>%
      filter(!str_detect(Description, regex(filter_regex, ignore_case = TRUE))) %>%
      arrange(p.adjust, pvalue) %>%
      slice_head(n = opt$top)
    kegg_top_path <- file.path(out_cluster_dir, "KEGG_top20_filtered.csv")
    save_df_csv(eke_top, kegg_top_path)

    # dotplot for filtered set (build from data.frame to avoid S4 method issues)
    if (nrow(eke_top) > 0) {
      eke_top2 <- eke_top %>%
        mutate(GeneRatioNum = ratio_to_numeric(GeneRatio)) %>%
        mutate(Description_wrapped = wrap_terms(Description, width = 50)) %>%
        arrange(p.adjust, pvalue)
      p_kegg <- ggplot(eke_top2, aes(x = GeneRatioNum, y = forcats::fct_reorder(Description_wrapped, GeneRatioNum))) +
        geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.8) +
        scale_color_viridis_c(name = "-log10(p.adjust)") +
        scale_size_continuous(name = "GeneRatio", range = c(1.5, 8)) +
        labs(x = "GeneRatio", y = "Pathway", title = glue("Cluster {cl} KEGG (Top{opt$top} filtered)")) +
        theme_bw(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
      save_plot_dual(p_kegg, file.path(out_cluster_dir, "KEGG_dotplot_top20_filtered"), width = 9, height = 7)
    }

    if (nrow(eke_top) > 0) {
      summary_kegg[[cl]] <- eke_top %>%
        transmute(cluster = cl, ID, Description, GeneRatio = ratio_to_numeric(GeneRatio), p.adjust)
    }
  } else {
    cat(glue("[WARN] Cluster {cl} has no significant KEGG terms.\n"))
    save_df_csv(eke_top, file.path(out_cluster_dir, "KEGG_top20_filtered.csv"))
  }
}

# -----------------------
# Summary bubble plots
# -----------------------
make_bubble <- function(df, title) {
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(term = Description)
  # order terms by mean -log10(p.adjust)
  df <- df %>% group_by(term) %>% mutate(score = -log10(p.adjust)) %>%
    mutate(avg_score = mean(score, na.rm = TRUE)) %>% ungroup()
  term_levels <- df %>% distinct(term, avg_score) %>% arrange(desc(avg_score)) %>% pull(term)
  df <- df %>% mutate(term = factor(wrap_terms(term, width = 50), levels = wrap_terms(term_levels, width = 50)))

  ggplot(df, aes(x = cluster, y = term)) +
    geom_point(aes(size = GeneRatio, color = -log10(p.adjust)), alpha = 0.8) +
    scale_color_viridis_c(name = "-log10(p.adjust)") +
    scale_size_continuous(name = "GeneRatio", range = c(1.5, 8)) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    ) +
    labs(x = "Cluster", y = "Term", title = title)
}

# GO summary
sum_go_df <- if (length(summary_go) > 0) bind_rows(summary_go) else tibble()
if (nrow(sum_go_df) > 0) {
  readr::write_csv(sum_go_df, file.path(summary_dir, "GO_BP_summary_top_terms.csv"))
  p_sum_go <- make_bubble(sum_go_df, title = glue("GO {opt$go_ont} summary (Top{opt$top} filtered across clusters)"))
  if (!is.null(p_sum_go)) {
    save_plot_dual(p_sum_go, file.path(summary_dir, "GO_BP_bubbleplot"), width = 10, height = 8)
  }
} else {
  cat("[WARN] No GO summary terms available.\n")
  readr::write_csv(tibble(), file.path(summary_dir, "GO_BP_summary_top_terms.csv"))
}

# KEGG summary
sum_kegg_df <- if (length(summary_kegg) > 0) bind_rows(summary_kegg) else tibble()
if (nrow(sum_kegg_df) > 0) {
  readr::write_csv(sum_kegg_df, file.path(summary_dir, "KEGG_summary_top_terms.csv"))
  p_sum_kegg <- make_bubble(sum_kegg_df, title = glue("KEGG summary (Top{opt$top} filtered across clusters)"))
  if (!is.null(p_sum_kegg)) {
    save_plot_dual(p_sum_kegg, file.path(summary_dir, "KEGG_bubbleplot"), width = 10, height = 8)
  }
} else {
  cat("[WARN] No KEGG summary terms available.\n")
  readr::write_csv(tibble(), file.path(summary_dir, "KEGG_summary_top_terms.csv"))
}

# -----------------------
# README
# -----------------------
readme_path <- file.path(outdir, "README.md")
readme_txt <- glue::glue(
  "# Cluster-wise Enrichment Report\n\n",
  "- Data source: {input}\n",
  "- Species: Mus musculus (OrgDb=org.Mm.eg.db; KEGG=mmu)\n",
  "- Marker filter: p_val_adj < 0.05 and avg_log2FC > 0\n",
  "- Excluded clusters: {ifelse(length(exclude_vec)>0, paste(exclude_vec, collapse=\", \") , \"none\")}\n",
  "- GO ontology: {opt$go_ont} (clusterProfiler::enrichGO)\n",
  "- KEGG: clusterProfiler::enrichKEGG (fallback to enricher+KEGGREST if needed)\n",
  "- Background (universe): union of successfully mapped ENTREZ across all clusters\n",
  "- Display filtering keywords (excluded from TopN only): {opt$filter_terms}\n",
  "- Top N shown per cluster after filtering: {opt$top}\n",
  "- Outputs: per-cluster CSV and dotplots; summary bubble plots and tables.\n",
  "- Log file: {log_file}\n\n",
  "## Reproduce\n",
  "Rscript {normalizePath(sys.frame(1)$input, mustWork = FALSE)} \\\n    --input {input} \\\n    --outdir {outdir} \\\n    --exclude {opt$exclude} \\\n    --top {opt$top} \\\n    --go_ont {opt$go_ont} \\\n    --filter_terms \"{opt$filter_terms}\"\n\n",
  "## Notes\n",
  "- Terms containing '{opt$filter_terms}' are removed only from displayed TopN and exported '*_top20_filtered.csv'; full tables remain unmodified.\n",
  "- If a cluster yields no significant terms, empty files are still emitted.\n",
  "- KEGG may use a fallback implementation depending on local package versions.\n"
)
writeLines(readme_txt, con = readme_path)

cat("[INFO] Finished.\n")

# close log sinks
sink(type = "message")
sink()
close(log_con)
