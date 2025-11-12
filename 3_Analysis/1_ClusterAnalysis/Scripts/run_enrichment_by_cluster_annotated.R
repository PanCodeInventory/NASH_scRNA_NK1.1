#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  # nothing yet
})

# -----------------------
# Helper: package installer
# -----------------------
# 包安装和加载函数：自动检查、安装并加载所需的R包
install_and_load <- function(pkgs_cran = c(), pkgs_bioc = c(), auto_install = TRUE) {
  # 内部函数：加载单个包
  load_one <- function(pkg) {
    suppressPackageStartupMessages(require(pkg, character.only = TRUE))
  }
  # 如果启用自动安装
  if (auto_install) {
    # 安装CRAN包
    for (p in pkgs_cran) {
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
      }
    }
    # 安装Bioconductor包
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
  # 加载所有包
  ok <- vapply(c(pkgs_cran, pkgs_bioc), load_one, logical(1))
  invisible(all(ok))
}

# -----------------------
# CLI options
# -----------------------
# 命令行参数解析相关包
optparse_pkgs <- c("optparse")
install_and_load(pkgs_cran = optparse_pkgs, pkgs_bioc = c(), auto_install = TRUE)

suppressPackageStartupMessages(library(optparse))

# 定义命令行选项列表
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

# 解析命令行参数
opt <- parse_args(OptionParser(option_list = option_list))

# 布尔值标准化函数：将字符串转换为逻辑值
to_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("1", "true", "t", "yes", "y")
}

# 确定是否自动安装包
auto_install <- to_bool(opt$install_missing_pkgs)

# -----------------------
# Load required packages
# -----------------------
# 定义所需的CRAN和Bioconductor包
cran_pkgs <- c("data.table", "dplyr", "stringr", "readr", "ggplot2", "tibble", "tidyr", "glue", "forcats")
bioc_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE", "KEGGREST")
# 安装并加载所有必需包
install_and_load(cran_pkgs, bioc_pkgs, auto_install = auto_install)

# 加载包，抑制启动消息
suppressPackageStartupMessages({
  library(data.table)    # 高效数据读取和处理
  library(dplyr)          # 数据操作
  library(stringr)        # 字符串处理
  library(readr)          # 文件读写
  library(ggplot2)        # 数据可视化
  library(tidyr)          # 数据整理
  library(tibble)         # 现代数据框
  library(glue)           # 字符串插值
  library(forcats)        # 因子操作
  library(clusterProfiler) # 功能富集分析
  library(org.Mm.eg.db)   # 小鼠基因注释数据库
  library(enrichplot)     # 富集结果可视化
  library(DOSE)           # 疾病本体语义分析
  library(KEGGREST)       # KEGG数据库访问
})

# -----------------------
# IO and dirs
# -----------------------
# 设置输入输出路径
input <- opt$input
outdir <- opt$outdir
logs_dir <- file.path(outdir, opt$logs)
summary_dir <- file.path(outdir, "summary")

# 创建必要的目录
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# 设置日志文件
stamp <- format(Sys.time(), "%Y%m%d-%H%M")
log_file <- file.path(logs_dir, paste0("enrichment_run_", stamp, ".log"))
log_con <- file(log_file, open = "wt")
# 重定向输出到日志文件
sink(log_con, split = TRUE)
sink(log_con, type = "message")

# 打印运行信息
cat("[INFO] Starting enrichment pipeline\n")
cat(glue("[INFO] Input: {input}\n"))
cat(glue("[INFO] Outdir: {outdir}\n"))
cat(glue("[INFO] Exclude clusters: {opt$exclude}\n"))
cat(glue("[INFO] GO ont: {opt$go_ont}, Top: {opt$top}\n"))
cat(glue("[INFO] Filter terms: {opt$filter_terms}\n"))

# -----------------------
# Read and preprocess markers
# -----------------------
# 检查输入文件是否存在
stopifnot(file.exists(input))

# 使用data.table高效读取标记基因数据
markers <- suppressWarnings(data.table::fread(input)) %>% as_tibble()

# 检查必需的列是否存在
required_cols <- c("gene", "cluster", "p_val_adj", "avg_log2FC")
missing_cols <- setdiff(required_cols, colnames(markers))
if (length(missing_cols) > 0) {
  stop(glue("Missing required columns in input: {paste(missing_cols, collapse=\", \")}"))
}

# 数据类型转换和列选择
markers <- markers %>%
  transmute(
    gene = as.character(.data$gene),
    cluster = as.character(.data$cluster),
    p_val_adj = as.numeric(.data$p_val_adj),
    avg_log2FC = as.numeric(.data$avg_log2FC)
  )

# 根据阈值过滤标记基因
markers_flt <- markers %>%
  filter(!is.na(gene), gene != "", !is.na(p_val_adj), !is.na(avg_log2FC)) %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0)

# 排除指定的cluster
exclude_vec <- unique(strsplit(opt$exclude, ",")[[1]])
exclude_vec <- trimws(exclude_vec)
exclude_vec <- exclude_vec[nzchar(exclude_vec)]
if (length(exclude_vec) > 0) {
  markers_flt <- markers_flt %>% filter(!(.data$cluster %in% exclude_vec))
}

# 检查过滤后是否还有数据
if (nrow(markers_flt) == 0) {
  stop("No markers remain after filtering and exclusion.")
}

# 获取要分析的cluster列表
clusters <- sort(unique(markers_flt$cluster))
cat(glue("[INFO] Clusters to analyze: {paste(clusters, collapse=\", \")}") %>% paste0("\n") )

# 按cluster分组基因列表
cluster_gene_sets <- markers_flt %>%
  distinct(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

# -----------------------
# ID mapping SYMBOL -> ENTREZ
# -----------------------
cat("[INFO] Mapping SYMBOL to ENTREZID...\n")
# 安全的基因ID映射函数
bitr_safely <- function(symbols) {
  if (length(symbols) == 0) return(tibble(SYMBOL = character(), ENTREZID = character()))
  res <- tryCatch({
    # 使用clusterProfiler进行基因ID转换
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>% as_tibble()
  }, error = function(e) {
    message("bitr error: ", conditionMessage(e))
    tibble(SYMBOL = character(), ENTREZID = character())
  })
  # 去重并过滤无效映射
  distinct(res, SYMBOL, ENTREZID) %>% filter(!is.na(ENTREZID), nzchar(ENTREZID))
}

# 对每个cluster进行基因ID映射
mapped_list <- lapply(seq_len(nrow(cluster_gene_sets)), function(i) {
  cl <- cluster_gene_sets$cluster[i]
  syms <- unlist(cluster_gene_sets$genes[i])
  mapdf <- bitr_safely(syms)
  n_in <- length(unique(syms)); n_map <- nrow(mapdf)
  cat(glue("[MAP] cluster={cl} input_symbols={n_in} mapped_ENTREZ={n_map} rate={sprintf('%.1f%%', 100*n_map/max(1,n_in))}\n"))
  list(cluster = cl, map = mapdf)
})

# 构建背景基因集：所有cluster映射成功的ENTREZID并集
bg_entrez <- mapped_list %>% lapply(function(x) x$map$ENTREZID) %>% unlist(use.names = FALSE) %>% unique()
cat(glue("[INFO] Background universe size (ENTREZ): {length(bg_entrez)}\n"))

# -----------------------
# Enrichment wrappers
# -----------------------
# GO富集分析包装函数
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

# KEGG富集分析：优先使用enrichKEGG，失败时使用KEGGREST + enricher
prepare_kegg_term2gene <- function() {
  cat("[INFO] Preparing KEGG TERM2GENE via KEGGREST...\n")
  # 基因到通路的映射
  link <- KEGGREST::keggLink("pathway", "mmu") # 返回命名字符向量：names如"mmu:1234"，values如"path:mmuXXXXX"
  gene_ids <- sub("^mmu:", "", names(link))
  pathway_ids <- sub("^path:", "", as.character(link))
  term2gene <- tibble(term = pathway_ids, gene = gene_ids)
  # 获取通路名称
  term2name <- KEGGREST::keggList("pathway", "mmu")
  term2name <- tibble(term = sub("^path:", "", names(term2name)), name = as.character(term2name))
  list(T2G = term2gene, T2N = term2name)
}

# KEGG富集分析包装函数
run_enrichKEGG_mmu <- function(entrez, t2g_bundle = NULL) {
  if (length(entrez) < 1) return(NULL)
  # 首先尝试使用enrichKEGG
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
  # 备用方案：使用enricher + TERM2GENE
  if (is.null(t2g_bundle)) t2g_bundle <- prepare_kegg_term2gene()
  t2g <- t2g_bundle$T2G; t2n <- t2g_bundle$T2N
  # 过滤到背景空间以保持universe一致
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
# 将GeneRatio格式（如"3/45"）转换为数值
ratio_to_numeric <- function(x) {
  ifelse(grepl("/", x), {
    sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
  }, as.numeric(x))
}

# 安全保存数据框为CSV
save_df_csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0) {
    readr::write_csv(tibble(), path)
  } else {
    readr::write_csv(df, path)
  }
}

# 双格式保存图形（PDF和PNG）
save_plot_dual <- function(p, basepath, width = 9, height = 7) {
  ggplot2::ggsave(filename = paste0(basepath, ".pdf"), plot = p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(filename = paste0(basepath, ".png"), plot = p, width = width, height = height, dpi = 300)
}

# 文本换行包装函数
wrap_terms <- function(x, width = 50) stringr::str_wrap(x, width = width)

# 构建过滤关键词的正则表达式
filter_keywords <- unique(trimws(unlist(strsplit(opt$filter_terms, ","))))
filter_keywords <- filter_keywords[nzchar(filter_keywords)]
filter_regex <- if (length(filter_keywords) > 0) paste(filter_keywords, collapse = "|") else "^$NEVERMATCH$"

# 初始化汇总容器
summary_go <- list()
summary_kegg <- list()

# 为备用路径准备KEGG bundle
kegg_bundle <- NULL

# -----------------------
# Per-cluster enrichment and outputs
# -----------------------
# 对每个cluster进行富集分析
for (cl in clusters) {
  # 清理cluster名称用于文件路径
  cl_sanit <- gsub("[^A-Za-z0-9_-]+", "_", cl)
  cat(glue("[INFO] Processing cluster {cl} ...\n"))
  out_cluster_dir <- file.path(outdir, paste0("cluster_", cl_sanit))
  if (!dir.exists(out_cluster_dir)) dir.create(out_cluster_dir, recursive = TRUE, showWarnings = FALSE)

  # 获取当前cluster的ENTREZID
  entrez <- mapped_list[[which(clusters == cl)]]$map$ENTREZID %>% unique()

  # GO富集分析
  ego <- run_enrichGO_bp(entrez)
  ego_df <- if (!is.null(ego)) as.data.frame(ego) else tibble()
  go_full_path <- file.path(out_cluster_dir, "GO_BP_full.csv")
  save_df_csv(ego_df, go_full_path)

  # 处理GO富集结果
  ego_top <- tibble()
  if (nrow(ego_df) > 0) {
    # 过滤并选择top N结果
    ego_top <- ego_df %>%
      filter(!str_detect(Description, regex(filter_regex, ignore_case = TRUE))) %>%
      arrange(p.adjust, pvalue) %>%
      slice_head(n = opt$top)
    go_top_path <- file.path(out_cluster_dir, "GO_BP_top20_filtered.csv")
    save_df_csv(ego_top, go_top_path)

    # 为过滤后的结果创建点图（从data.frame构建以避免S4方法问题）
    if (nrow(ego_top) > 0) {
      ego_top2 <- ego_top %>%
        mutate(GeneRatioNum = ratio_to_numeric(GeneRatio)) %>%
        mutate(Description_wrapped = wrap_terms(Description, width = 50)) %>%
        arrange(p.adjust, pvalue)
      # 创建ggplot点图
      p_go <- ggplot(ego_top2, aes(x = GeneRatioNum, y = forcats::fct_reorder(Description_wrapped, GeneRatioNum))) +
        geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.8) +
        scale_color_viridis_c(name = "-log10(p.adjust)") +
        scale_size_continuous(name = "GeneRatio", range = c(1.5, 8)) +
        labs(x = "GeneRatio", y = "Term", title = glue("Cluster {cl} GO {opt$go_ont} (Top{opt$top} filtered)")) +
        theme_bw(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
      save_plot_dual(p_go, file.path(out_cluster_dir, "GO_BP_dotplot_top20_filtered"), width = 9, height = 7)
    }

    # 添加到汇总容器
    if (nrow(ego_top) > 0) {
      summary_go[[cl]] <- ego_top %>%
        transmute(cluster = cl, ID, Description, GeneRatio = ratio_to_numeric(GeneRatio), p.adjust)
    }
  } else {
    cat(glue("[WARN] Cluster {cl} has no significant GO terms.\n"))
    # 仍然发出空的top文件
    save_df_csv(ego_top, file.path(out_cluster_dir, "GO_BP_top20_filtered.csv"))
  }

  # KEGG富集分析
  eke <- run_enrichKEGG_mmu(entrez, t2g_bundle = kegg_bundle)
  if (is.null(eke) && is.null(kegg_bundle)) {
    # 如果尚未准备，再次尝试准备
    kegg_bundle <- prepare_kegg_term2gene()
    eke <- run_enrichKEGG_mmu(entrez, t2g_bundle = kegg_bundle)
  }
  eke_df <- if (!is.null(eke)) as.data.frame(eke) else tibble()
  kegg_full_path <- file.path(out_cluster_dir, "KEGG_full.csv")
  save_df_csv(eke_df, kegg_full_path)

  # 处理KEGG富集结果
  eke_top <- tibble()
  if (nrow(eke_df) > 0) {
    # 过滤并选择top N结果
    eke_top <- eke_df %>%
      filter(!str_detect(Description, regex(filter_regex, ignore_case = TRUE))) %>%
      arrange(p.adjust, pvalue) %>%
      slice_head(n = opt$top)
    kegg_top_path <- file.path(out_cluster_dir, "KEGG_top20_filtered.csv")
    save_df_csv(eke_top, kegg_top_path)

    # 为过滤后的结果创建点图（从data.frame构建以避免S4方法问题）
    if (nrow(eke_top) > 0) {
      eke_top2 <- eke_top %>%
        mutate(GeneRatioNum = ratio_to_numeric(GeneRatio)) %>%
        mutate(Description_wrapped = wrap_terms(Description, width = 50)) %>%
        arrange(p.adjust, pvalue)
      # 创建ggplot点图
      p_kegg <- ggplot(eke_top2, aes(x = GeneRatioNum, y = forcats::fct_reorder(Description_wrapped, GeneRatioNum))) +
        geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.8) +
        scale_color_viridis_c(name = "-log10(p.adjust)") +
        scale_size_continuous(name = "GeneRatio", range = c(1.5, 8)) +
        labs(x = "GeneRatio", y = "Pathway", title = glue("Cluster {cl} KEGG (Top{opt$top} filtered)")) +
        theme_bw(base_size = 11) +
        theme(axis.text.y = element_text(size = 8))
      save_plot_dual(p_kegg, file.path(out_cluster_dir, "KEGG_dotplot_top20_filtered"), width = 9, height = 7)
    }

    # 添加到汇总容器
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
# 创建气泡图函数
make_bubble <- function(df, title) {
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(term = Description)
  # 按平均-log10(p.adjust)排序terms
  df <- df %>% group_by(term) %>% mutate(score = -log10(p.adjust)) %>%
    mutate(avg_score = mean(score, na.rm = TRUE)) %>% ungroup()
  term_levels <- df %>% distinct(term, avg_score) %>% arrange(desc(avg_score)) %>% pull(term)
  df <- df %>% mutate(term = factor(wrap_terms(term, width = 50), levels = wrap_terms(term_levels, width = 50)))

  # 创建气泡图
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

# GO汇总
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

# KEGG汇总
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
# 创建README文件
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

# 关闭日志重定向
sink(type = "message")
sink()
close(log_con)