#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(fgsea)
  library(msigdbr)
})

# ------------------------------------------------------------
# 中文内联注释版（在原脚本基础上补充注释）：
# 目的：基于差异分析输出（results/tables 与 results/data），按配置生成每个比较与“整体(all_cells)”及各簇的可视化：
# - 火山图（Volcano）：展示 log2FC 与 -log10(padj)，并标注最显著的若干基因
# - MA 图：展示表达量（baseMean）与 log2FC 的关系
# - 热图：展示 Top 基因的标准化计数（normalized counts），可按行缩放
# - FGSEA 富集（可选）：对 ranked_genes 的评分向量做快速 GSEA（如 HALLMARK 通路）
#
# 重要约定：
# - 差异分析阶段可能输出两类对象：
#   1) "all_cells.*"（整体分析，不分簇）
#   2) "<cluster>.*"（簇分析，按 seurat_clusters 分层）
#   本脚本会统一识别 *.DESeq2.results.csv 文件名中的前缀（包含 all_cells 或具体簇号），并逐一作图与富集。
# - 各种阈值与显示参数来自 scripts/config.yaml 的 viz 与 deg 配置。
# ------------------------------------------------------------

# ---------------------------
# CLI args：允许通过 --config 指定配置文件路径；若未提供则使用默认路径
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- NULL
if (length(args) >= 2) {
  for (i in seq_along(args)) {
    if (args[i] == "--config" && i + 1 <= length(args)) {
      cfg_path <- args[i + 1]
      break
    }
  }
}
if (is.null(cfg_path)) {
  cfg_path <- "3_Analysis/2_DifferetialAnalysis/scripts/config.yaml"
  message(sprintf("No --config provided. Falling back to default: %s", cfg_path))
}

# ---------------------------
# Helpers：目录创建、时间戳、空值合并、安全读取
# ---------------------------
ensure_dir <- function(p) {
  # 若目录不存在则递归创建
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

now_str <- function() format(Sys.time(), "%Y%m%d-%H%M%S")

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_fread <- function(path, ...) {
  # 安全读取：文件不存在或解析失败时返回NULL，避免流程中断
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s", path))
    return(NULL)
  }
  tryCatch({
    data.table::fread(path, ...)
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", path, e$message))
    NULL
  })
}

# ---------------------------
# Load config
# ---------------------------
cfg <- yaml::read_yaml(cfg_path)

out_root <- cfg$paths$out_root
viz_cfg  <- cfg$viz
deg_cfg  <- cfg$deg

# 显著性阈值（用于 sig 标注与图中辅助线）
lfc_threshold <- as.numeric(deg_cfg$lfc_threshold %||% 0.25)
padj_cutoff   <- as.numeric(deg_cfg$p_adj_cutoff %||% 0.05)

# 火山图设置
vol_top_labels <- as.integer(viz_cfg$volcano$top_labels %||% 10)
vol_label_by   <- viz_cfg$volcano$label_by %||% "padj"  # 可选: "padj" 或 "pvalue"

# MA 图设置（本脚本直接使用结果中的 log2FoldChange，不在此阶段重算收缩）
use_shrunk_lfc <- isTRUE(viz_cfg$ma$use_shrunk_lfc %||% TRUE)

# 热图设置
heat_top_n     <- as.integer(viz_cfg$heatmap$top_n %||% 30)
heat_scale     <- viz_cfg$heatmap$scale %||% "row"      # row/none

# 富集设置
enrich_enabled <- isTRUE(viz_cfg$enrichment$enabled %||% FALSE)
enrich_dbs     <- viz_cfg$enrichment$databases %||% c("HALLMARK")
enrich_min     <- as.integer(viz_cfg$enrichment$min_size %||% 10)
enrich_max     <- as.integer(viz_cfg$enrichment$max_size %||% 500)
rank_metric    <- viz_cfg$enrichment$rank_metric %||% "stat"   # 从 ranked_genes.csv 中取用
enrich_species <- viz_cfg$enrichment$species %||% "Homo sapiens"

comparisons <- cfg$comparisons

# ---------------------------
# Prepare dirs and logging
# ---------------------------
ensure_dir(file.path(out_root, "results", "plots", "volcano"))
ensure_dir(file.path(out_root, "results", "plots", "ma"))
ensure_dir(file.path(out_root, "results", "plots", "heatmap"))
ensure_dir(file.path(out_root, "results", "plots", "enrichment"))
ensure_dir(file.path(out_root, "results", "enrichment"))
ensure_dir(file.path(out_root, "logs", "viz"))

log_path <- file.path(out_root, "logs", "viz", "2_run_viz.log")
zz <- file(log_path, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat(sprintf("[INFO] %s | Start VIZ\n", now_str()))
cat(sprintf("[INFO] Config: %s\n", cfg_path))
cat(sprintf("[INFO] out_root: %s\n", out_root))

# ---------------------------
# Load gene sets (once)：当启用富集时，使用 msigdbr 载入指定数据库的基因集
# ---------------------------
msig_list <- list()
if (enrich_enabled) {
  cat("[INFO] Loading gene sets from msigdbr...\n")
  # Default species: Homo sapiens
  gs_df <- NULL
  db_map <- c(HALLMARK = "H")
  for (db in enrich_dbs) {
    catg <- db_map[[db]] %||% db
    x <- msigdbr::msigdbr(species = enrich_species, category = catg)
    gs_df <- dplyr::bind_rows(gs_df, x)
  }
  if (!is.null(gs_df) && nrow(gs_df) > 0) {
    # 拆分为命名列表：names=通路名，值=gene_symbol 向量
    msig_list <- split(gs_df$gene_symbol, gs_df$gs_name)
    cat(sprintf("[INFO] Loaded %d gene sets from databases: %s\n", length(msig_list), paste(enrich_dbs, collapse = ", ")))
  } else {
    warning("[WARN] No gene sets loaded; enrichment will be skipped.")
    enrich_enabled <- FALSE
  }
}

# ---------------------------
# Iterate comparisons：对每个比较，分别处理 all_cells 与各簇的结果。
# ---------------------------
for (cmp in comparisons) {
  cmp_id <- cmp$id
  cat(sprintf("\n[INFO] Comparison: %s\n", cmp_id))

  table_dir <- file.path(out_root, "results", "tables", cmp_id)
  data_dir  <- file.path(out_root, "results", "data", cmp_id)

  # 子目录（按比较分别存放输出）
  vol_dir <- file.path(out_root, "results", "plots", "volcano", cmp_id); ensure_dir(vol_dir)
  ma_dir  <- file.path(out_root, "results", "plots", "ma", cmp_id); ensure_dir(ma_dir)
  ht_dir  <- file.path(out_root, "results", "plots", "heatmap", cmp_id); ensure_dir(ht_dir)
  enr_plot_dir <- file.path(out_root, "results", "plots", "enrichment", cmp_id); ensure_dir(enr_plot_dir)
  enr_tbl_dir  <- file.path(out_root, "results", "enrichment", cmp_id); ensure_dir(enr_tbl_dir)

  if (!dir.exists(table_dir)) {
    warning(sprintf("[WARN] tables dir not found for %s: %s", cmp_id, table_dir))
    next
  }

  # 发现可用的“簇/整体”标识：依据 *.DESeq2.results.csv 文件名前缀
  res_files <- list.files(table_dir, pattern = "\\.DESeq2\\.results\\.csv$", full.names = FALSE)
  clusters <- unique(gsub("\\.DESeq2\\.results\\.csv$", "", res_files))
  # clusters 含 "all_cells"（整体）或具体簇名（例如 0, 1, 2）
  if (length(clusters) == 0) {
    warning(sprintf("[WARN] No results csv found for %s", cmp_id))
    next
  }

  for (cl in clusters) {
    cat(sprintf("[INFO]   Cluster: %s\n", cl))
    res_path <- file.path(table_dir, sprintf("%s.DESeq2.results.csv", cl))
    res_df <- safe_fread(res_path)
    if (is.null(res_df) || nrow(res_df) == 0) {
      warning(sprintf("[WARN]   Empty or missing results: %s", res_path))
      next
    }

    # 确保关键列存在：gene, log2FoldChange, padj, pvalue（缺失时仍尽量绘图但会发警告）
    needed_cols <- c("gene", "log2FoldChange", "padj", "pvalue")
    missing_cols <- setdiff(needed_cols, colnames(res_df))
    if (length(missing_cols) > 0) {
      warning(sprintf("[WARN]   Missing columns in results for %s: %s", cl, paste(missing_cols, collapse = ", ")))
    }

    # 若未包含 sig 列，则按配置阈值重算
    if (!"sig" %in% colnames(res_df)) {
      res_df <- res_df %>%
        mutate(sig = ifelse(!is.na(padj) & abs(log2FoldChange) >= lfc_threshold & padj <= padj_cutoff, TRUE, FALSE))
    }

    # ---------------------
    # 火山图：y 轴为 -log10(padj)，并用水平/垂直虚线标示阈值；
    # 顶部标注若干基因（由 label_by 决定排序基准）
    # ---------------------
    res_df <- res_df %>% mutate(nlog10_padj = -log10(padj))

    # 处理 Inf/NA 以避免绘图报错：将非有限值设为最大有限值+1
    finite_mask <- is.finite(res_df$nlog10_padj)
    max_finite <- if (any(finite_mask)) max(res_df$nlog10_padj[finite_mask], na.rm = TRUE) else 0
    res_df$nlog10_padj[!finite_mask] <- max_finite + 1

    label_col <- if (tolower(vol_label_by) == "pvalue") "pvalue" else "padj"
    ord <- order(res_df[[label_col]], na.last = TRUE)          # 较小更显著，排前
    top_lbl_df <- res_df[ord, ][seq_len(min(vol_top_labels, nrow(res_df))), , drop = FALSE]

    p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = nlog10_padj, color = sig)) +
      geom_point(alpha = 0.6, size = 1.2) +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D55E00")) +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey50") +
      ggrepel::geom_text_repel(
        data = top_lbl_df %>% filter(!is.na(gene)),
        aes(label = gene),
        size = 3, max.overlaps = 50, min.segment.length = 0
      ) +
      labs(title = sprintf("Volcano | %s | %s", cmp_id, cl),
           x = "log2FoldChange",
           y = "-log10(padj)") +
      theme_minimal(base_size = 12)

    vol_out <- file.path(vol_dir, sprintf("%s.volcano.png", cl))
    ggsave(vol_out, plot = p_vol, width = 7, height = 5, dpi = 300)
    cat(sprintf("[INFO]     Saved volcano: %s\n", vol_out))

    # ---------------------
    # MA 图：x 轴为 log10(baseMean+1)，y 轴为 log2FoldChange（颜色由 sig 决定）
    # 若缺失 baseMean，则以 |log2FC| 的秩作为替代横轴（尽量给出可视化）
    # ---------------------
    if (!"baseMean" %in% colnames(res_df)) {
      warning("[WARN]   'baseMean' not found; MA plot x-axis will use rank of mean instead.")
      res_df$baseMean <- rank(abs(res_df$log2FoldChange), ties.method = "average")
    }
    p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
      geom_point(alpha = 0.6, size = 1.2) +
      scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#1B9E77")) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = sprintf("MA | %s | %s", cmp_id, cl),
           x = "log10(baseMean + 1)",
           y = "log2FoldChange") +
      theme_minimal(base_size = 12)

    ma_out <- file.path(ma_dir, sprintf("%s.MA.png", cl))
    ggsave(ma_out, plot = p_ma, width = 7, height = 5, dpi = 300)
    cat(sprintf("[INFO]     Saved MA: %s\n", ma_out))

    # ---------------------
    # 热图：选取 Top_n 基因（优先 padj 更小，若 padj 相同则 |lfc| 更大），
    # 使用标准化计数（normalized counts），行可按 heat_scale 进行 z-score 等缩放。
    # ---------------------
    nc_path <- file.path(data_dir, sprintf("%s.normalized_counts.csv", cl))
    ncounts_df <- safe_fread(nc_path)
    if (!is.null(ncounts_df) && nrow(ncounts_df) > 0) {
      res_df_ord <- res_df %>%
        arrange(is.na(padj), padj, desc(abs(log2FoldChange)))
      top_genes <- unique(res_df_ord$gene[!is.na(res_df_ord$gene)])[seq_len(min(heat_top_n, sum(!is.na(res_df_ord$gene))))]
      if (!"gene" %in% colnames(ncounts_df)) {
        warning(sprintf("[WARN]   normalized_counts missing 'gene' column: %s", nc_path))
      } else {
        # 转为数据框并设置行名为 gene，避免 data.table 行名赋值不一致问题
        ncounts_df <- as.data.frame(ncounts_df)
        rn <- ncounts_df$gene
        mat <- as.matrix(ncounts_df[, setdiff(colnames(ncounts_df), "gene"), drop = FALSE])
        # 去除 NA 与重复基因名并对齐
        rownames(mat) <- rn
        valid_rows <- !is.na(rownames(mat)) & !duplicated(rownames(mat))
        mat <- mat[valid_rows, , drop = FALSE]
        keep <- intersect(rownames(mat), top_genes)
        if (length(keep) >= 2) {
          hm_path <- file.path(ht_dir, sprintf("%s.topGenes.heatmap.png", cl))
          pheatmap::pheatmap(mat[keep, , drop = FALSE],
                             scale = heat_scale,
                             show_rownames = TRUE,
                             show_colnames = TRUE,
                             clustering_method = "complete",
                             filename = hm_path,
                             width = 8, height = 8)
          cat(sprintf("[INFO]     Saved heatmap: %s\n", hm_path))
        } else {
          warning(sprintf("[WARN]   Not enough top genes for heatmap in %s", cl))
        }
      }
    } else {
      warning(sprintf("[WARN]   normalized_counts not available: %s", nc_path))
    }

    # ---------------------
    # 富集分析：FGSEA
    # 使用 ranked_genes.csv 中的 rank_metric（默认 stat）作为评分向量，
    # 与 msig_list 路径进行快速 GSEA，输出结果表与 Top 条形图。
    # ---------------------
    if (enrich_enabled && length(msig_list) > 0) {
      rk_path <- file.path(data_dir, sprintf("%s.ranked_genes.csv", cl))
      rk_df <- safe_fread(rk_path)
      if (!is.null(rk_df) && nrow(rk_df) > 0) {
        if (!all(c("gene", rank_metric) %in% colnames(rk_df))) {
          warning(sprintf("[WARN]   ranked_genes missing required columns: %s or %s", "gene", rank_metric))
        } else {
          stats <- rk_df[[rank_metric]]
          names(stats) <- rk_df$gene
          # 清理评分：去除非有限值与 NA，去重保留首次出现
          ok <- is.finite(stats) & !is.na(stats) & !is.na(names(stats))
          stats <- stats[ok]
          stats <- stats[!duplicated(names(stats))]
          # 执行 FGSEA
          fg <- tryCatch({
            fgsea::fgsea(pathways = msig_list, stats = stats,
                         minSize = enrich_min, maxSize = enrich_max)
          }, error = function(e) {
            warning(sprintf("[WARN]   FGSEA failed for %s: %s", cl, e$message))
            NULL
          })
          if (!is.null(fg) && nrow(fg) > 0) {
            # 保存结果表
            enr_tbl_out <- file.path(enr_tbl_dir, sprintf("%s.fgsea.results.csv", cl))
            data.table::fwrite(as.data.frame(fg), enr_tbl_out)
            cat(sprintf("[INFO]     Saved enrichment table: %s\n", enr_tbl_out))
            # 另存一份到 results/tables/<cmp_id>/ 便于统一查阅
            enr_tbl_out_tables <- file.path(out_root, "results", "tables", cmp_id, sprintf("%s.fgsea.results.csv", cl))
            data.table::fwrite(as.data.frame(fg), enr_tbl_out_tables)
            cat(sprintf("[INFO]     Copied enrichment table to tables: %s\n", enr_tbl_out_tables))
            # 绘制 Top 15（按 padj 与 |NES| 排序）
            topn <- fg %>%
              arrange(padj, desc(abs(NES))) %>%
              head(15) %>%
              mutate(label = factor(pathway, levels = rev(pathway)))
            if (nrow(topn) > 0) {
              p_enr <- ggplot(topn, aes(x = label, y = NES, fill = -log10(padj))) +
                geom_col() +
                coord_flip() +
                scale_fill_gradient(low = "#fee8c8", high = "#e34a33") +
                labs(title = sprintf("FGSEA | %s | %s", cmp_id, cl),
                     x = "Pathway", y = "NES", fill = "-log10(padj)") +
                theme_minimal(base_size = 12)
              enr_plot_out <- file.path(enr_plot_dir, sprintf("%s.fgsea.hallmark.png", cl))
              ggsave(enr_plot_out, plot = p_enr, width = 7, height = 5, dpi = 300)
              cat(sprintf("[INFO]     Saved enrichment plot: %s\n", enr_plot_out))
            }
          }
        }
      } else {
        warning(sprintf("[WARN]   ranked_genes not available: %s", rk_path))
      }
    } # end enrichment
  } # end cluster loop
} # end comparison loop

# Session info：记录会话信息以便复现
sess_path <- file.path(out_root, "logs", sprintf("sessionInfo_%s.txt", now_str()))
writeLines(c(capture.output(sessionInfo())), con = sess_path)
cat(sprintf("[INFO] Wrote sessionInfo -> %s\n", sess_path))

cat(sprintf("[INFO] %s | VIZ finished\n", now_str()))

# 关闭日志 sink
sink(type = "message")
sink(type = "output")
close(zz)

invisible(TRUE)
