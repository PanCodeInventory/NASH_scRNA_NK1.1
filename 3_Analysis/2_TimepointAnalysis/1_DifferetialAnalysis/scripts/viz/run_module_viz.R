#!/usr/bin/env Rscript

# NASH NK module visualization script
# - Volcano (per contrast)
# - MA plot (per contrast)
# - Module heatmap across contrasts
# - Enrichment bar/dot for selected pathways (per contrast)

suppressPackageStartupMessages({
  req <- c(
    "yaml","readr","dplyr","tibble","stringr","ggplot2","ggrepel",
    "scales","purrr","tidyr","glue","fs","pheatmap","ggtext"
  )
  missing <- req[!vapply(req, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    msg <- paste0(
      "缺少R包: ", paste(missing, collapse = ", "),
      "\n请先安装，例如: install.packages(c(\"",
      paste(missing, collapse = "\",\""),
      "\"))\n然后重试: Rscript scripts/viz/run_module_viz.R"
    )
    stop(msg)
  }
})

library(yaml)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(scales)
library(purrr)
library(tidyr)
library(glue)
library(fs)
library(pheatmap)
library(ggtext)

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg
}

ensure_dir <- function(path) {
  if (!fs::dir_exists(path)) fs::dir_create(path, recurse = TRUE)
}

# Load DEG table and standardize column names
load_deg_for_contrast <- function(contrast, cfg) {
  base <- cfg$input$tables_dir
  candidates <- c(cfg$input$deg_filename, cfg$input$deg_alt_filename)
  paths <- file.path(base, contrast, candidates)
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop(glue("未找到DEG文件: {paste(paths, collapse=' , ')}"))
  }
  path <- existing[1]
  message(glue("读取DEG: {path}"))
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  # Standardize
  nm <- names(df)
  # Accept common variants
  if ("log2FoldChange" %in% nm) df <- dplyr::rename(df, log2FC = log2FoldChange)
  if (!"log2FC" %in% names(df)) stop("DE表缺少 log2FC/log2FoldChange")
  if (!"padj" %in% names(df)) {
    # try qval or adjusted_pvalue variants
    cand <- c("qval","qvalue","padj")
    hit <- cand[cand %in% names(df)][1]
    if (!is.na(hit)) df <- dplyr::rename(df, padj = !!sym(hit)) else stop("DE表缺少 padj/qval")
  }
  if (!"gene" %in% names(df)) {
    alt <- c("Gene","symbol","SYMBOL","gene_name")
    hit <- alt[alt %in% names(df)][1]
    if (!is.na(hit)) df <- dplyr::rename(df, gene = !!sym(hit)) else stop("DE表缺少 gene 列")
  }
  if (!"baseMean" %in% names(df)) {
    # allow baseMean missing; create pseudo if absent
    df$baseMean <- NA_real_
  }
  df <- df %>% mutate(contrast = contrast)
  df
}

# Flatten module definitions
get_modules <- function(cfg) {
  modules <- cfg$modules
  # Build a tidy frame for genes and colors
  mod_list <- purrr::imap(modules, function(m, name) {
    tibble(
      module_id = name,
      label = m$label %||% name,
      color = m$color %||% "#888888",
      gene = m$genes %||% character()
    )
  }) %>% bind_rows()
  # unique genes per module
  mod_list <- mod_list %>% distinct(module_id, gene, .keep_all = TRUE)
  # pathways (flatten)
  pathway_tbl <- purrr::imap_dfr(modules, function(m, name) {
    p1 <- m$pathways %||% character()
    p2 <- character()
    det <- m$details
    if (!is.null(det)) {
      # recursive extraction of $pathways in nested lists
      extract_paths <- function(x) {
        out <- character()
        if (is.list(x)) {
          for (k in names(x)) {
            if (k == "pathways" && !is.null(x[[k]])) {
              out <- c(out, as.character(x[[k]]))
            } else {
              out <- c(out, extract_paths(x[[k]]))
            }
          }
        }
        out
      }
      p2 <- extract_paths(det)
    }
    tibble(module_id = name, pathway = unique(c(p1, p2)))
  }) %>% filter(!is.na(pathway), pathway != "")

  list(genes = mod_list, pathways = pathway_tbl)
}

# Volcano plot
plot_volcano <- function(df, mod_genes, mod_colors, thresholds, outpath, y_max_override = NULL) {
  d <- df %>% mutate(
    neglog10padj = -log10(padj),
    gene_upper = toupper(gene),
    module = dplyr::case_when(
      gene %in% mod_genes[[1]] ~ names(mod_genes)[1],
      gene %in% mod_genes[[2]] ~ names(mod_genes)[2],
      TRUE ~ "Other"
    )
  )
  cols <- c("Other" = "#BFBFBF")
  cols[names(mod_colors)] <- unname(mod_colors)

  # 计算用于绘图的Y轴上限，并将超过上限的点与标签钉在顶部
  if (!is.null(y_max_override)) {
    y_max <- y_max_override
  } else {
    y_max <- suppressWarnings(max(d$neglog10padj, na.rm = TRUE))
    if (!is.finite(y_max)) y_max <- 5
    y_max <- y_max + 1.0
  }
  d <- d %>% mutate(neglog10padj_plot = pmin(neglog10padj, y_max - 0.5))

  p <- ggplot(d, aes(x = log2FC, y = neglog10padj_plot, color = module)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = cols) +
    theme_minimal(base_size = 12) +
    labs(x = "log2FC", y = "-log10(padj)", color = "Module") +
    geom_vline(xintercept = c(-thresholds$log2fc, thresholds$log2fc), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(thresholds$padj), linetype = "dashed", color = "grey50")

  # Label all specially annotated genes (module genes) with bold black text
  label_genes <- unique(unlist(mod_genes))
  lab_df <- d %>% filter(gene %in% label_genes)
  if (nrow(lab_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      aes(y = neglog10padj_plot, label = gene),
      size = 3.2, max.overlaps = 500, fontface = "bold", color = "black", show.legend = FALSE
    )
  }
  # 保持X轴左右范围一致（以|log2FC|的最大值对称），并提升Y轴上限以容纳全部标签
  max_abs <- suppressWarnings(max(abs(d$log2FC), na.rm = TRUE))
  if (!is.finite(max_abs)) max_abs <- thresholds$log2fc
  max_abs <- max(max_abs, thresholds$log2fc + 0.1)
  # 使用已计算的 y_max，并设置坐标范围
  p <- p + coord_cartesian(xlim = c(-max_abs, max_abs), ylim = c(0, y_max))
  # 输出画布比例横向 2:1，确保纵向空间充足
  ggsave(outpath, p, width = 12, height = 6, dpi = 300)
}

# MA plot
plot_ma <- function(df, mod_genes, mod_colors, thresholds, outpath) {
  d <- df %>% mutate(
    A = log10(baseMean + 1),
    module = dplyr::case_when(
      gene %in% mod_genes[[1]] ~ names(mod_genes)[1],
      gene %in% mod_genes[[2]] ~ names(mod_genes)[2],
      TRUE ~ "Other"
    )
  )
  cols <- c("Other" = "#BFBFBF")
  cols[names(mod_colors)] <- unname(mod_colors)

  p <- ggplot(d, aes(x = A, y = log2FC, color = module)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = cols) +
    theme_minimal(base_size = 12) +
    labs(x = "log10(baseMean+1)", y = "log2FC", color = "Module") +
    geom_hline(yintercept = c(-thresholds$log2fc, thresholds$log2fc), linetype = "dashed", color = "grey50")
  ggsave(outpath, p, width = 6.5, height = 5.2, dpi = 300)
}

# Heatmap across contrasts
plot_heatmap <- function(deg_list, mod_genes_all, mod_info, outpath) {
  # Build matrix genes x contrasts (robust)
  contrasts <- names(deg_list)
  all_genes <- unique(mod_genes_all)
  if (length(all_genes) == 0) {
    message("模块基因集为空，跳过热图")
    return(invisible(NULL))
  }
  cols_list <- lapply(contrasts, function(ct){
    d <- deg_list[[ct]] %>% select(gene, log2FC)
    v <- setNames(d$log2FC, d$gene)
    vals <- v[all_genes]
    as.numeric(vals)
  })
  mat <- do.call(cbind, cols_list)
  colnames(mat) <- contrasts
  rownames(mat) <- all_genes
  # skip if all NA
  if (all(is.na(mat))) {
    message("所有模块基因在对比中均缺失log2FC，跳过热图")
    return(invisible(NULL))
  }
  # Row annotation for modules
  ann <- mod_info %>%
    filter(gene %in% all_genes) %>%
    group_by(gene) %>%
    summarise(module_label = dplyr::first(label), .groups = "drop") %>%
    right_join(tibble(gene = all_genes), by = "gene") %>%
    mutate(module_label = ifelse(is.na(module_label), "Other", module_label))
  row_ann <- data.frame(Module = ann$module_label, row.names = ann$gene, check.names = FALSE)
  # Colors
  breaks <- seq(-2, 2, length.out = 101)
  palette <- colorRampPalette(c("#3B4CC0","white","#B40426"))(length(breaks)-1)
  # filter rows that are all NA and align annotation
  keep <- rowSums(is.na(mat)) < ncol(mat)
  if (!all(keep)) {
    mat <- mat[keep, , drop = FALSE]
    row_ann <- row_ann[rownames(mat), , drop = FALSE]
  }
  pheatmap::pheatmap(
    mat,
    color = palette,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    annotation_row = row_ann,
    fontsize_row = 6,
    fontsize_col = 10,
    filename = outpath,
    width = 6.5,
    height = max(4, min(10, 0.2 * nrow(mat) + 2))
  )
}

# Enrichment plot (bar)
plot_enrichment <- function(contrast, cfg, mod_paths_tbl, outpath) {
  base <- cfg$input$tables_dir
  f <- file.path(base, contrast, cfg$input$gsea_filename)
  if (!file.exists(f)) {
    message(glue("GSEA结果缺失: {f}, 跳过"))
    return(invisible(NULL))
  }
  eg <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  # standardize column names
  nm <- names(eg)
  if (!"pathway" %in% nm) {
    cand <- c("term","Pathway","pathway")
    hit <- cand[cand %in% nm][1]
    if (!is.na(hit)) eg <- dplyr::rename(eg, pathway = !!sym(hit)) else stop("GSEA表缺少 pathway 列")
  }
  if (!"NES" %in% names(eg)) {
    cand <- c("nes","NES")
    hit <- cand[cand %in% names(eg)][1]
    if (!is.na(hit)) eg <- dplyr::rename(eg, NES = !!sym(hit)) else stop("GSEA表缺少 NES 列")
  }
  if (!"padj" %in% names(eg)) {
    cand <- c("padj","qval","qvalue")
    hit <- cand[cand %in% names(eg)][1]
    if (!is.na(hit)) eg <- dplyr::rename(eg, padj = !!sym(hit)) else stop("GSEA表缺少 padj/qval 列")
  }
  # 排序依据（padj升序或NES降序）
  rank_by <- cfg$enrichment$rank_by %||% "padj"
  if (tolower(rank_by) == "nes") {
    egs <- eg %>% arrange(desc(NES), padj)
  } else {
    egs <- eg %>% arrange(padj, desc(NES))
  }
  # 焦点通路 + 自定义颜色通路
  focus_terms <- cfg$enrichment$focus_pathways %||% character()
  custom_colors <- cfg$enrichment$custom_colors %||% list()
  colors_vec <- if (length(custom_colors) > 0) unlist(custom_colors) else character()
  focus_terms <- unique(c(focus_terms, names(colors_vec)))
  # 计算焦点通路在排序后的名次，取最大名次作为top数量
  pos <- match(focus_terms, egs$pathway)
  pos <- pos[!is.na(pos)]
  k <- if (length(pos) > 0) max(pos) else min(10L, nrow(egs))
  top_terms <- head(egs$pathway, k)
  sel <- unique(c(top_terms, focus_terms))
  egf <- egs %>% filter(pathway %in% sel) %>% mutate(label_raw = stringr::str_replace(pathway, "HALLMARK_", ""))
  if (nrow(egf) == 0) {
    message(glue("{contrast}: 未找到焦点或Top通路，跳过绘图"))
    return(invisible(NULL))
  }
  # 构造富文本标签：焦点/自定义颜色通路加粗并着色
  focus_color <- cfg$enrichment$focus_color %||% "#D62728"
  custom_map <- if (length(colors_vec) > 0) tibble(pathway = names(colors_vec), custom_color = unname(colors_vec)) else tibble(pathway = character(), custom_color = character())
  egf <- egf %>%
    left_join(custom_map, by = "pathway") %>%
    mutate(
      color_use = ifelse(!is.na(custom_color), custom_color, ifelse(pathway %in% focus_terms, focus_color, NA_character_)),
      display_label = ifelse(!is.na(color_use), glue("<span style='color:{color_use}'><b>{label_raw}</b></span>"), label_raw)
    )
  # 保持当前排序的因子水平（翻转坐标时自上而下顺序）
  egf$display_label <- factor(egf$display_label, levels = rev(egf$display_label))
  p <- ggplot(egf, aes(x = display_label, y = NES, fill = NES > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2C7BB6", "FALSE" = "#D7191C"), guide = "none") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = ggtext::element_markdown()) +
    labs(x = NULL, y = "NES", title = contrast)
  ggsave(outpath, p, width = 6.5, height = max(3.5, 0.35 * nrow(egf)), dpi = 300)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

main <- function() {
  # config path default
  cfg_path <- "3_Analysis/2_DifferetialAnalysis/scripts/viz/modules.yaml"
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) cfg_path <- args[1]
  cfg <- read_config(cfg_path)

  mod <- get_modules(cfg)
  mod_info <- mod$genes # module_id, label, color, gene
  # two primary modules assumed
  modules_order <- unique(mod_info$module_id)
  if (length(modules_order) < 2) stop("配置中模块数少于2，请检查 modules.yaml")

  # per-module colors and gene sets (only first two modules used for colored highlight)
  mod_colors <- mod_info %>% group_by(module_id) %>% summarise(color = dplyr::first(color)) %>% deframe()
  mod_colors <- mod_colors[modules_order]
  mod_genes_named <- lapply(modules_order[1:2], function(mid){
    mod_info %>% filter(module_id == mid) %>% pull(gene)
  })
  names(mod_genes_named) <- modules_order[1:2]

  thresholds <- list(padj = cfg$thresholds$padj %||% 0.05,
                     log2fc = cfg$thresholds$log2fc %||% 0.5,
                     top_n = cfg$labeling$top_n_per_module %||% 10)

  out_base <- cfg$output$base_dir %||% "3_Analysis/2_DifferetialAnalysis/results/plots/module_viz"
  dirs <- list(
    volcano = file.path(out_base, "volcano"),
    ma = file.path(out_base, "ma"),
    heatmap = file.path(out_base, "heatmap"),
    enrichment = file.path(out_base, "enrichment")
  )
  purrr::walk(dirs, ensure_dir)

  contrasts <- cfg$contrasts$include %||% character()
  if (length(contrasts) == 0) stop("配置 contrasts.include 为空")

  # load all deg
  deg_list <- purrr::map(contrasts, ~load_deg_for_contrast(.x, cfg))
  names(deg_list) <- contrasts

  # plots per contrast
  purrr::walk(contrasts, function(ct){
    d <- deg_list[[ct]]
    # volcano
    out_v <- file.path(dirs$volcano, glue("{ct}_module_volcano.png"))
    plot_volcano(d, mod_genes_named, mod_colors, thresholds, out_v, cfg$volcano$y_max %||% NULL)
    message(glue("保存: {out_v}"))
    # MA
    out_m <- file.path(dirs$ma, glue("{ct}_module_ma.png"))
    plot_ma(d, mod_genes_named, mod_colors, thresholds, out_m)
    message(glue("保存: {out_m}"))
    # enrichment per contrast
    out_e <- file.path(dirs$enrichment, glue("{ct}_module_enrichment.png"))
    plot_enrichment(ct, cfg, mod$pathways, out_e)
    if (file.exists(out_e)) message(glue("保存: {out_e}"))
  })

  # heatmap across contrasts
  all_module_genes <- unique(unlist(mod_genes_named))
  out_h <- file.path(dirs$heatmap, "modules_log2FC_heatmap.png")
  plot_heatmap(deg_list, all_module_genes, mod_info, out_h)
  message(glue("保存: {out_h}"))

  message("完成。")
}

if (identical(environment(), globalenv())) {
  tryCatch(main(), error = function(e){
    message("运行失败: ", e$message)
    quit(status = 1)
  })
}
