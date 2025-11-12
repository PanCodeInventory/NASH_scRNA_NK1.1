#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# 输入/输出路径
input_rds_path <- "2_DataProcessing/RDS/nk.integrated.v1.rds"
output_dir <- "2_DataProcessing/RDS"
output_filename <- "nk.integrated.v4.rds"

main <- function() {
  message("=== 筛选(NK cells/ILC) + 清理(metadata) 开始 ===")

  # 1) 读取对象
  if (!file.exists(input_rds_path)) {
    stop(sprintf("RDS文件不存在: %s", input_rds_path))
  }
  obj <- readRDS(input_rds_path)
  if (!inherits(obj, "Seurat")) stop("读取对象非 Seurat 类型")
  message(sprintf("原始对象: 细胞=%d, 基因=%d", ncol(obj), nrow(obj)))
  message(sprintf("DefaultAssay(原): %s", DefaultAssay(obj)))

  # 2) 选择 SingleR 列并筛选细胞
  md <- obj@meta.data
  singleR_col <- if ("singleR.label" %in% colnames(md)) {
    "singleR.label"
  } else if ("singleR.pruned" %in% colnames(md)) {
    "singleR.pruned"
  } else {
    stop("未找到 SingleR 注释列（singleR.label 或 singleR.pruned）")
  }
  message(sprintf("使用 SingleR 注释列: %s", singleR_col))

  message("原始 SingleR 细胞类型分布：")
  print(table(md[[singleR_col]]))

  keep_types <- c("NK cells", "ILC")
  target_cells <- md %>%
    filter(.data[[singleR_col]] %in% keep_types) %>%
    rownames()

  if (length(target_cells) == 0) {
    stop("错误：未找到符合条件的细胞（NK cells 或 ILC）")
  }

  obj2 <- subset(obj, cells = target_cells)
  message(sprintf("筛选后: 细胞=%d (保留%.1f%%)", ncol(obj2), 100 * ncol(obj2) / ncol(obj)))
  message(sprintf("DefaultAssay(筛选后): %s", DefaultAssay(obj2)))

  # 清理因子水平
  if (is.factor(obj2@meta.data[[singleR_col]])) {
    obj2@meta.data[[singleR_col]] <- droplevels(obj2@meta.data[[singleR_col]])
  }

  # 3) 清理 metadata（删除 orig.ident 与 integrated_snn_res.*）
  md2 <- obj2@meta.data
  res_cols <- grep("^integrated_snn_res\\.", colnames(md2), value = TRUE)
  cols_to_delete <- intersect(c("orig.ident", res_cols), colnames(md2))

  if (length(cols_to_delete) > 0) {
    message(sprintf("删除 metadata 列: %s", paste(cols_to_delete, collapse = ", ")))
    # 关键修复：drop = FALSE 避免退化为向量
    obj2@meta.data <- md2[, !(colnames(md2) %in% cols_to_delete), drop = FALSE]
  } else {
    message("无可删除的 metadata 列")
  }

  # 4) 结构校验与 Idents 处理
  stopifnot(is.data.frame(obj2@meta.data))
  stopifnot(identical(rownames(obj2@meta.data), colnames(obj2)))

  if ("seurat_clusters" %in% colnames(obj2@meta.data)) {
    Idents(obj2) <- factor(obj2@meta.data$seurat_clusters)
    message("已将 Idents 设为 seurat_clusters")
  } else {
    message("未发现 seurat_clusters，保留原 Idents")
  }

  # 5) 保存
  output_path <- file.path(output_dir, output_filename)
  saveRDS(obj2, file = output_path)
  message(sprintf("已保存筛选+清理后的 RDS: %s", output_path))

  # 6) 摘要
  summary_path <- file.path(output_dir, "nk.integrated.v4_filter_clean_summary.txt")
  filtered_summary <- table(obj2@meta.data[[singleR_col]])
  del_lines <- if (length(cols_to_delete) > 0) paste("-", cols_to_delete) else "- 无"
  clean_cols_lines <- if (ncol(obj2@meta.data) > 0) paste("-", colnames(obj2@meta.data)) else "- 无"

  summary_text <- c(
    sprintf("筛选+清理摘要 - %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("输入文件: %s", basename(input_rds_path)),
    sprintf("输出文件: %s", output_filename),
    sprintf("原始细胞数: %d", ncol(obj)),
    sprintf("筛选后细胞数: %d", ncol(obj2)),
    sprintf("保留比例: %.1f%%", 100 * ncol(obj2) / ncol(obj)),
    sprintf("筛选条件: %s 列中值为 'NK cells' 或 'ILC'", singleR_col),
    "",
    "筛选后 SingleR 分布：",
    paste(names(filtered_summary), filtered_summary, sep = ": "),
    "",
    "删除的 metadata 列：",
    del_lines,
    "",
    "清理后 metadata 列：",
    clean_cols_lines
  )
  writeLines(summary_text, summary_path)
  message(sprintf("摘要保存: %s", summary_path))
  message("=== 处理完成 ===")
}

tryCatch({
  main()
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})
