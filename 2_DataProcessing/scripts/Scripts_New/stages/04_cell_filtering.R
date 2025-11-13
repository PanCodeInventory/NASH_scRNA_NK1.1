#!/usr/bin/env Rscript

#' 阶段3.2：特定细胞类型去除
#' 
#' 功能：
#' 1. 细胞类型过滤
#' 2. 规则配置
#' 3. 过滤验证
#' 
#' 输入：阶段3.1输出的Seurat对象
#' 输出：过滤后的Seurat对象、细胞类型组成变化图、过滤效果报告

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(patchwork)
})

# 加载工具函数
source("utils/seurat_utils.R")
source("utils/plotting_utils.R")
source("utils/validation_utils.R")

#' 主函数
main <- function() {
  # 解析命令行参数
  args <- commandArgs(trailingOnly = TRUE)
  
  # 默认配置文件路径
  config_path <- if (length(args) > 0) args[1] else "config/parameters.yaml"
  paths_path <- if (length(args) > 1) args[2] else "config/paths.yaml"
  cluster_mapping_path <- if (length(args) > 2) args[3] else "config/cluster_mapping.yaml"
  input_path <- if (length(args) > 3) args[4] else "2_DataProcessing/RDS/nk.integrated.annotated.rds"
  
  # 加载配置
  message("加载配置文件...")
  config <- load_config(config_path)
  paths_config <- load_config(paths_path)
  cluster_config <- load_config(cluster_mapping_path)
  
  # 设置随机种子
  set.seed(config$computation$seed)
  
  # 创建输出目录
  output_dirs <- c(paths_config$output$rds_dir, 
                   paths_config$output$plots_dir, 
                   paths_config$output$reports_dir,
                   paths_config$output$data_dir)
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # 步骤1：加载输入数据
  message("[1/4] 加载输入数据...")
  if (!file.exists(input_path)) {
    stop(sprintf("输入文件不存在: %s", input_path))
  }
  
  obj <- readRDS(input_path)
  message(sprintf("成功加载Seurat对象: %d个细胞, %d个基因", ncol(obj), nrow(obj)))
  
  # 验证输入对象
  validation_results <- validate_seurat_object(obj)
  if (validation_results$cells$status == "FAIL") {
    stop("输入Seurat对象验证失败: 无细胞数据")
  }
  
  # 步骤2：细胞类型过滤
  message("[2/4] 细胞类型过滤...")
  
  # 选择标签列
  label_col <- if ("singleR.pruned" %in% colnames(obj@meta.data)) {
    "singleR.pruned"
  } else if ("singleR.label" %in% colnames(obj@meta.data)) {
    "singleR.label"
  } else {
    stop("未找到SingleR注释列 (singleR.pruned/singleR.label)")
  }
  
  # 获取过滤配置
  remove_types <- config$cell_filtering$remove_cell_types
  keep_types <- config$cell_filtering$keep_cell_types
  filter_by_score <- config$cell_filtering$filter_by_score
  score_threshold <- config$cell_filtering$score_threshold
  
  message(sprintf("过滤配置: 去除类型=%s, 保留类型=%s", 
                 paste(remove_types, collapse = ", "), 
                 paste(keep_types, collapse = ", ")))
  
  # 保存处理前的对象用于比较
  obj_before <- obj
  
  # 应用细胞类型过滤
  obj_filtered <- filter_cells_by_type(
    obj, 
    remove_types = remove_types,
    keep_types = keep_types,
    label_col = label_col
  )
  
  # 步骤3：生成过滤效果可视化
  message("[3/4] 生成过滤效果可视化...")
  
  # 计算过滤前后的细胞类型统计
  labels_before <- obj_before@meta.data[[label_col]]
  labels_after <- obj_filtered@meta.data[[label_col]]
  
  # 过滤前统计
  stats_before <- data.frame(
    label = names(table(labels_before)),
    count_before = as.numeric(table(labels_before)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(percent_before = 100 * count_before / sum(count_before))
  
  # 过滤后统计
  stats_after <- data.frame(
    label = names(table(labels_after)),
    count_after = as.numeric(table(labels_after)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(percent_after = 100 * count_after / sum(count_after))
  
  # 合并统计
  comparison_df <- dplyr::full_join(stats_before, stats_after, by = "label") %>%
    dplyr::mutate(
      count_before = ifelse(is.na(count_before), 0, count_before),
      count_after = ifelse(is.na(count_after), 0, count_after),
      percent_before = ifelse(is.na(percent_before), 0, percent_before),
      percent_after = ifelse(is.na(percent_after), 0, percent_after),
      removed = count_before - count_after,
      removed_percent = 100 * (count_before - count_after) / count_before
    ) %>%
    dplyr::arrange(dplyr::desc(count_before))
  
  # 保存比较统计
  comparison_path <- file.path(paths_config$output$data_dir, "cell_type_filtering_comparison.csv")
  write.csv(comparison_df, comparison_path, row.names = FALSE)
  
  # 生成过滤前后对比图
  p_before <- ggplot(comparison_df, aes(x = reorder(label, -count_before), y = count_before)) +
    geom_col(fill = "#6baed6") +
    coord_flip() +
    theme_bw() +
    labs(title = sprintf("Cell Types Before Filtering (n=%d)", ncol(obj_before)), 
         x = "Cell Type", y = "Cell Count")
  
  p_after <- ggplot(comparison_df, aes(x = reorder(label, -count_after), y = count_after)) +
    geom_col(fill = "#3182bd") +
    coord_flip() +
    theme_bw() +
    labs(title = sprintf("Cell Types After Filtering (n=%d)", ncol(obj_filtered)), 
         x = "Cell Type", y = "Cell Count")
  
  # 组合对比图
  combined_plot <- p_before + p_after + 
    plot_annotation(title = "Cell Type Filtering: Before vs After")
  
  save_publication_plots(
    combined_plot,
    filename = "Cell_Type_Filtering_Comparison",
    output_dir = paths_config$output$plots_dir,
    width = 12, height = 10,
    formats = config$visualization$formats
  )
  
  # 步骤4：保存结果和生成报告
  message("[4/4] 保存结果和生成报告...")
  
  # 保存过滤后的对象
  output_rds_path <- file.path(paths_config$output$rds_dir, "nk.integrated.filtered.rds")
  saveRDS(obj_filtered, file = output_rds_path)
  message(sprintf("保存过滤后的对象到: %s", output_rds_path))
  
  # 比较处理前后的数据
  comparison_results <- compare_before_after(obj_before, obj_filtered, step_name = "cell_type_filtering")
  
  # 生成过滤报告
  report_path <- file.path(paths_config$output$reports_dir, "cell_filtering_report.md")
  
  # 创建报告内容
  report_lines <- c(
    "# 细胞类型过滤报告",
    "",
    sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## 处理参数",
    sprintf("- 使用的标签列: %s", label_col),
    sprintf("- 去除的细胞类型: %s", paste(remove_types, collapse = ", ")),
    sprintf("- 保留的细胞类型: %s", paste(keep_types, collapse = ", ")),
    sprintf("- 基于评分过滤: %s", if (filter_by_score) "是" else "否"),
    if (filter_by_score) sprintf("- 评分阈值: %.2f", score_threshold) else "",
    "",
    "## 输入数据",
    sprintf("- 输入文件: %s", basename(input_path)),
    sprintf("- 过滤前细胞数: %d", ncol(obj_before)),
    sprintf("- 过滤前基因数: %d", nrow(obj_before)),
    sprintf("- 过滤前细胞类型数: %d", length(unique(labels_before))),
    "",
    "## 过滤结果",
    sprintf("- 过滤后细胞数: %d", ncol(obj_filtered)),
    sprintf("- 移除细胞数: %d", ncol(obj_before) - ncol(obj_filtered)),
    sprintf("- 移除比例: %.2f%%", 100 * (ncol(obj_before) - ncol(obj_filtered)) / ncol(obj_before)),
    sprintf("- 过滤后细胞类型数: %d", length(unique(labels_after))),
    "",
    "## 细胞类型变化（前10个）",
    "```",
    paste(utils::capture.output(print(utils::head(comparison_df, 10))), collapse = "\n"),
    "```",
    "",
    "## 数据完整性验证",
    sprintf("- 细胞验证: %s", validation_results$cells$status),
    sprintf("- 基因验证: %s", validation_results$genes$status),
    sprintf("- Assay验证: %s", validation_results$assays$status),
    sprintf("- 元数据验证: %s", validation_results$metadata$status),
    "",
    "## 输出文件",
    sprintf("- 过滤后对象: %s", basename(output_rds_path)),
    sprintf("- 过滤比较数据: %s", basename(comparison_path)),
    sprintf("- 过滤效果图件: %s", paths_config$output$plots_dir),
    "",
    "## 下一步",
    "运行 05_parameter_tuning.R 进行参数优化"
  )
  
  writeLines(report_lines, report_path)
  message(sprintf("过滤报告保存到: %s", report_path))
  
  # 打印摘要
  message("=== 细胞类型过滤完成 ===")
  message(sprintf("过滤前细胞数: %d", ncol(obj_before)))
  message(sprintf("过滤后细胞数: %d", ncol(obj_filtered)))
  message(sprintf("移除细胞数: %d (%.2f%%)", 
                 ncol(obj_before) - ncol(obj_filtered),
                 100 * (ncol(obj_before) - ncol(obj_filtered)) / ncol(obj_before)))
  message(sprintf("过滤前细胞类型数: %d", length(unique(labels_before))))
  message(sprintf("过滤后细胞类型数: %d", length(unique(labels_after))))
  
  return(obj_filtered)
}

# 错误处理
tryCatch({
  result <- main()
  message("脚本执行成功!")
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})