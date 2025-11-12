#!/usr/bin/env Rscript

#' 阶段5：可视化
#' 
#' 功能：
#' 1. 最终UMAP生成
#' 2. 多维度可视化
#' 3. 组合图生成
#' 4. 多格式输出
#' 
#' 输入：阶段4输出的Seurat对象和最佳参数
#' 输出：最终UMAP图件、簇组成趋势图、堆叠柱状图、组合展示图、可视化报告

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(patchwork)
  library(tidyr)
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
  input_path <- if (length(args) > 2) args[3] else "2_DataProcessing/RDS/nk.integrated.tuned.rds"
  
  # 加载配置
  message("加载配置文件...")
  config <- load_config(config_path)
  paths_config <- load_config(paths_path)
  
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
  message("[1/6] 加载输入数据...")
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
  
  # 确保timepoint因子顺序
  if ("timepoint" %in% colnames(obj@meta.data)) {
    obj$timepoint <- factor(obj$timepoint, levels = paths_config$timepoints)
  }
  
  # 选择标签列
  label_col <- if ("singleR.pruned" %in% colnames(obj@meta.data)) {
    "singleR.pruned"
  } else if ("singleR.label" %in% colnames(obj@meta.data)) {
    "singleR.label"
  } else {
    NULL
  }
  
  # 步骤2：确保降维结果存在
  message("[2/6] 确保降维结果存在...")
  
  # 确保PCA存在
  if (!"pca" %in% names(obj@reductions)) {
    message("未发现PCA，运行PCA...")
    obj <- RunPCA(obj, verbose = FALSE)
  }
  
  # 确保UMAP存在
  if (!"umap" %in% names(obj@reductions)) {
    message("未发现UMAP，运行UMAP...")
    # 获取当前使用的维度数
    dims <- 30  # 默认值
    if ("pca" %in% names(obj@reductions)) {
      dims <- ncol(Embeddings(obj, "pca"))
    }
    obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:dims, verbose = FALSE)
  }
  
  # 步骤3：生成基础UMAP可视化
  message("[3/6] 生成基础UMAP可视化...")
  
  # 按簇着色的UMAP总览图
  overview_plot <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    label.size = config$visualization$label_size,
    repel = TRUE
  ) + ggtitle("Final NK Cells Clustering")
  
  save_publication_plots(
    overview_plot,
    filename = "UMAP_Final_Overview",
    output_dir = paths_config$output$plots_dir,
    width = config$visualization$plot_width,
    height = config$visualization$plot_height,
    formats = config$visualization$formats
  )
  
  # 按时间点分面的UMAP图
  timepoint_plot <- plot_umap_by_timepoint(
    obj,
    group_by = "seurat_clusters",
    split_by = "timepoint",
    label_size = config$visualization$label_size
  )
  
  save_publication_plots(
    timepoint_plot,
    filename = "UMAP_Final_by_Timepoint",
    output_dir = paths_config$output$plots_dir,
    width = config$visualization$plot_width,
    height = config$visualization$plot_height,
    formats = config$visualization$formats
  )
  
  # 按细胞类型着色的UMAP图（如果有SingleR注释）
  if (!is.null(label_col)) {
    celltype_plot <- DimPlot(
      obj,
      reduction = "umap",
      group.by = label_col,
      label = FALSE
    ) + ggtitle(paste("Final NK Cells colored by", label_col))
    
    save_publication_plots(
      celltype_plot,
      filename = paste0("UMAP_Final_by_", label_col),
      output_dir = paths_config$output$plots_dir,
      width = config$visualization$plot_width,
      height = config$visualization$plot_height,
      formats = config$visualization$formats
    )
  }
  
  # 步骤4：生成按时间点的独立UMAP图
  message("[4/6] 生成按时间点的独立UMAP图...")
  
  for (tp in paths_config$timepoints) {
    if (!tp %in% levels(obj$timepoint)) next
    
    obj_tp <- subset(obj, subset = timepoint == tp)
    p_tp <- DimPlot(
      obj_tp,
      reduction = "umap",
      group.by = "seurat_clusters",
      label = TRUE,
      repel = TRUE
    ) + ggtitle(paste("NK Cells UMAP -", tp))
    
    tp_safe <- gsub("[^A-Za-z0-9_-]", "_", tp)
    save_publication_plots(
      p_tp,
      filename = paste0("UMAP_Final_Clusters_", tp_safe),
      output_dir = paths_config$output$plots_dir,
      width = 10, height = 8,
      formats = config$visualization$formats
    )
  }
  
  # 步骤5：生成簇组成分析图
  message("[5/6] 生成簇组成分析图...")
  
  # 计算簇组成数据
  composition_df <- obj@meta.data %>%
    dplyr::select(seurat_clusters, timepoint) %>%
    dplyr::mutate(
      seurat_clusters = as.factor(seurat_clusters),
      timepoint = factor(timepoint, levels = paths_config$timepoints)
    ) %>%
    dplyr::count(timepoint, seurat_clusters, name = "n") %>%
    dplyr::group_by(timepoint) %>%
    dplyr::mutate(percent = 100 * n / sum(n)) %>%
    dplyr::ungroup()
  
  # 保存簇组成数据
  composition_path <- file.path(paths_config$output$data_dir, "cluster_composition_by_timepoint.csv")
  write.csv(composition_df, composition_path, row.names = FALSE)
  
  # 簇组成趋势图
  trend_plot <- plot_cluster_composition_trends(
    composition_df,
    timepoint_col = "timepoint",
    cluster_col = "seurat_clusters",
    count_col = "n"
  )
  
  save_publication_plots(
    trend_plot,
    filename = "Cluster_Composition_Trends",
    output_dir = paths_config$output$plots_dir,
    width = 12, height = 8,
    formats = config$visualization$formats
  )
  
  # 堆叠柱状图
  stacked_plot <- plot_stacked_bar(
    composition_df,
    timepoint_col = "timepoint",
    cluster_col = "seurat_clusters",
    count_col = "n",
    position = "fill"
  )
  
  save_publication_plots(
    stacked_plot,
    filename = "Cluster_Composition_Stacked",
    output_dir = paths_config$output$plots_dir,
    width = 12, height = 8,
    formats = config$visualization$formats
  )
  
  # 步骤6：生成组合图和报告
  message("[6/6] 生成组合图和报告...")
  
  # 创建组合图
  plot_list <- list(overview_plot, timepoint_plot)
  if (!is.null(label_col)) {
    plot_list[[3]] <- celltype_plot
  }
  
  combined_plot <- create_combined_plot(
    plot_list[1:2],  # 只组合前两个图避免过于复杂
    layout = "vertical",
    heights = c(1, 1)
  )
  
  save_publication_plots(
    combined_plot,
    filename = "UMAP_Final_Combined",
    output_dir = paths_config$output$plots_dir,
    width = 12, height = 16,
    formats = config$visualization$formats
  )
  
  # 保存最终对象
  final_rds_path <- file.path(paths_config$output$rds_dir, "nk.integrated.final.rds")
  saveRDS(obj, file = final_rds_path)
  message(sprintf("保存最终对象到: %s", final_rds_path))
  
  # 生成可视化报告
  report_path <- file.path(paths_config$output$reports_dir, "visualization_report.md")
  
  # 创建报告内容
  report_lines <- c(
    "# 可视化报告",
    "",
    sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## 处理参数",
    sprintf("- 图件宽度: %d", config$visualization$plot_width),
    sprintf("- 图件高度: %d", config$visualization$plot_height),
    sprintf("- 分辨率: %d", config$visualization$dpi),
    sprintf("- 输出格式: %s", paste(config$visualization$formats, collapse = ", ")),
    sprintf("- 配色方案: %s", config$visualization$color_palette),
    "",
    "## 输入数据",
    sprintf("- 输入文件: %s", basename(input_path)),
    sprintf("- 细胞数: %d", ncol(obj)),
    sprintf("- 基因数: %d", nrow(obj)),
    sprintf("- 簇数: %d", length(levels(obj$seurat_clusters))),
    sprintf("- 时间点数: %d", length(levels(obj$timepoint))),
    "",
    "## 生成的图件",
    "- UMAP总览图 (按簇着色)",
    "- UMAP分面图 (按时间点分面)",
    if (!is.null(label_col)) paste("- UMAP细胞类型图 (按", label_col, "着色)") else "- UMAP细胞类型图 (无注释数据)",
    "- 按时间点的独立UMAP图",
    "- 簇组成趋势图",
    "- 簇组成堆叠柱状图",
    "- 组合展示图",
    "",
    "## 数据完整性验证",
    sprintf("- 细胞验证: %s", validation_results$cells$status),
    sprintf("- 基因验证: %s", validation_results$genes$status),
    sprintf("- Assay验证: %s", validation_results$assays$status),
    sprintf("- 元数据验证: %s", validation_results$metadata$status),
    "",
    "## 输出文件",
    sprintf("- 最终对象: %s", basename(final_rds_path)),
    sprintf("- 簇组成数据: %s", basename(composition_path)),
    sprintf("- 可视化图件: %s", paths_config$output$plots_dir),
    "",
    "## 分析摘要",
    sprintf("- 主要簇: %s (%.1f%%)", 
            levels(obj$seurat_clusters)[which.max(table(obj$seurat_clusters))],
            100 * max(table(obj$seurat_clusters)) / length(obj$seurat_clusters)),
    sprintf("- 最丰富时间点: %s", levels(obj$timepoint)[which.max(table(obj$timepoint))]),
    if (!is.null(label_col)) sprintf("- 主要细胞类型: %s", 
                                   levels(obj@meta.data[[label_col]])[which.max(table(obj@meta.data[[label_col]]))]) else "- 主要细胞类型: 无数据",
    "",
    "## 流程完成",
    "所有核心阶段脚本已完成，生成的数据可用于下游分析。"
  )
  
  writeLines(report_lines, report_path)
  message(sprintf("可视化报告保存到: %s", report_path))
  
  # 打印摘要
  message("=== 可视化完成 ===")
  message(sprintf("可视化细胞数: %d", ncol(obj)))
  message(sprintf("识别簇数: %d", length(levels(obj$seurat_clusters))))
  message(sprintf("时间点数: %d", length(levels(obj$timepoint))))
  message(sprintf("生成图件数: %d", length(list.files(paths_config$output$plots_dir, pattern = "\\.(png|pdf|svg)$"))))
  
  return(obj)
}

# 错误处理
tryCatch({
  result <- main()
  message("脚本执行成功!")
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})