#!/usr/bin/env Rscript

#' 阶段1：数据导入与预处理
#' 
#' 功能：
#' 1. 读取10x数据
#' 2. 创建Seurat对象
#' 3. SCTransform归一化
#' 4. 数据整合
#' 5. 基础降维
#' 6. 初始可视化
#' 
#' 输入：四个时间点的10x数据目录
#' 输出：整合后的Seurat对象、基础可视化图件、数据质量报告

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(future)
})

# 加载工具函数
source("utils/seurat_utils.R")
source("utils/plotting_utils.R")
source("utils/validation_utils.R")

# 设置全局选项
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)  # 8GB
future::plan("sequential")

#' 主函数
main <- function() {
  # 解析命令行参数
  args <- commandArgs(trailingOnly = TRUE)
  
  # 默认配置文件路径
  config_path <- if (length(args) > 0) args[1] else "config/parameters.yaml"
  paths_path <- if (length(args) > 1) args[2] else "config/paths.yaml"
  
  # 加载配置
  message("加载配置文件...")
  config <- load_config(config_path)
  paths_config <- load_config(paths_path)
  
  # 设置随机种子
  set.seed(config$computation$seed)
  
  # 创建输出目录
  output_dirs <- c(paths_config$output$rds_dir, 
                   paths_config$output$plots_dir, 
                   paths_config$output$reports_dir)
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # 步骤1：读取10x数据并创建Seurat对象
  message("[1/6] 读取10x数据并创建Seurat对象...")
  obj_list <- list()
  
  for (i in seq_along(paths_config$input$sample_dirs)) {
    sample_dir <- paths_config$input$sample_dirs[i]
    timepoint <- paths_config$timepoints[i]
    sample_name <- paths_config$sample_mapping[[timepoint]]
    
    message(sprintf("  处理样本: %s (%s)", sample_name, timepoint))
    
    # 检查目录存在性
    if (!dir.exists(sample_dir)) {
      stop(sprintf("数据目录不存在: %s", sample_dir))
    }
    
    # 安全读取10x数据
    counts <- read_10x_data_safe(sample_dir)
    
    # 创建Seurat对象并添加元数据
    obj <- create_seurat_with_metadata(counts, 
                                      project = "NK", 
                                      timepoint = timepoint)
    
    # 添加样本名
    obj$sample_name <- sample_name
    
    obj_list[[timepoint]] <- obj
  }
  
  # 步骤2：SCTransform归一化
  message("[2/6] SCTransform归一化...")
  obj_list <- lapply(obj_list, function(obj) {
    normalize_with_sctransform(obj, vst_flavor = "v2")
  })
  
  # 步骤3：数据整合
  message("[3/6] 数据整合...")
  nfeatures <- config$preprocessing$nfeatures
  integrated_obj <- integrate_seurat_objects(obj_list, nfeatures = nfeatures)
  
  # 添加时间点因子
  integrated_obj$timepoint <- factor(integrated_obj$timepoint, 
                                    levels = paths_config$timepoints)
  
  # 步骤4：基础降维
  message("[4/6] 基础降维...")
  dims <- config$umap$n_neighbors  # 使用配置中的邻居数作为默认维度
  resolution <- 0.5  # 默认分辨率
  
  integrated_obj <- run_dimensionality_reduction(integrated_obj, 
                                               dims = dims, 
                                               resolution = resolution)
  
  # 步骤5：生成基础可视化
  message("[5/6] 生成基础可视化...")
  
  # ElbowPlot
  elbow_plot <- ElbowPlot(integrated_obj, ndims = dims)
  save_publication_plots(elbow_plot, 
                        filename = "ElbowPlot",
                        output_dir = paths_config$output$plots_dir,
                        width = 8, height = 6,
                        formats = config$visualization$formats)
  
  # UMAP总览图
  overview_plot <- DimPlot(
    integrated_obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    label.size = config$visualization$label_size,
    repel = TRUE
  ) + ggtitle("All NK Cells Clustered")
  
  save_publication_plots(overview_plot,
                        filename = "UMAP_Overview",
                        output_dir = paths_config$output$plots_dir,
                        width = config$visualization$plot_width,
                        height = config$visualization$plot_height,
                        formats = config$visualization$formats)
  
  # 按时间点分面的UMAP图
  timepoint_plot <- plot_umap_by_timepoint(
    integrated_obj,
    group_by = "seurat_clusters",
    split_by = "timepoint",
    label_size = config$visualization$label_size
  )
  
  save_publication_plots(timepoint_plot,
                        filename = "UMAP_by_Timepoint",
                        output_dir = paths_config$output$plots_dir,
                        width = config$visualization$plot_width,
                        height = config$visualization$plot_height,
                        formats = config$visualization$formats)
  
  # 步骤6：保存结果和生成报告
  message("[6/6] 保存结果和生成报告...")
  
  # 保存整合后的Seurat对象
  output_rds_path <- file.path(paths_config$output$rds_dir,
                               paste0("nk.integrated.imported.rds"))
  saveRDS(integrated_obj, file = output_rds_path)
  message(sprintf("保存整合对象到: %s", output_rds_path))
  
  # 验证对象
  validation_results <- validate_seurat_object(integrated_obj)
  
  # 生成数据质量报告
  report_path <- file.path(paths_config$output$reports_dir, 
                           "data_import_report.md")
  
  # 创建报告内容
  report_lines <- c(
    "# 数据导入与预处理报告",
    "",
    sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## 处理参数",
    sprintf("- 可变基因数: %d", nfeatures),
    sprintf("- 归一化方法: %s", config$preprocessing$normalization_method),
    sprintf("- 整合方法: %s", config$preprocessing$integration_method),
    sprintf("- 主成分数: %d", dims),
    sprintf("- 聚类分辨率: %.2f", resolution),
    "",
    "## 输入数据",
    sprintf("- 样本数: %d", length(paths_config$input$sample_dirs)),
    sprintf("- 时间点: %s", paste(paths_config$timepoints, collapse = ", ")),
    "",
    "## 处理结果",
    sprintf("- 总细胞数: %d", ncol(integrated_obj)),
    sprintf("- 总基因数: %d", nrow(integrated_obj)),
    sprintf("- 簇数: %d", length(levels(integrated_obj$seurat_clusters))),
    "",
    "## 验证结果",
    sprintf("- 细胞验证: %s", validation_results$cells$status),
    sprintf("- 基因验证: %s", validation_results$genes$status),
    sprintf("- Assay验证: %s", validation_results$assays$status),
    sprintf("- 元数据验证: %s", validation_results$metadata$status),
    "",
    "## 输出文件",
    sprintf("- Seurat对象: %s", basename(output_rds_path)),
    sprintf("- 可视化图件: %s", paths_config$output$plots_dir),
    "",
    "## 下一步",
    "运行 02_basic_qc.R 进行质量控制"
  )
  
  writeLines(report_lines, report_path)
  message(sprintf("数据质量报告保存到: %s", report_path))
  
  # 打印摘要
  message("=== 数据导入与预处理完成 ===")
  message(sprintf("处理样本数: %d", length(obj_list)))
  message(sprintf("总细胞数: %d", ncol(integrated_obj)))
  message(sprintf("总基因数: %d", nrow(integrated_obj)))
  message(sprintf("识别簇数: %d", length(levels(integrated_obj$seurat_clusters))))
  
  return(integrated_obj)
}

# 错误处理
tryCatch({
  result <- main()
  message("脚本执行成功!")
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})