#!/usr/bin/env Rscript

#' 阶段4：参数优化
#' 
#' 功能：
#' 1. 网格搜索（dims×resolution）
#' 2. 多指标评估
#' 3. 候选选择
#' 4. 热力图可视化
#' 
#' 输入：阶段3.2输出的Seurat对象
#' 输出：最佳参数组合、参数评估热力图、候选参数列表、调优报告

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(cluster)
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
  input_path <- if (length(args) > 2) args[3] else "2_DataProcessing/RDS/nk.integrated.annotated.rds"
  
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
  message("[1/5] 加载输入数据...")
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
  
  # 步骤2：确保PCA存在
  message("[2/5] 确保PCA存在...")
  
  # 确保PCA存在且至少覆盖最大dims
  max_dims <- max(config$clustering$dims_range)
  ensure_pca <- function(obj, target_dims = max_dims) {
    has_pca <- "pca" %in% names(obj@reductions)
    if (has_pca) {
      emb <- tryCatch(Embeddings(obj, "pca"), error = function(e) NULL)
      if (!is.null(emb) && ncol(emb) >= target_dims) {
        message(sprintf("检测到现有PCA (%d PCs)，跳过重算", ncol(emb)))
        return(obj)
      }
    }
    
    # 选择用于降维的assay
    assays_avail <- Seurat::Assays(obj)
    dr_assay <- if ("SCT" %in% assays_avail) "SCT" else "RNA"
    if (dr_assay == "SCT") {
      vf <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
      if (length(vf) == 0) dr_assay <- "RNA"
    }
    DefaultAssay(obj) <- dr_assay
    
    if (dr_assay == "RNA") {
      message("PCA保障：使用RNA路线")
      obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      if (length(VariableFeatures(obj)) == 0) {
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      }
      feats <- VariableFeatures(obj)
    } else {
      message("PCA保障：使用SCT路线")
      feats <- VariableFeatures(obj)
    }
    
    obj <- ScaleData(obj, features = feats, verbose = FALSE)
    obj <- RunPCA(obj, features = feats, npcs = max(50, target_dims), verbose = FALSE)
    return(obj)
  }
  
  obj <- ensure_pca(obj, target_dims = max_dims)
  
  # 步骤3：参数调优网格搜索
  message("[3/5] 参数调优网格搜索...")
  
  dims_grid <- config$clustering$dims_range
  res_grid <- seq(config$clustering$resolution_range[1], 
                  config$clustering$resolution_range[2], 
                  by = config$clustering$resolution_step)
  
  message(sprintf("调参范围: dims=%s; res=[%.1f..%.1f] by %.1f",
                 paste(dims_grid, collapse = ","), 
                 min(res_grid), max(res_grid), config$clustering$resolution_step))
  
  # 使用工具函数进行参数调优
  metrics_df <- tune_clustering_parameters(obj, dims_grid = dims_grid, res_grid = res_grid)
  
  # 保存指标表
  metrics_path <- file.path(paths_config$output$data_dir, "parameter_tuning_metrics.csv")
  write.csv(metrics_df, metrics_path, row.names = FALSE)
  message(sprintf("参数调优指标保存到: %s", metrics_path))
  
  # 步骤4：生成热力图可视化
  message("[4/5] 生成热力图可视化...")
  
  # 轮廓系数热力图
  silhouette_heatmap <- plot_parameter_heatmap(
    metrics_df,
    x_var = "res",
    y_var = "dims", 
    fill_var = "median_silhouette",
    title = "Median Silhouette across Dims × Resolution"
  )
  
  save_publication_plots(
    silhouette_heatmap,
    filename = "Parameter_Tuning_Silhouette_Heatmap",
    output_dir = paths_config$output$plots_dir,
    width = 10, height = 6,
    formats = config$visualization$formats
  )
  
  # 簇数热力图
  clusters_heatmap <- plot_parameter_heatmap(
    metrics_df,
    x_var = "res",
    y_var = "dims", 
    fill_var = "n_clusters",
    title = "Number of Clusters across Dims × Resolution"
  )
  
  save_publication_plots(
    clusters_heatmap,
    filename = "Parameter_Tuning_Clusters_Heatmap",
    output_dir = paths_config$output$plots_dir,
    width = 10, height = 6,
    formats = config$visualization$formats
  )
  
  # 步骤5：选择最佳参数并生成报告
  message("[5/5] 选择最佳参数并生成报告...")
  
  # 选择候选组合
  select_candidates <- function(df, per_dims_top = config$tuning$top_candidates_per_dims, 
                               global_top = config$tuning$global_top_candidates) {
    df2 <- df %>%
      dplyr::filter(!is.na(median_silhouette)) %>%
      dplyr::filter(n_clusters >= config$tuning$silhouette_min_clusters, 
                    n_clusters <= config$tuning$max_clusters) %>%
      dplyr::filter(is.na(min_cluster_frac) | min_cluster_frac >= config$tuning$min_cluster_fraction) %>%
      dplyr::arrange(dplyr::desc(median_silhouette),
                     dplyr::across(pct_silhouette_negative, ~ ifelse(is.na(.x), Inf, .x)),
                     dims)
    
    best_per_dims <- df2 %>%
      dplyr::group_by(dims) %>%
      dplyr::slice_head(n = per_dims_top) %>%
      dplyr::ungroup()
    
    best_global <- df2 %>%
      dplyr::slice_head(n = global_top)
    
    list(best_per_dims = best_per_dims, best_global = best_global)
  }
  
  candidates <- select_candidates(metrics_df)
  
  # 保存候选组合
  best_per_dims_path <- file.path(paths_config$output$data_dir, "best_parameters_per_dims.csv")
  best_global_path <- file.path(paths_config$output$data_dir, "top_parameter_candidates.csv")
  
  write.csv(candidates$best_per_dims, best_per_dims_path, row.names = FALSE)
  write.csv(candidates$best_global, best_global_path, row.names = FALSE)
  
  message(sprintf("候选参数保存到: %s 和 %s", best_per_dims_path, best_global_path))
  
  # 选择最佳参数
  dims_best <- NA_integer_
  res_best <- NA_real_
  
  if (nrow(candidates$best_global) > 0) {
    dims_best <- candidates$best_global$dims[1]
    res_best <- candidates$best_global$res[1]
    message(sprintf("选择最佳参数: dims=%d, res=%.1f", dims_best, res_best))
  } else if (nrow(candidates$best_per_dims) > 0) {
    dims_best <- candidates$best_per_dims$dims[1]
    res_best <- candidates$best_per_dims$res[1]
    message(sprintf("选择最佳参数: dims=%d, res=%.1f", dims_best, res_best))
  } else {
    # 回退策略
    dims_best <- 20
    res_best <- 0.5
    warning(sprintf("未找到有效候选组合，使用回退参数: dims=%d, res=%.1f", dims_best, res_best))
  }
  
  # 使用最佳参数重新计算并保存对象
  obj <- FindNeighbors(obj, dims = 1:dims_best, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:dims_best, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res_best, verbose = FALSE)
  
  # 保存调优后的对象
  output_rds_path <- file.path(paths_config$output$rds_dir, "nk.integrated.tuned.rds")
  saveRDS(obj, file = output_rds_path)
  message(sprintf("保存调优后的对象到: %s", output_rds_path))
  
  # 生成调优报告
  report_path <- file.path(paths_config$output$reports_dir, "parameter_tuning_report.md")
  
  # 创建报告内容
  report_lines <- c(
    "# 参数优化报告",
    "",
    sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## 处理参数",
    sprintf("- 主成分候选: %s", paste(dims_grid, collapse = ", ")),
    sprintf("- 分辨率范围: %.1f - %.1f (步长: %.1f)", 
            min(res_grid), max(res_grid), config$clustering$resolution_step),
    sprintf("- 最小簇数: %d", config$tuning$silhouette_min_clusters),
    sprintf("- 最大簇数: %d", config$tuning$max_clusters),
    sprintf("- 最小簇比例: %.3f", config$tuning$min_cluster_fraction),
    "",
    "## 输入数据",
    sprintf("- 输入文件: %s", basename(input_path)),
    sprintf("- 细胞数: %d", ncol(obj)),
    sprintf("- 基因数: %d", nrow(obj)),
    "",
    "## 调优结果",
    sprintf("- 评估的组合数: %d", nrow(metrics_df)),
    sprintf("- 最佳参数: dims=%d, res=%.1f", dims_best, res_best),
    sprintf("- 最佳轮廓系数: %.4f", 
            if (nrow(candidates$best_global) > 0) candidates$best_global$median_silhouette[1] else NA),
    sprintf("- 最佳簇数: %d", 
            if (nrow(candidates$best_global) > 0) candidates$best_global$n_clusters[1] else NA),
    "",
    "## 前5个候选参数",
    "```",
    paste(utils::capture.output(print(utils::head(candidates$best_global, 5))), collapse = "\n"),
    "```",
    "",
    "## 数据完整性验证",
    sprintf("- 细胞验证: %s", validation_results$cells$status),
    sprintf("- 基因验证: %s", validation_results$genes$status),
    sprintf("- Assay验证: %s", validation_results$assays$status),
    sprintf("- 元数据验证: %s", validation_results$metadata$status),
    "",
    "## 输出文件",
    sprintf("- 调优后对象: %s", basename(output_rds_path)),
    sprintf("- 调优指标: %s", basename(metrics_path)),
    sprintf("- 候选参数: %s, %s", basename(best_per_dims_path), basename(best_global_path)),
    sprintf("- 热力图件: %s", paths_config$output$plots_dir),
    "",
    "## 下一步",
    "运行 06_visualization.R 生成最终可视化"
  )
  
  writeLines(report_lines, report_path)
  message(sprintf("调优报告保存到: %s", report_path))
  
  # 打印摘要
  message("=== 参数优化完成 ===")
  message(sprintf("评估参数组合数: %d", nrow(metrics_df)))
  message(sprintf("选择最佳参数: dims=%d, res=%.1f", dims_best, res_best))
  message(sprintf("最终簇数: %d", length(levels(obj$seurat_clusters))))
  
  return(list(
    tuned_object = obj,
    best_params = list(dims = dims_best, resolution = res_best),
    metrics = metrics_df,
    candidates = candidates
  ))
}

# 错误处理
tryCatch({
  result <- main()
  message("脚本执行成功!")
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})