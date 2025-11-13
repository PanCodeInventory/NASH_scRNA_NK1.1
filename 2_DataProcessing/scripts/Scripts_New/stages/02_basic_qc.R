#!/usr/bin/env Rscript

#' 阶段2：基础质量控制
#' 
#' 功能：
#' 1. 双胞检测（scDblFinder）
#' 2. 双胞去除
#' 3. 质量控制报告生成
#' 
#' 输入：阶段1输出的Seurat对象
#' 输出：去除双胞后的Seurat对象、双胞检测结果图件、QC报告

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(SingleCellExperiment)
  library(scDblFinder)
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
  input_path <- if (length(args) > 2) args[3] else "2_DataProcessing/RDS/nk.integrated.imported.rds"
  
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
  
  # 步骤2：双胞检测
  message("[2/4] 双胞检测...")
  
  # 确定样本分组列
  sample_col <- NULL
  if ("orig.ident" %in% colnames(obj@meta.data) && length(unique(obj$orig.ident)) > 1) {
    sample_col <- "orig.ident"
  } else if ("timepoint" %in% colnames(obj@meta.data) && length(unique(obj$timepoint)) > 1) {
    sample_col <- "timepoint"
  } else if ("sample_name" %in% colnames(obj@meta.data) && length(unique(obj$sample_name)) > 1) {
    sample_col <- "sample_name"
  } else {
    message("未找到可用的样本分组列，scDblFinder将不分组运行")
  }
  
  # 运行双胞检测
  obj_before <- obj  # 保存处理前的对象用于比较
  
  obj <- detect_doublets_scDblFinder(
    obj, 
    sample_col = sample_col,
    expected = config$quality_control$doublet_expected
  )
  
  # 步骤3：生成双胞检测结果可视化
  message("[3/4] 生成双胞检测结果可视化...")
  
  # 生成双胞检测图
  doublet_plots <- generate_doublet_plots(
    obj, 
    score_col = "scDblFinder.score",
    class_col = "scDblFinder.class"
  )
  
  # 保存双胞检测图
  for (plot_name in names(doublet_plots)) {
    filename <- sprintf("Doublet_%s", plot_name)
    save_publication_plots(
      doublet_plots[[plot_name]],
      filename = filename,
      output_dir = paths_config$output$plots_dir,
      width = 8, height = 6,
      formats = config$visualization$formats
    )
  }
  
  # 步骤4：去除双胞并保存结果
  message("[4/4] 去除双胞并保存结果...")
  
  # 统计双胞检测前的细胞数
  n_before <- ncol(obj)
  
  # 去除双胞
  obj_filtered <- subset(obj, subset = scDblFinder.class == "singlet")
  n_after <- ncol(obj_filtered)
  n_removed <- n_before - n_after
  removal_percent <- 100 * n_removed / n_before
  
  message(sprintf("双胞去除: %d -> %d 细胞 (移除 %d, %.2f%%)", 
                 n_before, n_after, n_removed, removal_percent))
  
  # 保存去除双胞后的对象
  output_rds_path <- file.path(paths_config$output$rds_dir, "nk.integrated.qc_filtered.rds")
  saveRDS(obj_filtered, file = output_rds_path)
  message(sprintf("保存QC过滤后的对象到: %s", output_rds_path))
  
  # 生成QC报告
  report_path <- file.path(paths_config$output$reports_dir, "basic_qc_report.md")
  
  # 比较处理前后的数据
  comparison_results <- compare_before_after(obj_before, obj_filtered, step_name = "doublet_removal")
  
  # 创建报告内容
  report_lines <- c(
    "# 基础质量控制报告",
    "",
    sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## 处理参数",
    sprintf("- 双胞检测方法: scDblFinder"),
    sprintf("- 期望双胞比例: %.2f", config$quality_control$doublet_expected),
    sprintf("- 样本分组列: %s", ifelse(is.null(sample_col), "无", sample_col)),
    "",
    "## 输入数据",
    sprintf("- 输入文件: %s", basename(input_path)),
    sprintf("- 处理前细胞数: %d", n_before),
    sprintf("- 处理前基因数: %d", nrow(obj_before)),
    "",
    "## 双胞检测结果",
    sprintf("- 双胞分类统计:"),
    "  ```",
    paste0("  ", capture.output(print(table(obj$scDblFinder.class))), collapse = "\n"),
    "  ```",
    sprintf("- 双胞评分范围: %.4f - %.4f", 
            range(obj$scDblFinder.score)[1], range(obj$scDblFinder.score)[2]),
    sprintf("- 双胞评分中位数: %.4f", median(obj$scDblFinder.score)),
    "",
    "## 过滤结果",
    sprintf("- 过滤后细胞数: %d", n_after),
    sprintf("- 移除细胞数: %d", n_removed),
    sprintf("- 移除比例: %.2f%%", removal_percent),
    "",
    "## 数据完整性验证",
    sprintf("- 细胞验证: %s", validation_results$cells$status),
    sprintf("- 基因验证: %s", validation_results$genes$status),
    sprintf("- Assay验证: %s", validation_results$assays$status),
    sprintf("- 元数据验证: %s", validation_results$metadata$status),
    "",
    "## 输出文件",
    sprintf("- QC过滤后对象: %s", basename(output_rds_path)),
    sprintf("- 双胞检测图件: %s", paths_config$output$plots_dir),
    "",
    "## 下一步",
    "运行 03_cell_annotation.R 进行细胞类型注释"
  )
  
  writeLines(report_lines, report_path)
  message(sprintf("QC报告保存到: %s", report_path))
  
  # 打印摘要
  message("=== 基础质量控制完成 ===")
  message(sprintf("处理前细胞数: %d", n_before))
  message(sprintf("处理后细胞数: %d", n_after))
  message(sprintf("移除双胞数: %d (%.2f%%)", n_removed, removal_percent))
  message(sprintf("检测到的双胞比例: %.2f%%", 100 * sum(obj$scDblFinder.class == "doublet") / n_before))
  
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