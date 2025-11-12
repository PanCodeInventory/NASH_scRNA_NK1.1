#!/usr/bin/env Rscript

#' 阶段3.1：细胞标注
#' 
#' 功能：
#' 1. SingleR细胞类型注释
#' 2. 注释验证
#' 3. 细胞类型统计
#' 
#' 输入：阶段2输出的Seurat对象
#' 输出：带注释的Seurat对象、细胞类型分布图、注释质量报告

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
  library(scater)
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
  input_path <- if (length(args) > 2) args[3] else "2_DataProcessing/RDS/nk.integrated.qc_filtered.rds"
  
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
  
  # 步骤2：SingleR细胞类型注释
  message("[2/4] SingleR细胞类型注释...")
  
  # 选择assay（优先RNA，其次SCT，再次integrated）
  assays_avail <- Seurat::Assays(obj)
  chosen_assay <- if ("RNA" %in% assays_avail) "RNA" else if ("SCT" %in% assays_avail) "SCT" else assays_avail[[1]]
  DefaultAssay(obj) <- chosen_assay
  message(sprintf("选用assay: %s", chosen_assay))
  
  # 稳健提取counts数据
  get_counts <- function(obj, assay_name) {
    # 尝试layer接口（Seurat v5）
    m <- tryCatch(SeuratObject::GetAssayData(object = obj[[assay_name]], layer = "counts"), error = function(e) NULL)
    if (!is.null(m) && !is.null(dim(m)) && all(dim(m) > 0)) return(m)
    
    # 回退到slot接口（兼容v4/v5）
    m <- tryCatch(Seurat::GetAssayData(obj[[assay_name]], slot = "counts"), error = function(e) NULL)
    if (!is.null(m) && !is.null(dim(m)) && all(dim(m) > 0)) return(m)
    
    return(NULL)
  }
  
  counts <- get_counts(obj, chosen_assay)
  if (is.null(counts) && chosen_assay != "RNA") counts <- get_counts(obj, "RNA")
  if (is.null(counts) && chosen_assay != "SCT") counts <- get_counts(obj, "SCT")
  if (is.null(counts) && chosen_assay != "integrated") counts <- get_counts(obj, "integrated")
  
  if (is.null(counts)) {
    stop("无法从Seurat对象提取非空的counts数据")
  }
  
  # 与Seurat对象细胞对齐
  cells <- colnames(obj)
  common <- intersect(cells, colnames(counts))
  if (length(common) == 0) stop("counts与Seurat对象细胞名不重叠")
  counts <- counts[, common, drop = FALSE]
  cd <- obj@meta.data[common, , drop = FALSE]
  
  # 构建SCE并计算logcounts
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = cd)
  sce <- scater::logNormCounts(sce)
  message("已使用scater::logNormCounts生成logcounts")
  
  # 加载参考数据集
  ref_dataset_name <- config$singler$ref_dataset
  message(sprintf("加载参考数据集: %s", ref_dataset_name))
  ref <- get(ref_dataset_name)()
  
  # 运行SingleR注释
  pred <- SingleR::SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main,
    assay.type.test = config$singler$assay_type,
    assay.type.ref = config$singler$assay_type
  )
  
  # 写回Seurat对象
  obj <- subset(obj, cells = common)
  obj$singleR.label <- pred$labels
  if (!is.null(pred$pruned.labels)) {
    obj$singleR.pruned <- pred$pruned.labels
  }
  message("SingleR注释已写入Seurat对象")
  
  # 步骤3：生成细胞类型统计和可视化
  message("[3/4] 生成细胞类型统计和可视化...")
  
  # 选择使用的标签列
  label_col <- if ("singleR.pruned" %in% colnames(obj@meta.data)) {
    "singleR.pruned"
  } else {
    "singleR.label"
  }
  
  # 生成细胞类型统计
  cell_type_stats <- obj@meta.data %>%
    dplyr::count(!!sym(label_col), name = "n") %>%
    dplyr::mutate(percent = 100 * n / sum(n)) %>%
    dplyr::arrange(desc(n))
  
  # 保存细胞类型统计
  stats_path <- file.path(paths_config$output$data_dir, "cell_type_statistics.csv")
  if (!dir.exists(paths_config$output$data_dir)) {
    dir.create(paths_config$output$data_dir, recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(cell_type_stats, stats_path, row.names = FALSE)
  message(sprintf("细胞类型统计保存到: %s", stats_path))
  
  # 生成细胞类型比例图
  cell_type_plot <- plot_cell_type_proportions(
    cell_type_stats,
    label_col = label_col,
    count_col = "n",
    top_n = 15
  )
  
  save_publication_plots(
    cell_type_plot,
    filename = "Cell_Type_Proportions",
    output_dir = paths_config$output$plots_dir,
    width = config$visualization$plot_width,
    height = config$visualization$plot_height,
    formats = config$visualization$formats
  )
  
  # 步骤4：保存结果和生成报告
  message("[4/4] 保存结果和生成报告...")
  
  # 保存注释后的对象
  output_rds_path <- file.path(paths_config$output$rds_dir, "nk.integrated.annotated.rds")
  saveRDS(obj, file = output_rds_path)
  message(sprintf("保存注释后的对象到: %s", output_rds_path))
  
  # 生成注释报告
  report_path <- file.path(paths_config$output$reports_dir, "cell_annotation_report.md")
  
  # 创建报告内容
  report_lines <- c(
    "# 细胞类型注释报告",
    "",
    sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## 处理参数",
    sprintf("- 注释方法: SingleR"),
    sprintf("- 参考数据集: %s", ref_dataset_name),
    sprintf("- 选用assay: %s", chosen_assay),
    sprintf("- 数据类型: %s", config$singler$assay_type),
    sprintf("- 标签类型: %s", config$singler$labels),
    "",
    "## 输入数据",
    sprintf("- 输入文件: %s", basename(input_path)),
    sprintf("- 细胞数: %d", ncol(obj)),
    sprintf("- 基因数: %d", nrow(obj)),
    "",
    "## 注释结果",
    sprintf("- 使用的标签列: %s", label_col),
    sprintf("- 识别的细胞类型数: %d", nrow(cell_type_stats)),
    sprintf("- 最主要的细胞类型: %s (%d个细胞, %.2f%%)", 
            cell_type_stats[[label_col]][1], cell_type_stats$n[1], cell_type_stats$percent[1]),
    "",
    "## 细胞类型统计（前10个）",
    "```",
    paste(utils::capture.output(print(utils::head(cell_type_stats, 10))), collapse = "\n"),
    "```",
    "",
    "## 数据完整性验证",
    sprintf("- 细胞验证: %s", validation_results$cells$status),
    sprintf("- 基因验证: %s", validation_results$genes$status),
    sprintf("- Assay验证: %s", validation_results$assays$status),
    sprintf("- 元数据验证: %s", validation_results$metadata$status),
    "",
    "## 输出文件",
    sprintf("- 注释后对象: %s", basename(output_rds_path)),
    sprintf("- 细胞类型统计: %s", basename(stats_path)),
    sprintf("- 细胞类型图件: %s", paths_config$output$plots_dir),
    "",
    "## 下一步",
    "运行 04_cell_filtering.R 进行特定细胞类型过滤"
  )
  
  writeLines(report_lines, report_path)
  message(sprintf("注释报告保存到: %s", report_path))
  
  # 打印摘要
  message("=== 细胞类型注释完成 ===")
  message(sprintf("注释细胞数: %d", ncol(obj)))
  message(sprintf("识别细胞类型数: %d", nrow(cell_type_stats)))
  message(sprintf("最主要细胞类型: %s (%.2f%%)", 
                 cell_type_stats[[label_col]][1], cell_type_stats$percent[1]))
  
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