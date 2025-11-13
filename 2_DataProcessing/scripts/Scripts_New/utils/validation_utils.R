#!/usr/bin/env Rscript

#' 验证工具函数
#' 用于scRNA-seq数据处理流程中的质量验证

# 加载必要的包
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

#' 验证Seurat对象完整性
#' @param obj Seurat对象
#' @param check_items 检查项目列表
#' @return 验证结果列表
validate_seurat_object <- function(obj, check_items = c("cells", "genes", "assays", "metadata")) {
  validation_results <- list()
  
  if ("cells" %in% check_items) {
    validation_results$cells <- list(
      count = ncol(obj),
      status = if (ncol(obj) > 0) "PASS" else "FAIL"
    )
  }
  
  if ("genes" %in% check_items) {
    validation_results$genes <- list(
      count = nrow(obj),
      status = if (nrow(obj) > 0) "PASS" else "FAIL"
    )
  }
  
  if ("assays" %in% check_items) {
    assays <- Seurat::Assays(obj)
    validation_results$assays <- list(
      available = assays,
      count = length(assays),
      status = if (length(assays) > 0) "PASS" else "FAIL"
    )
  }
  
  if ("metadata" %in% check_items) {
    metadata_cols <- colnames(obj@meta.data)
    validation_results$metadata <- list(
      columns = metadata_cols,
      count = length(metadata_cols),
      status = if (length(metadata_cols) > 0) "PASS" else "FAIL"
    )
  }
  
  return(validation_results)
}

#' 检查数据完整性
#' @param obj Seurat对象
#' @param required_columns 必需的元数据列
#' @param required_assays 必需的assay
#' @return 检查结果列表
check_data_completeness <- function(obj, required_columns = c("timepoint", "seurat_clusters"), 
                               required_assays = c("RNA")) {
  completeness_results <- list()
  
  # 检查元数据列
  missing_columns <- setdiff(required_columns, colnames(obj@meta.data))
  completeness_results$metadata <- list(
    required = required_columns,
    available = colnames(obj@meta.data),
    missing = missing_columns,
    complete = length(missing_columns) == 0,
    status = if (length(missing_columns) == 0) "PASS" else "FAIL"
  )
  
  # 检查assay
  available_assays <- Seurat::Assays(obj)
  missing_assays <- setdiff(required_assays, available_assays)
  completeness_results$assays <- list(
    required = required_assays,
    available = available_assays,
    missing = missing_assays,
    complete = length(missing_assays) == 0,
    status = if (length(missing_assays) == 0) "PASS" else "FAIL"
  )
  
  return(completeness_results)
}

#' 比较处理前后的数据
#' @param before_obj 处理前的Seurat对象
#' @param after_obj 处理后的Seurat对象
#' @param step_name 处理步骤名称
#' @return 比较结果列表
compare_before_after <- function(before_obj, after_obj, step_name = "processing") {
  comparison_results <- list()
  
  # 细胞数比较
  before_cells <- ncol(before_obj)
  after_cells <- ncol(after_obj)
  cells_removed <- before_cells - after_cells
  cells_removed_pct <- if (before_cells > 0) 100 * cells_removed / before_cells else 0
  
  comparison_results$cell_counts <- list(
    before = before_cells,
    after = after_cells,
    removed = cells_removed,
    removed_percent = cells_removed_pct,
    step = step_name
  )
  
  # 基因数比较
  before_genes <- nrow(before_obj)
  after_genes <- nrow(after_obj)
  genes_lost <- before_genes - after_genes
  
  comparison_results$gene_counts <- list(
    before = before_genes,
    after = after_genes,
    lost = genes_lost,
    step = step_name
  )
  
  # 元数据列比较
  before_metadata <- colnames(before_obj@meta.data)
  after_metadata <- colnames(after_obj@meta.data)
  metadata_added <- setdiff(after_metadata, before_metadata)
  metadata_lost <- setdiff(before_metadata, after_metadata)
  
  comparison_results$metadata_changes <- list(
    before = before_metadata,
    after = after_metadata,
    added = metadata_added,
    lost = metadata_lost,
    step = step_name
  )
  
  return(comparison_results)
}

#' 生成QC报告
#' @param validation_results 验证结果列表
#' @param output_path 输出路径
#' @param step_name 处理步骤名称
generate_qc_report <- function(validation_results, output_path, step_name = "QC") {
  report_lines <- c(
    sprintf("# %s Report", step_name),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Validation Results",
    ""
  )
  
  # 添加各项验证结果
  for (item_name in names(validation_results)) {
    item_result <- validation_results[[item_name]]
    report_lines <- c(report_lines,
      sprintf("### %s", str_to_title(item_name)),
      sprintf("- Status: %s", item_result$status),
      sprintf("- Details: %s", paste(names(item_result), collapse = ", ")),
      ""
    )
  }
  
  # 写入报告文件
  writeLines(report_lines, output_path)
  message(sprintf("QC report saved to: %s", output_path))
}

#' 验证下游兼容性
#' @param obj Seurat对象
#' @param expected_clusters 期望的簇列表
#' @param marker_genes 标记基因列表
#' @return 兼容性验证结果
validate_downstream_compatibility <- function(obj, expected_clusters = c(0, 1, 2, 3, 4, 5, 6),
                                          marker_genes = list(
                                            cluster_5 = c("Mki67", "H2afx", "Cdca3", "Hist1h3c"),
                                            cluster_6 = c("Cd79a", "Iglc3", "Ebf1", "Ms4a1")
                                          )) {
  compatibility_results <- list()
  
  # 检查簇完整性
  actual_clusters <- sort(unique(as.numeric(as.character(obj$seurat_clusters))))
  missing_clusters <- setdiff(expected_clusters, actual_clusters)
  extra_clusters <- setdiff(actual_clusters, expected_clusters)
  
  compatibility_results$clusters <- list(
    expected = expected_clusters,
    actual = actual_clusters,
    missing = missing_clusters,
    extra = extra_clusters,
    complete = length(missing_clusters) == 0 && length(extra_clusters) == 0,
    status = if (length(missing_clusters) == 0 && length(extra_clusters) == 0) "PASS" else "FAIL"
  )
  
  # 检查标记基因表达
  if ("RNA" %in% Seurat::Assays(obj)) {
    marker_validation <- list()
    
    for (cluster_name in names(marker_genes)) {
      cluster_num <- as.numeric(gsub("cluster_", "", cluster_name))
      if (cluster_num %in% actual_clusters) {
        genes <- marker_genes[[cluster_name]]
        expression_values <- c()
        
        for (gene in genes) {
          if (gene %in% rownames(obj)) {
            cluster_cells <- WhichCells(obj, expression = seurat_clusters == cluster_num)
            if (length(cluster_cells) > 0) {
              expr <- AverageExpression(obj, features = gene, cells = cluster_cells)
              expression_values <- c(expression_values, expr[[gene]])
            }
          }
        }
        
        marker_validation[[cluster_name]] <- list(
          cluster = cluster_num,
          marker_genes = genes,
          expression_values = expression_values,
          detected = length(expression_values) > 0 && mean(expression_values) > 0
        )
      }
    }
    
    compatibility_results$markers <- marker_validation
  }
  
  # 检查时间点完整性
  if ("timepoint" %in% colnames(obj@meta.data)) {
    timepoints <- unique(obj$timepoint)
    compatibility_results$timepoints <- list(
      available = timepoints,
      count = length(timepoints),
      expected = c("0W_NCD", "1W_MCD", "2W_MCD", "6W_MCD"),
      complete = all(c("0W_NCD", "1W_MCD", "2W_MCD", "6W_MCD") %in% timepoints)
    )
  }
  
  return(compatibility_results)
}

#' 生成兼容性验证报告
#' @param compatibility_results 兼容性验证结果
#' @param output_path 输出路径
generate_compatibility_report <- function(compatibility_results, output_path) {
  expected_clusters <- compatibility_results$clusters$expected
  actual_clusters <- compatibility_results$clusters$actual
  missing_clusters <- compatibility_results$clusters$missing
  extra_clusters <- compatibility_results$clusters$extra
  cluster_status <- compatibility_results$clusters$status
  
  report_lines <- c(
    "# Downstream Compatibility Validation Report",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Cluster Validation",
    sprintf("- Expected clusters: %s", paste(expected_clusters, collapse = ", ")),
    sprintf("- Actual clusters: %s", paste(actual_clusters, collapse = ", ")),
    sprintf("- Missing clusters: %s", paste(missing_clusters, collapse = ", ")),
    sprintf("- Extra clusters: %s", paste(extra_clusters, collapse = ", ")),
    sprintf("- Status: %s", cluster_status),
    ""
  )
  
  # 添加标记基因验证结果
  if ("markers" %in% names(compatibility_results)) {
    report_lines <- c(report_lines,
      "## Marker Gene Validation",
      ""
    )
    
    for (cluster_name in names(compatibility_results$markers)) {
      marker_info <- compatibility_results$markers[[cluster_name]]
      detected_status <- if (marker_info$detected) "YES" else "NO"
      expr_values <- if (length(marker_info$expression_values) > 0) {
        paste(round(marker_info$expression_values, 3), collapse = ", ")
      } else {
        "Not detected"
      }
      
      report_lines <- c(report_lines,
        sprintf("### Cluster %s", marker_info$cluster),
        sprintf("- Marker genes: %s", paste(marker_info$marker_genes, collapse = ", ")),
        sprintf("- Detected: %s", detected_status),
        sprintf("- Average expression: %s", expr_values),
        ""
      )
    }
  }
  
  # 添加时间点验证结果
  if ("timepoints" %in% names(compatibility_results)) {
    timepoint_status <- if (compatibility_results$timepoints$complete) "YES" else "NO"
    report_lines <- c(report_lines,
      "## Timepoint Validation",
      sprintf("- Available timepoints: %s", paste(compatibility_results$timepoints$available, collapse = ", ")),
      sprintf("- Expected timepoints: %s", paste(compatibility_results$timepoints$expected, collapse = ", ")),
      sprintf("- Complete: %s", timepoint_status),
      ""
    )
  }
  
  # 写入报告文件
  writeLines(report_lines, output_path)
  message(sprintf("Compatibility report saved to: %s", output_path))
}

#' 字符串转换为标题格式
#' @param str 字符串
#' @return 格式化后的字符串
str_to_title <- function(str) {
  # 将下划线转换为空格，首字母大写
  str <- gsub("_", " ", str)
  str <- paste(toupper(substring(str, 1, 1)), 
              tolower(substring(str, 2)), sep = "")
  return(str)
}

#' 验证文件路径存在性
#' @param file_paths 文件路径列表
#' @param file_descriptions 文件描述列表
#' @return 验证结果列表
validate_file_paths <- function(file_paths, file_descriptions = NULL) {
  validation_results <- list()
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[[i]]
    file_exists <- file.exists(file_path)
    
    result <- list(
      path = file_path,
      exists = file_exists,
      status = if (file_exists) "PASS" else "FAIL"
    )
    
    if (!is.null(file_descriptions) && i <= length(file_descriptions)) {
      result$description <- file_descriptions[[i]]
    }
    
    validation_results[[i]] <- result
  }
  
  return(validation_results)
}

#' 生成验证摘要
#' @param validation_results 验证结果列表
#' @return 摘要数据框
generate_validation_summary <- function(validation_results) {
  summary_df <- data.frame(
    item = names(validation_results),
    status = sapply(validation_results, function(x) x$status),
    stringsAsFactors = FALSE
  )
  
  overall_status <- if (all(summary_df$status == "PASS")) "PASS" else "FAIL"
  
  summary_list <- list(
    summary = summary_df,
    overall_status = overall_status,
    total_items = nrow(summary_df),
    passed_items = sum(summary_df$status == "PASS"),
    failed_items = sum(summary_df$status == "FAIL")
  )
  
  return(summary_list)
}

message("Validation utility functions loaded successfully!")