#!/usr/bin/env Rscript

#' 绘图工具函数
#' 用于scRNA-seq数据处理流程中的可视化

# 加载必要的包
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(viridis)
  library(scales)
})

#' 生成UMAP按时间点分面图
#' @param obj Seurat对象
#' @param group_by 分组变量
#' @param split_by 分面变量
#' @param label_size 标签大小
#' @param ncol 分面列数
#' @return ggplot对象
plot_umap_by_timepoint <- function(obj, group_by = "seurat_clusters", 
                               split_by = "timepoint", 
                               label_size = 4, ncol = 2) {
  p <- DimPlot(
    obj,
    reduction = "umap",
    group.by = group_by,
    split.by = split_by,
    label = TRUE,
    repel = TRUE,
    label.size = label_size,
    ncol = ncol
  ) + 
    theme(strip.text.x = element_text(size = 12)) +
    ggtitle(sprintf("UMAP by %s (split=%s)", group_by, split_by))
  
  return(p)
}

#' 生成簇组成趋势图
#' @param composition_df 簇组成数据框
#' @param timepoint_col 时间点列名
#' @param cluster_col 簇列名
#' @param count_col 计数列名
#' @return ggplot对象
plot_cluster_composition_trends <- function(composition_df, 
                                       timepoint_col = "timepoint",
                                       cluster_col = "cluster", 
                                       count_col = "n") {
  # 计算百分比
  composition_df <- composition_df %>%
    group_by(!!sym(timepoint_col)) %>%
    mutate(percent = 100 * !!sym(count_col) / sum(!!sym(count_col))) %>%
    ungroup()
  
  p <- ggplot(composition_df, 
               aes(x = !!sym(timepoint_col), 
                   y = percent, 
                   group = !!sym(cluster_col), 
                   color = !!sym(cluster_col))) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    labs(title = "Cluster Composition Trends Across Timepoints",
         x = "Timepoint", 
         y = "Percentage of cells",
         color = "Cluster") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    scale_color_brewer(palette = "Set2")
  
  return(p)
}

#' 生成堆叠柱状图
#' @param composition_df 簇组成数据框
#' @param timepoint_col 时间点列名
#' @param cluster_col 簇列名
#' @param count_col 计数列名
#' @param position 位置类型("fill"或"stack")
#' @return ggplot对象
plot_stacked_bar <- function(composition_df, 
                          timepoint_col = "timepoint",
                          cluster_col = "cluster", 
                          count_col = "n",
                          position = "fill") {
  # 计算百分比（如果是fill位置）
  if (position == "fill") {
    composition_df <- composition_df %>%
      group_by(!!sym(timepoint_col)) %>%
      mutate(percent = 100 * !!sym(count_col) / sum(!!sym(count_col))) %>%
      ungroup()
    y_var <- "percent"
    y_label <- "Percentage of cells"
  } else {
    y_var <- count_col
    y_label <- "Cell Count"
  }
  
  p <- ggplot(composition_df, 
               aes(x = !!sym(timepoint_col), 
                   y = !!sym(y_var), 
                   fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", position = position) +
    labs(title = sprintf("Cluster Composition (%s)", 
                     if (position == "fill") "Stacked Percentage" else "Stacked Count"),
         x = "Timepoint", 
         y = y_label,
         fill = "Cluster") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}

#' 生成参数调优热力图
#' @param metrics_df 指标数据框
#' @param x_var x轴变量名
#' @param y_var y轴变量
#' @param fill_var 填充变量名
#' @param title 图标题
#' @return ggplot对象
plot_parameter_heatmap <- function(metrics_df, 
                            x_var = "res", 
                            y_var = "dims", 
                            fill_var = "median_silhouette",
                            title = "Parameter Tuning Heatmap") {
  # 转换数据类型
  metrics_df[[x_var]] <- as.numeric(metrics_df[[x_var]])
  metrics_df[[y_var]] <- as.factor(metrics_df[[y_var]])
  
  p <- ggplot(metrics_df, aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(fill_var))) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "C", na.value = "#f0f0f0") +
    labs(title = title,
         x = str_to_title(x_var),
         y = str_to_title(y_var),
         fill = str_to_title(fill_var)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

#' 生成细胞类型比例图
#' @param cell_type_df 细胞类型数据框
#' @param label_col 标签列名
#' @param count_col 计数列名
#' @param top_n 显示前N个细胞类型
#' @return ggplot对象
plot_cell_type_proportions <- function(cell_type_df, 
                                   label_col = "label", 
                                   count_col = "n",
                                   top_n = 10) {
  # 按计数排序并取前N个
  cell_type_df <- cell_type_df %>%
    arrange(desc(!!sym(count_col))) %>%
    slice_head(n = top_n)
  
  # 获取排序后的标签顺序
  label_levels <- cell_type_df[[label_col]]
  
  p <- ggplot(cell_type_df, aes(x = factor(!!sym(label_col), levels = label_levels), 
                               y = !!sym(count_col))) +
    geom_col(fill = "#3182bd") +
    coord_flip() +
    labs(title = sprintf("Top %d Cell Types by Count", top_n),
         x = "Cell Type",
         y = "Cell Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  return(p)
}

#' 保存多种格式的图件
#' @param plot ggplot对象
#' @param filename 文件名（不含扩展名）
#' @param output_dir 输出目录
#' @param width 图宽度
#' @param height 图高度
#' @param dpi 分辨率
#' @param formats 输出格式列表
save_publication_plots <- function(plot, filename, output_dir, 
                              width = 10, height = 8, dpi = 300,
                              formats = c("png", "pdf", "svg")) {
  for (format in formats) {
    output_path <- file.path(output_dir, sprintf("%s.%s", filename, format))
    
    if (format == "png") {
      ggsave(output_path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
    } else if (format == "pdf") {
      ggsave(output_path, plot = plot, width = width, height = height, device = "pdf")
    } else if (format == "svg") {
      ggsave(output_path, plot = plot, width = width, height = height, device = "svg")
    }
    
    message(sprintf("Saved: %s", output_path))
  }
}

#' 生成组合图
#' @param plot_list 图对象列表
#' @param layout 布局类型
#' @param heights 高度比例
#' @return patchwork对象
create_combined_plot <- function(plot_list, layout = "vertical", heights = NULL) {
  if (length(plot_list) == 2) {
    if (layout == "vertical") {
      combined <- plot_list[[1]] / plot_list[[2]]
    } else if (layout == "horizontal") {
      combined <- plot_list[[1]] | plot_list[[2]]
    }
  } else if (length(plot_list) > 2) {
    combined <- wrap_plots(plot_list, ncol = 1)
  } else {
    stop("Need at least 2 plots to combine")
  }
  
  if (!is.null(heights)) {
    combined <- combined + plot_layout(heights = heights)
  }
  
  return(combined)
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

#' 生成质量控制图
#' @param obj Seurat对象
#' @param qc_metrics QC指标列表
#' @return 图对象列表
generate_qc_plots <- function(obj, qc_metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
  plot_list <- list()
  
  for (metric in qc_metrics) {
    if (metric %in% colnames(obj@meta.data)) {
      p <- VlnPlot(obj, features = metric, pt.size = 0.1) +
        ggtitle(str_to_title(metric)) +
        theme_minimal()
      
      plot_list[[metric]] <- p
    }
  }
  
  return(plot_list)
}

#' 生成双胞检测结果图
#' @param obj Seurat对象
#' @param score_col 双胞评分列名
#' @param class_col 双胞分类列名
#' @return 图对象列表
generate_doublet_plots <- function(obj, score_col = "scDblFinder.score", 
                               class_col = "scDblFinder.class") {
  plot_list <- list()
  
  # 双胞评分密度图
  if (score_col %in% colnames(obj@meta.data)) {
    p_score <- ggplot(obj@meta.data, aes_string(x = score_col)) +
      geom_density(fill = "#2c7fb8", alpha = 0.4) +
      labs(title = "Doublet Score Distribution",
           x = "Doublet Score", y = "Density") +
      theme_minimal()
    
    plot_list[["score_density"]] <- p_score
  }
  
  # 双胞分类柱状图
  if (class_col %in% colnames(obj@meta.data)) {
    class_counts <- table(obj@meta.data[[class_col]])
    class_df <- data.frame(
      class = names(class_counts),
      count = as.numeric(class_counts)
    )
    
    p_class <- ggplot(class_df, aes(x = class, y = count, fill = class)) +
      geom_col() +
      labs(title = "Doublet Classification",
           x = "Cell Class", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
    
    plot_list[["class_counts"]] <- p_class
  }
  
  return(plot_list)
}

message("Plotting utility functions loaded successfully!")