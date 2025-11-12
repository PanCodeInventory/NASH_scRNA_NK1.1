#!/usr/bin/env Rscript
# 基于已保存的 RDS（nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds）
# 生成按时间点分割（split.by = "timepoint"）且按 seurat_clusters 着色的 UMAP 分面图
# 输出：PNG 与 SVG 保存到 2_DataProcessing/3_UMAP-Tuning/plots

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

# 路径设置
rds_path <- "/home/harry/NASH/scRNA-seq/2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds"
plots_dir <- "/home/harry/NASH/scRNA-seq/2_DataProcessing/3_UMAP-Tuning/plots"

# 确保输出目录存在
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

message("[1/3] 读取 RDS 对象：", rds_path)
if (!file.exists(rds_path)) {
  stop("未找到 RDS 文件：", rds_path)
}
obj <- readRDS(rds_path)

# 检查必须的元数据列
message("[2/3] 检查元数据列与降维结果...")
if (!("timepoint" %in% colnames(obj@meta.data))) {
  stop("对象的 meta.data 中不包含 'timepoint' 列，无法按时间点分割绘图。")
}
if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
  message("未发现 'seurat_clusters'，尝试基于 PCA 计算邻居与聚类（resolution=0.3）...")
  if (!("pca" %in% Reductions(obj))) {
    obj <- RunPCA(obj, verbose = FALSE)
  }
  obj <- FindNeighbors(obj, dims = 1:10)
  obj <- FindClusters(obj, resolution = 0.3)
}

# 确保存在 UMAP 降维
if (!("umap" %in% Reductions(obj))) {
  message("未发现 'umap' 降维，尝试运行 RunUMAP(dims=1:10)...")
  if (!("pca" %in% Reductions(obj))) {
    obj <- RunPCA(obj, verbose = FALSE)
  }
  obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
}

# 设定时间点的因子顺序（若存在）
tp_levels <- c("0W_NCD","1W_MCD","2W_MCD","6W_MCD")
if (is.character(obj$timepoint) || is.factor(obj$timepoint)) {
  obj$timepoint <- factor(as.character(obj$timepoint), levels = tp_levels)
}

# 生成分面 UMAP 图
message("[3/3] 生成按时间点分割的 UMAP 分面图并保存 PNG/SVG...")
p_split <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "timepoint",
  label = TRUE,
  repel = TRUE,
  ncol = 2
) + theme(strip.text.x = element_text(size = 12))

# 保存出版级别图件
output_filename <- file.path(plots_dir, "UMAP_NoCluster6_Clusters_by_Timepoint")
plot_width <- 12
plot_height <- 10

ggsave(
  filename = paste0(output_filename, ".png"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  dpi = 500
)

ggsave(
  filename = paste0(output_filename, ".svg"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  device = "svg"
)

# 额外保存一个与命名约定一致的最终文件名
ggsave(
  filename = file.path(plots_dir, "UMAP_noCluster6_final_byTimepoint.png"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  dpi = 500
)

# 新增：按时间点分别导出独立 UMAP 图（PNG/SVG）
for (tp in tp_levels) {
  obj_tp <- subset(obj, subset = timepoint == tp)
  p_tp <- DimPlot(
    obj_tp,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) + ggtitle(paste("NK Cells UMAP -", tp))
  tp_safe <- gsub("[^A-Za-z0-9_-]", "_", tp)
  out_tp <- file.path(plots_dir, paste0("UMAP_NoCluster6_Clusters_", tp_safe))
  ggsave(paste0(out_tp, ".png"), p_tp, width = 10, height = 8, dpi = 500)
  ggsave(paste0(out_tp, ".svg"), p_tp, width = 10, height = 8, device = "svg")
}

# 新增：不同簇在不同时间点之间的变化趋势图（折线）
cluster_df <- obj@meta.data %>%
  dplyr::select(seurat_clusters, timepoint) %>%
  dplyr::mutate(
    seurat_clusters = as.factor(seurat_clusters),
    timepoint = factor(timepoint, levels = tp_levels)
  ) %>%
  dplyr::count(timepoint, seurat_clusters, name = "n") %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::ungroup()

p_trend <- ggplot(cluster_df, aes(x = timepoint, y = percent, group = seurat_clusters, color = seurat_clusters)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "NK Cluster Composition Trends Across Timepoints (NoCluster6)",
       x = "Timepoint", y = "Percentage of cells") +
  theme_bw() +
  theme(legend.position = "right")

trend_out <- file.path(plots_dir, "UMAP_NoCluster6_Cluster_Trends")
ggsave(paste0(trend_out, ".png"), p_trend, width = 12, height = 8, dpi = 500)
ggsave(paste0(trend_out, ".svg"), p_trend, width = 12, height = 8, device = "svg")

# 新增：堆叠柱状图（按 timepoint 汇总各簇百分比）
p_bar_stacked <- ggplot(cluster_df, aes(x = timepoint, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "NK Cluster Composition (Stacked Bar, NoCluster6)",
       x = "Timepoint", y = "Percentage of cells") +
  theme_bw() +
  theme(legend.position = "right")

bar_out <- file.path(plots_dir, "UMAP_NoCluster6_Cluster_Composition_StackedBar")
ggsave(paste0(bar_out, ".png"), p_bar_stacked, width = 12, height = 8, dpi = 500)
ggsave(paste0(bar_out, ".svg"), p_bar_stacked, width = 12, height = 8, device = "svg")

# 导出簇构成统计表（每时间点各簇的细胞数与百分比）
write.csv(cluster_df, file.path(plots_dir, "NoCluster6_cluster_composition_by_timepoint.csv"), row.names = FALSE)

message("已保存分面 UMAP 图至：", plots_dir)
