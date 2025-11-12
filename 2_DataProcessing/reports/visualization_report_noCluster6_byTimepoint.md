# 可视化报告：NoCluster6 时间点分面 UMAP

- 数据对象：`2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds`
- 生成脚本：`2_DataProcessing/Scripts/plot_umap_from_rds_noCluster6_timepoint_split.R`
- 方法概述：
  - 直接读取已保存的 Seurat 对象，不进行重新整合。
  - 若缺失聚类或降维结果（`seurat_clusters`/`umap`），脚本会自动运行 `RunPCA`、`FindNeighbors`、`FindClusters(resolution=0.3)`、`RunUMAP(dims=1:10)` 以补齐。
  - 时间点因子顺序设为：`0W_NCD` → `1W_MCD` → `2W_MCD` → `6W_MCD`，保证分面面板顺序一致。
  - 可视化调用：`DimPlot(obj, reduction="umap", group.by="seurat_clusters", split.by="timepoint", label=TRUE, repel=TRUE, ncol=2)`。

## 输出图件

下列图件已由脚本生成，可用于展示 NK（去除 cluster6 后）的时间点分面 UMAP：

1) 分面 UMAP（PNG，高分辨率）
   - 路径：`2_DataProcessing/3_UMAP-Tuning/plots/UMAP_NoCluster6_Clusters_by_Timepoint.png`
   - 预览：
   ![](../3_UMAP-Tuning/plots/UMAP_NoCluster6_Clusters_by_Timepoint.png)

2) 分面 UMAP（SVG，适合矢量编辑）
   - 路径：`2_DataProcessing/3_UMAP-Tuning/plots/UMAP_NoCluster6_Clusters_by_Timepoint.svg`
   - 预览（SVG 在部分查看器中需专门打开）：
   ![](../3_UMAP-Tuning/plots/UMAP_NoCluster6_Clusters_by_Timepoint.svg)

3) 最终命名约定版本（PNG）
   - 路径：`2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noCluster6_final_byTimepoint.png`
   - 预览：
   ![](../3_UMAP-Tuning/plots/UMAP_noCluster6_final_byTimepoint.png)

## 解读要点（建议按下列维度检查）

- 面板顺序与标签
  - 面板顺序是否符合 0W_NCD → 1W_MCD → 2W_MCD → 6W_MCD？
  - 分面标签字体与可读性（`strip.text.x` 已设定为 size = 12）。
- 聚类稳定性与空间形态
  - `seurat_clusters` 在不同时间点的空间位置与形态是否保持一致？
  - 是否存在某些簇在 MCD 条件下（1W/2W/6W）明显扩张或收缩的趋势？
- 细胞数量与稠密度
  - 面板间细胞数量是否显著差异（可从点密度大致判断）？
  - 是否出现孤立小团或可能的技术噪声簇（需结合 marker 进一步验证）？
- 标注与可读性
  - 簇标签是否覆盖或遮挡关键区域；如需可考虑调整 `label.size` 或开启 `repel`（已开启）。

## 复现方式

- 在项目根目录运行：
  ```bash
  Rscript 2_DataProcessing/Scripts/plot_umap_from_rds_noCluster6_timepoint_split.R
  ```
- 输出文件将保存至：`2_DataProcessing/3_UMAP-Tuning/plots/`

## 相关文件与路径

- RDS：`2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds`
- 脚本：`2_DataProcessing/Scripts/plot_umap_from_rds_noCluster6_timepoint_split.R`
- 图件目录：`2_DataProcessing/3_UMAP-Tuning/plots/`

## 后续建议（可选）

- 若需与先前分析保持完全一致的输出路径，可将脚本中的 `plots_dir` 修改为：`/home/harry/NASH/scRNA-seq/Files/UMAP/Results/plots` 并重新运行。
- 如需更细粒度的展示：
  - 按每个时间点分别导出独立 UMAP 图（与 `generate_umap_nk.R` 的扩展一致）。
  - 计算并导出各时间点的簇构成百分比、趋势折线与堆叠柱状图，以量化时间序列变化。
- 如需美术化（论文级别）：
  - 调整颜色方案、标签字号、面板间距，并在 SVG 上进行后期排版。

（本报告用于汇总与定位上述三张图件，便于审阅与后续图件挑选/修订。）
