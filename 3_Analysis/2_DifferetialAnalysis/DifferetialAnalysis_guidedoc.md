## 时间点驱动的差异分析（整体 + 按簇）

### 概览
- 不依赖 orig.ident 或 condition；以 `meta.data$timepoint` 作为唯一处理条件，比较组由 `comparisons` 指定（示例：`0W_NCD_vs_1W_MCD`、`0W_NCD_vs_2W_MCD`、`0W_NCD_vs_6W_MCD`）。
- 整体分析（overall）：在选定比较下合并所有细胞，按 `timepoint` 做伪样本聚合并进行 DESeq2 差异分析，输出前缀 `all_cells.*`。
- 按簇分析（by_cluster）：按 `meta.data$seurat_clusters` 分层，每个簇内再按 `timepoint` 做伪样本聚合并进行 DESeq2 差异分析，输出前缀 `<cluster>.*`。
- 伪重复（pseudoreplicates）：由于每组只有一个时间点，启用 `replicates_per_group=3` 将每个时间点的细胞拆分为 K=3 个子集（伪样本），以满足 DESeq2 对离散度估计的需求（可在配置中调整或关闭）。
- 物种与富集：`viz.enrichment.species` 使用 `Mus musculus`（小鼠），数据库映射 HALLMARK→`H`，FGSEA 以 `ranked_genes.csv` 的 `stat` 作为评分。

### 配置关键项（3_Analysis/2_DifferetialAnalysis/scripts/config.yaml）
```
columns:
  sample_id: "timepoint"
  cluster: "seurat_clusters"
  timepoint: "timepoint"

group_label_format: "{timepoint}"

deg:
  method: "DESeq2"
  scope: "both"                    # overall + by_cluster
  lfc_threshold: 0.25
  p_adj_cutoff: 0.05
  pseudobulk:
    aggregation: "sum"
    min_cells_per_cluster: 50
    replicates_per_group: 3        # 伪重复K（可调）
  filter_expr:
    type: "counts"
    min_count: 10
    min_samples: 2
  shrinkage:
    enabled: true
    type: "apeglm"

comparisons:
  - id: "0W_NCD_vs_1W_MCD"
    control: "0W_NCD"
    case: "1W_MCD"
  - id: "0W_NCD_vs_2W_MCD"
    control: "0W_NCD"
    case: "2W_MCD"
  - id: "0W_NCD_vs_6W_MCD"
    control: "0W_NCD"
    case: "6W_MCD"

viz:
  enrichment:
    enabled: true
    databases: ["HALLMARK"]        # 映射为类别 "H"
    species: "Mus musculus"        # 小鼠
    min_size: 10
    max_size: 500
    rank_metric: "stat"
```


### 输出结构示例
- 差异表达结果表：`results/tables/<cmp_id>/`
  - `all_cells.DESeq2.results.csv`、`<cluster>.DESeq2.results.csv`
- 标准化计数与排名基因：`results/data/<cmp_id>/`
  - `all_cells.normalized_counts.csv`、`all_cells.ranked_genes.csv` 等
- 图形：`results/plots/<type>/<cmp_id>/`（`type = volcano|ma|heatmap|enrichment`）
  - 例：`results/plots/volcano/0W_NCD_vs_1W_MCD/0.volcano.png`
- 富集：`results/enrichment/<cmp_id>/` 与 `results/plots/enrichment/<cmp_id>/`
  - 例：`0.fgsea.results.csv`、`0.fgsea.hallmark.png`

### 调参建议
- 显著性阈值：`deg.lfc_threshold`（默认 0.25）、`deg.p_adj_cutoff`（默认 0.05）。
- 表达过滤：`deg.filter_expr.min_count/min_samples`（控制低表达基因过滤强度）。
- 伪重复：`deg.pseudobulk.replicates_per_group`（默认 3；若某簇细胞较少，可降至 2）。
- 簇最小细胞数：`deg.pseudobulk.min_cells_per_cluster=50`（不足则跳过并记录到 `data/cell_counts_threshold_failures.tsv`）。
- 物种与数据库：`viz.enrichment.species` 与 `viz.enrichment.databases`（已映射 HALLMARK→`H`）。

## 故障排查（FAQ）

- 报错 “missing value where TRUE/FALSE needed”
  - 原因：因子与空字符串比较、或顺序拼接含 NA/空值导致 if/while 接收到 NA；
  - 处理：脚本已改为“先字符过滤、再因子化”，并清洗与回退时间点顺序；可通过 `--timepoint-order` 显式指定顺序。
- Seurat v5 警告 “slot 已废弃、请用 layer”
  - 属正常版本提示；分析脚本兼容 layer/slot 接口，已在关键节点做回退与检查。
- FindAllMarkers 报错或返回 0 行
  - 请确认所选 assay 的 data/VariableFeatures 非空；可通过参数化脚本的 integrated→SCT→RNA 回退策略或在 RNA 上自动准备；
  - 可调整 `--min-pct`、`--logfc-threshold` 或 `--test-use "MAST"`。
- 运行缓慢
  - 建议安装 presto 包，并根据机器资源配置并行；当前参数化脚本已默认启用稳定执行策略。

## 主要脚本（当前有效）
- `2_DataProcessing/Scripts/remove_doublets_and_contaminants.R`：去双胞 + 注释 + UCell + 去污染 + 重分析主流程
- `2_DataProcessing/Scripts/singleR_annotation_fix.R`：SingleR 空 data 层修复（counts→logNormCounts→SingleR）
- `2_DataProcessing/Scripts/remove_NKT_cells.R`：在已注释对象上剔除 NKT 并生成报告与图件
- `2_DataProcessing/Scripts/tune_noNKT_dims_resolution.R`：基于 noNKT 对象进行 dims × resolution 调参、选择并重跑生成最终产物
- `3_Analysis/Scripts/export_cluster_proportions.R`：簇比例导出与绘图
- `3_Analysis/Scripts/find_cluster_markers.R`：参数化差异基因（integrated→SCT→RNA 回退）
- `3_Analysis/Scripts/find_markers_simple.R`：简化版差异基因（快速产出）
- 历史脚本（仍可参考）：`Files/UMAP/scripts/*`

## 主要结果（样例）
- 去双胞/清理整体：
  - `2_DataProcessing/2_Doublet_Removed/reports/cleaning_report.md`
  - `2_DataProcessing/2_Doublet_Removed/plots/UMAP_filtered_by_SingleR.png`
  - `2_DataProcessing/2_Doublet_Removed/plots/UMAP_filtered_clusters_by_Timepoint.png`
- NKT 剔除：
  - `2_DataProcessing/2_Doublet_Removed/reports/noNKT_removal_report.md`
  - `2_DataProcessing/2_Doublet_Removed/plots/NKT_removal_label_counts_before_after.png`
- 调参与最终：
  - `2_DataProcessing/3_UMAP-Tuning/data/nk_noNKT_tuning_metrics.csv`
  - `2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_final_byTimepoint.png`
  - `2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_final_bySingleR.png`
- 下游分析（3_Analysis）：
  - `3_Analysis/1.ClusterAnalysis/data/cluster_counts_by_timepoint.csv`
  - `3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv`
  - `3_Analysis/1.ClusterAnalysis/plots/cluster_proportion_lineplot.(png|pdf)`
  - `3_Analysis/1.ClusterAnalysis/plots/cluster_composition_stackedbar.(png|pdf)`
  - `3_Analysis/1.ClusterAnalysis/data/markers_all_clusters.csv`
  - `3_Analysis/1.ClusterAnalysis/data/markers_top10_per_cluster.csv`