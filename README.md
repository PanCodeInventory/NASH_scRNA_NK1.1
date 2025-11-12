# scRNA-seq NASH NK细胞分析项目

## 项目概述

本项目专注于 NASH（非酒精性脂肪性肝炎）疾病模型中 NK 细胞的单细胞 RNA 测序分析。通过分析不同时间点（0 周、1 周、2 周、6 周）的 NK 细胞样本，研究 NASH 疾病进程中 NK 细胞的变化规律。

## 项目结构

```
scRNA-seq/
├── README.md
├── .gitignore
├── 1_Files/                          # 原始/预处理数据（按分组）
│   ├── NK1.1/
│   └── CD45.2/
├── 2_DataProcessing/                 # 主数据处理管线（脚本与产物）
│   ├── 1_Samples_Merging/
│   │   ├── Scripts & guidedoc*.md
│   │   └── Results/{data,plots,rds}
│   ├── 2_Doublet_Removed/            # 去双胞/注释/清理产出（报告与图件）
│   │   ├── RDS/                      # R 对象（清理后 / 单细胞注释 / noNKT / tuned）
│   │   ├── plots/                    # UMAP 等图件
│   │   └── reports/                  # 报告（cleaning / singleR 修复 / noNKT）
│   ├── 3_UMAP-Tuning/                # UMAP 调参与选择产出（metrics/plots/logs）
│   │   ├── data/                     # 调参指标与候选 CSV
│   │   ├── plots/                    # 调参热力图与候选 UMAP
│   │   └── logs/                     # 运行配置与会话信息
│   └── Scripts/                      # 脚本（生成、清理、调参）
│       ├── generate_umap_nk.R
│       ├── generate_umap_nk_post.R
│       ├── remove_doublets_and_contaminants.R
│       ├── singleR_annotation_fix.R
│       ├── remove_NKT_cells.R                  # 新增：在已注释对象上剔除 NKT
│       ├── tune_noNKT_dims_resolution.R        # 新增：基于 noNKT 对象进行 dims × resolution 调参并重跑
│       ├── remove_clusters_and_recompute.R     # 新增：移除指定簇并重新计算UMAP/聚类
│       └── tune_noCluster6_dims_resolution.R   # 新增：基于无Cluster6对象的调参脚本
├── 2_Filter/                         # 可选镜像产出目录（按你的偏好保留）
│   └── 2_Doublet_Removed/{RDS,plots,reports}
└── 3_Analysis/                       # 下游分析
    ├── 1.ClusterAnalysis/            # 簇比例与差异基因分析产出
    │   ├── data/                     # CSV 表格（比例/markers）
    │   ├── plots/                    # 图件（折线/堆叠）
    │   └── logs/                     # 运行日志与会话信息
    └── Scripts/                      # 下游分析脚本
        ├── export_cluster_proportions.R        # 按时间点×簇统计并绘图
        ├── find_cluster_markers.R             # 每簇差异基因（CLI 参数版，含回退策略）
        └── find_markers_simple.R              # 每簇差异基因（简化版，快速产出）
```

## 样本与分组信息
- 分组：NCD（0W）与 MCD（1W/2W/6W）
- 细胞类型：NK1.1（自然杀伤细胞）

## 分析流程概览

1) 样本合并与整合（SCTransform + Anchors）  
2) 去双胞（scDblFinder，按样本/时间点分组）  
3) 自动注释（SingleR，logcounts 修复策略）  
4) UCell 基因签名评分（NK/T/B/Myeloid/DC/Plasma/Endothelium/Fibroblast/Hepatocyte）  
5) 去污染（细胞级阈值 + 簇级非 NK 占比阈值）  
6) 重跑降维/聚类/UMAP（兼容 SCT/RNA 多模型，必要时回退）  
7) NKT 剔除（基于 SingleR 标签严格规则）  
8) dims × resolution 调参（UMAP/聚类）与最终参数选择  
9) 按最终参数生成分面 UMAP（timepoint）与标签 UMAP（SingleR）

## 技术栈
- R、Seurat、SingleCellExperiment、scDblFinder、SingleR、celldex、scater、UCell、ggplot2、patchwork

## 更新

- 2025-10-20
  - 新增脚本：`Files/UMAP/scripts/remove_doublets_and_contaminants.R`（集成 scDblFinder 去双胞、SingleR 自动注释、UCell 签名评分、去污染规则与重分析的主流程）
  - 修复 Seurat v5 多层 assay 转换为 SCE 的问题，增强元数据行名对齐与日志/报告目录创建的健壮性

- 2025-10-20 深夜
  - 新增特异性 SingleR 修复脚本：`Files/UMAP/scripts/singleR_annotation_fix.R`（从 counts 构建 SCE，scater::logNormCounts 生成 logcounts，显式以 logcounts 作为 SingleR 输入，规避 data 层为空告警）
  - 产出对象与文档：
    - `Files/Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds`
    - `Files/Doublet_Removed/reports/singleR_fix_report.md`
    - `Files/Doublet_Removed/plots/SingleR_label_barplot.png`

- 2025-10-21（清理与重分析）
  - 将 SingleR 修复策略集成至主流程并完成全流程清理与重分析（保留 NK/ILC）
  - 去双胞结果：移除 312 个细胞（约 1.61%）
  - 生成清理后对象与报告（历史路径 Files/*）：
    - `Files/Doublet_Removed/RDS/nk.integrated.filtered.rds`
    - `Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds`
    - `Files/Doublet_Removed/reports/cleaning_report.md`
  - 图件（历史路径 Files/*）：DoubletScore、Cluster_nonNK_fraction、UMAP 等

- 2025-10-21（NKT 剔除 + UMAP 调参与最终）
  - 新增脚本：`2_DataProcessing/Scripts/remove_NKT_cells.R`（严格规则排除 NKT）
  - 新增脚本：`2_DataProcessing/Scripts/tune_noNKT_dims_resolution.R`（调参与最终降维/聚类/UMAP）
  - 产出：
    - NKT 剔除：19126 → 19007（移除 119，0.62%）
    - 调参指标与候选：`2_DataProcessing/3_UMAP-Tuning/data/*`
    - 最终参数：dims=10、res=0.3（见 `selected_params.txt`）
    - 最终对象与图：`2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.noNKT.tuned.rds`；`3_UMAP-Tuning/plots/*`

- 2025-10-21（下游分析 3_Analysis）
  - 新增脚本：`3_Analysis/Scripts/export_cluster_proportions.R`、`3_Analysis/Scripts/find_cluster_markers.R`、`3_Analysis/Scripts/find_markers_simple.R`
  - 产出（示例）：
    - `3_Analysis/1.ClusterAnalysis/data/cluster_counts_by_timepoint.csv`
    - `3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv`
    - `3_Analysis/1.ClusterAnalysis/plots/cluster_proportion_lineplot.(png|pdf)`
    - `3_Analysis/1.ClusterAnalysis/plots/cluster_composition_stackedbar.(png|pdf)`
    - `3_Analysis/1.ClusterAnalysis/data/markers_all_clusters.csv`
    - `3_Analysis/1.ClusterAnalysis/data/markers_top10_per_cluster.csv`
  - 兼容与性能：
    - 时间点顺序与因子/字符比较的稳健处理；integrated→SCT→RNA 的差异分析回退；可选安装 presto 提升速度

- 2025-10-23（Cluster 6污染清理 + 重新调优）
  - 新增脚本：`2_DataProcessing/Scripts/remove_clusters_and_recompute.R`（移除指定簇并重新计算UMAP/聚类）
  - 新增脚本：`2_DataProcessing/Scripts/tune_noCluster6_dims_resolution.R`（基于无Cluster6对象的参数优化）
  - 污染清理：移除596个细胞（19,126→18,530，移除3.1%），主要是B细胞污染
  - 重新调优：最佳参数dims=10, resolution=0.3，轮廓系数0.276
  - 重新聚类：获得7个生物学意义明确的NK细胞亚群
  - 更新分析：
    - `3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv`
    - `3_Analysis/1.ClusterAnalysis/data/markers_top10_per_cluster.csv`
    - `3_Analysis/1.ClusterAnalysis/plots/cluster_proportion_lineplot.(png|pdf)`
  - 核心对象：`2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds`

- 2025-10-23（簇重命名和隐藏策略）
  - 新增脚本：`2_DataProcessing/Scripts/rename_and_hide_clusters.R`（簇重命名和RDS对象更新）
  - 新增脚本：`2_DataProcessing/Scripts/correct_cluster_renaming.R`（分析文件重命名校正）
  - 新增脚本：`2_DataProcessing/Scripts/visualize_with_hidden_clusters.R`（隐藏簇6的可视化）
  - 重命名策略：
    - 原簇5（B细胞污染：Iglc3、Cd79a等）→ 簇6（在可视化中隐藏）
    - 原簇6（增殖：H2afx、Mki67等）→ 簇5（可见）
  - 更新文件：
    - `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds`
    - `3_Analysis/1.ClusterAnalysis/data/*`（所有标记基因和比例文件）
    - `3_Analysis/1.ClusterAnalysis/plots/*_hidden_cluster6.(png|pdf)`（隐藏簇6的图件）
  - 报告文档：
    - `2_DataProcessing/reports/cluster_renaming_final_report.md`（完整执行报告）
    - `2_DataProcessing/reports/visualization_report_hidden_cluster6.md`（可视化报告）
  - 质量控制：所有原始文件完整备份至`2_DataProcessing/reports/backup/`


