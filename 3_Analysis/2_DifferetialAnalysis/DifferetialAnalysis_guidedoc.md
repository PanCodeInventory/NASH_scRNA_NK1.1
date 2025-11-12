# 差异表达分析统一化指南（配置驱动，单脚本遍历全部比较）

本指南将原先“每个比较一个脚本”的差异分析流程，重构为“单一脚本 + 配置驱动”的统一方案；可视化与富集分析保持独立脚本，仅消费分析阶段导出的 CSV 结果。

## 1. 目标与范围
- 目标：以配置文件（config.yaml）定义输入、列名映射、DE 参数与比较列表，由一个主脚本遍历全部比较，完成 pseudobulk + DESeq2 差异分析；可视化阶段读取结果 CSV 绘图与富集。
- 对比组别（默认覆盖，来自 comparisons 列表）：
  - NCD_0w_vs_MCD_1w
  - NCD_0w_vs_MCD_2w
  - NCD_0w_vs_MCD_6w
- 分析语言：R

## 2. 输入数据与前置假设
- 输入 Seurat 对象：2_DataProcessing/RDS 中最终对象（默认：`2_DataProcessing/RDS/nk.integrated.v4.rds`，由 config 指定）。
- 元数据列要求（由 config 指定映射）：
  - sample_id：样本ID（本项目以 timepoint 作为样本名；不同 timepoint 标签即为不同“样本/组别”的聚合）
  - cell_type：细胞类型/簇（示例：singleR.label）
  - timepoint：时间点（0w、1w、2w、6w）
  - condition：处理/饮食（NCD、MCD）
  - 组合分组（group）：推荐统一分组标签，如 0W_NCD、1W_MCD、2W_MCD、6W_MCD。若对象中尚无该列，可由脚本按 `"{timepoint}_{condition}"` 动态构建。
- 假设：上游合并/QC/去双细胞/重聚类/命名等已完成，且对象可直接用于下游。

## 3. 方法路线与判定标准（Pseudobulk + DESeq2）
- Pseudobulk：按“样本（=timepoint）× 细胞类型/簇”聚合计数，生成样本级 count 矩阵。
- DESeq2：以样本为单位拟合差异模型。
  - 设计：~ group（MCD vs NCD）；可选添加批次/线粒体比例等协变量（在代码中设置，不纳入 YAML）。
  - 归一化：DESeq2 sizeFactors。
  - 低表达过滤：如在 ≥ min_samples 个样本中 counts ≥ min_count。
- 显著性：padj ≤ p_adj_cutoff 且 |log2FC| ≥ lfc_threshold（默认 0.25）。
- 可选收缩：lfcShrink(type="apeglm"/"ashr")，用于更稳健的 log2FC（建议用于可视化）。

## 4. 目录结构与产物组织（分析与可视化分离）
在 `3_Analysis/2_DifferetialAnalysis/` 下组织如下（脚本运行后自动生成缺失目录）：
- scripts/
  - config.yaml（路径/列名映射/DE 参数/比较定义）
  - 1_run_deg_analysis.R（读取 config，遍历 comparisons，完成全部 DE 分析与数据导出）
  - viz/
    - 2_run_viz.R（读取分析阶段 CSV，批量生成火山图/MA 图/热图与（可选）富集）
- data/（理解性与中间数据）
  - metadata_snapshot.tsv（元数据列与唯一值快照）
  - counts_by_sample_cluster.tsv（样本×簇细胞计数；此处样本= timepoint）
  - cell_counts_threshold_failures.tsv（未满足阈值的组合）
- results/
  - tables/（DE 结果表）
    - NCD_0w_vs_MCD_1w/
      - clusterX.DESeq2.results.csv（gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj）
      - clusterX.DESeq2.results.shrunk.csv（可选：收缩后 lfc）
    - NCD_0w_vs_MCD_2w/
    - NCD_0w_vs_MCD_6w/
  - data/（供可视化/富集消费的CSV）
    - NCD_0w_vs_MCD_1w/
      - clusterX.normalized_counts.csv（样本×基因归一化计数或 rlog/VST）
      - clusterX.ranked_genes.csv（排序分数，默认使用 stat）
  - plots/（仅由 viz 阶段生成）
    - volcano/NCD_0w_vs_MCD_1w/clusterX.volcano.png
    - ma/NCD_0w_vs_MCD_1w/clusterX.MA.png
    - heatmap/NCD_0w_vs_MCD_1w/clusterX.topGenes.heatmap.png
    - enrichment/NCD_0w_vs_MCD_1w/clusterX.fgsea.hallmark.png
  - enrichment/（富集结果表）
    - NCD_0w_vs_MCD_1w/clusterX.fgsea.results.csv
- logs/
  - deg/1_run_deg_analysis.log
  - viz/2_run_viz.log
  - sessionInfo_YYYYMMDD-HHMMSS.txt
- reports/
  - run_summary.md（本次分析概览）
  - qc_summary.md（阈值命中与剔除统计）
  - methods.md（方法与参数记录）

命名规范：
- 比较ID：采用 `<control>_vs_<case>`，如 `NCD_0w_vs_MCD_1w`。
- 输出命名：`<stage>/<comparison>/<cluster>.<artifact>.<ext>`。

## 5. 配置文件（scripts/config.yaml）示例与说明
示例（建议作为起始模板）：
```yaml
paths:
  seurat_rds: "2_DataProcessing/RDS/nk.integrated.v4.rds"
  out_root: "3_Analysis/2_DifferetialAnalysis"

columns:
  sample_id: "timepoint"
  cell_type: "singleR.label"
  timepoint: "timepoint"
  condition: "condition"
  # 如果对象中已有组合分组列，如 "group" 存储 0W_NCD/1W_MCD 等，可直接设置：
  group: "group"
  # 若无组合列，脚本将按下面的格式动态生成（大小写需与元数据实际一致）
group_label_format: "{timepoint}_{condition}"   # 可选；当 columns.group 不存在时启用

deg:
  method: "DESeq2"
  lfc_threshold: 0.25
  p_adj_cutoff: 0.05
  pseudobulk:
    aggregation: "sum"
    min_cells_per_cluster: 50
  filter_expr:
    type: "counts"
    min_count: 10
    min_samples: 2
  shrinkage:
    enabled: true
    type: "apeglm"   # 可选 ashr

comparisons:
  - id: "NCD_0w_vs_MCD_1w"
    control: "0W_NCD"
    case: "1W_MCD"
  - id: "NCD_0w_vs_MCD_2w"
    control: "0W_NCD"
    case: "2W_MCD"
  - id: "NCD_0w_vs_MCD_6w"
    control: "0W_NCD"
    case: "6W_MCD"

viz:
  volcano:
    top_labels: 10
    label_by: "padj"     # 或 "pvalue"
  ma:
    use_shrunk_lfc: true
  heatmap:
    top_n: 30
    scale: "row"         # row / none
  enrichment:
    enabled: true
    databases: ["HALLMARK"]
    min_size: 10
    max_size: 500
    rank_metric: "stat"  # 用于 FGSEA 的排序分数
```

说明：
- paths/columns：输入 RDS 与元数据列映射。若 `columns.group` 缺失，将使用 `group_label_format` 由 `timepoint` 与 `condition` 组合生成。
- deg：DE 参数与阈值（可在项目中统一管理，确保可复现）。
- comparisons：定义全部比较，主脚本将逐一遍历，无需额外脚本。
- viz：仅被可视化脚本使用，均从 results/tables 与 results/data 读取 CSV。
- 降维参数（如 PCA/UMAP）不在此 YAML 中维护；相关旧代码将移除并按当前方案重写。

## 6. 运行顺序与命令示例
- 依赖（分析阶段）：Seurat、DESeq2、Matrix、dplyr、data.table、yaml
- 依赖（可视化阶段）：ggplot2、pheatmap/ComplexHeatmap、fgsea、msigdbr、dplyr、data.table、yaml

运行顺序：
1) 统一差异分析（一次性完成全部比较）
```
Rscript 3_Analysis/2_DifferetialAnalysis/scripts/1_run_deg_analysis.R \
  --config 3_Analysis/2_DifferetialAnalysis/scripts/config.yaml
```

2) 批量可视化与（可选）富集（消费 CSV）
```
Rscript 3_Analysis/2_DifferetialAnalysis/scripts/viz/2_run_viz.R \
  --config 3_Analysis/2_DifferetialAnalysis/scripts/config.yaml
```

说明：
- 运行后将生成/更新 logs 下的对应日志与 sessionInfo。
- 所有图与富集结果均由 viz 阶段生成；分析阶段不产出任何图。

## 7. 分析步骤与产物细节
A. 数据理解与阈值检查（由主脚本自动记录到 data/ 与 logs/）
- 快照：元数据列与唯一值（metadata_snapshot.tsv）。
- 统计：按 timepoint×condition×cell_type 的细胞数、每 timepoint 细胞数（将 timepoint 作为样本）、每簇样本覆盖度（counts_by_sample_cluster.tsv）。
- 阈值检查：`min_cells_per_cluster`、`filter_expr` 等命中/未命中组合（cell_counts_threshold_failures.tsv）。

B. 差异分析（Pseudobulk + DESeq2）
- 对每个比较、每个 cell_type：
  - 构建样本级计数矩阵（pseudobulk）。
  - 拟合 DESeq2 模型并导出 results CSV。
  - 可选导出收缩后 lfc 的 CSV。
  - 导出 normalized_counts.csv 与 ranked_genes.csv，供热图与 FGSEA 使用。

C. 可视化与富集（独立脚本）
- 火山图、MA 图、热图基于 CSV 生成，富集使用 msigdbr 集合与 FGSEA（输出 CSV + 图）。

## 8. 质量控制与判定细则
- 数据进入条件：
  - 每簇细胞数 ≥ min_cells_per_cluster。
  - 低表达过滤（counts-based）满足最小样本/计数阈值。
- 显著性：
  - padj ≤ p_adj_cutoff 且 |log2FoldChange| ≥ lfc_threshold。
- 稳健性：
  - 建议在可视化中使用收缩后的 lfc。

## 9. FAQ
- Q：某簇在某比较中样本不足怎么办？
  - A：在 cell_counts_threshold_failures.tsv 中记录，并在该簇跳过 DEG（在 qc_summary.md 汇总）。
- Q：热图/火山图分布异常？
  - A：检查预过滤阈值、sizeFactors 与 lfcShrink 设置；确认使用的 lfc 是否为收缩版本。
- Q：无显著通路富集？
  - A：调整 rank_metric 与数据库，或放宽/收紧阈值，并在 methods.md 中记录变更。

## 10. 版本与可复现性
- 每次运行写入 `logs/sessionInfo_*.txt`。
- 所有关键参数与比较列表在 `scripts/config.yaml` 固化。
- `reports/run_summary.md` 记录运行日期、输入 RDS、比较列表与主要统计概览。

## 10.5 迁移与清理
- 移除旧差异分析脚本：scripts/deg/2_deg_*.R 等模块化脚本不再使用；统一由 scripts/1_run_deg_analysis.R 执行。
- 删除/忽略任何与降维（PCA/UMAP）相关的 YAML 配置；降维相关代码不在本阶段维护，若后续可视化需要，在 viz/ 脚本内部定义。
- 如有旧的运行命令、Makefile 或 Runner 仍指向分散脚本，请同步替换为“1_run_deg_analysis.R + viz/2_run_viz.R”的新命令。

## 11. 变更记录
- v0.3：重构为单一脚本 + 配置驱动；比较由 comparisons 列表定义；以 timepoint 作为样本标识进行 pseudobulk；移除降维参数的 YAML 外置化并重写相关代码；可视化保持独立并消费 CSV（YYYY-MM-DD）。
- v0.2：DESeq2-only；分析与可视化分离；新增 MA 图。
- v0.1：定义初始目录结构、配置与运行规范。
