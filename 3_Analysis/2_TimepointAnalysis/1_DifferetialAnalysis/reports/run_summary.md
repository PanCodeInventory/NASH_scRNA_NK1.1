# 差异表达与可视化运行摘要（时间点驱动的整体与簇分析）

本摘要描述基于时间点（timepoint）进行样本比对的两套分析路径：
- 整体分析（overall）：不区分簇与细胞类型，直接在当前比较的所有细胞上进行伪样本聚合与差异分析。
- 簇分析（by_cluster）：按 Seurat 的 seurat_clusters 分层，在每个簇内分别进行伪样本聚合与差异分析（不使用 SingleR 的细胞类型标签）。

已完成的脚本内联更新与注释：
- 3_Analysis/2_DifferetialAnalysis/scripts/1_run_deg_analysis.R
  - 支持 cfg$deg$scope ∈ {overall, by_cluster, both}，默认 by_cluster。
  - 当没有条件列时，分组标签仅用 timepoint 构造（建议 config.yaml 中 group_label_format = "{timepoint}"）。
  - sample_id 使用 columns.sample_id（推荐设置为 "timepoint"），用于伪样本聚合（同一时间点内的细胞计数求和）。
  - 簇分析使用 columns.cluster（推荐新增为 "seurat_clusters"）作为分层依据，不再使用 columns.cell_type 做分层。
  - 整体分析输出以 all_cells.* 命名；簇分析输出以 <cluster>.* 命名。
  - 关键步骤（聚合、低表达过滤、DESeq2、log2FC 收缩、显著性判定）均加入中文注释。

- 3_Analysis/2_DifferetialAnalysis/scripts/viz/2_run_viz.R
  - 自动识别 *.DESeq2.results.csv 的前缀（含 all_cells 或具体簇号），逐个生成：火山图、MA 图、热图；当启用时，执行 FGSEA 富集并输出条形图与结果表。
  - 关键参数（top_labels、top_n、scale、rank_metric 等）读取自 config.yaml，并在绘图逻辑处加入中文注释。

配置文件（示例建议）
- scripts/config.yaml 中推荐设置（关键项）：
  - columns:
    - sample_id: "timepoint"
    - cluster: "seurat_clusters"  # 新增用于簇分析
    - timepoint: "timepoint"
    - condition:（可空）
    - group:（可空，缺失时自动用 group_label_format 构造）
  - group_label_format: "{timepoint}"
  - deg:
    - scope: "both"  # 同时做整体与簇分析；也可设为 "overall" 或 "by_cluster"
    - lfc_threshold: 0.25
    - p_adj_cutoff: 0.05
    - pseudobulk:
      - aggregation: "sum"
      - min_cells_per_cluster: 50  # 仅用于簇分析；整体分析不受此阈值限制
    - filter_expr:
      - type: "counts"
      - min_count: 10
      - min_samples: 2
    - shrinkage:
      - enabled: true
      - type: "apeglm"
  - comparisons:
    - id: "0W_vs_1W"; control: "0W"; case: "1W"
    - id: "0W_vs_2W"; control: "0W"; case: "2W"
    - id: "0W_vs_6W"; control: "0W"; case: "6W"
  - viz:
    - volcano: { top_labels: 10, label_by: "padj" }
    - ma: { use_shrunk_lfc: true }
    - heatmap: { top_n: 30, scale: "row" }
    - enrichment: { enabled: true, databases: ["HALLMARK"], min_size: 10, max_size: 500, rank_metric: "stat" }

运行方式
- 差异分析：
  - Rscript 3_Analysis/2_DifferetialAnalysis/scripts/1_run_deg_analysis.R --config 3_Analysis/2_DifferetialAnalysis/scripts/config.yaml
- 可视化：
  - Rscript 3_Analysis/2_DifferetialAnalysis/scripts/viz/2_run_viz.R --config 3_Analysis/2_DifferetialAnalysis/scripts/config.yaml

输出结构（位于 paths.out_root）
- results/tables/<cmp_id>/
  - all_cells.DESeq2.results.csv（整体分析，收缩版）
  - all_cells.DESeq2.results.raw.csv（整体分析，未收缩）
  - <cluster>.DESeq2.results.csv（簇分析，收缩版）
  - <cluster>.DESeq2.results.raw.csv（簇分析，未收缩）
- results/data/<cmp_id>/
  - all_cells.normalized_counts.csv / all_cells.ranked_genes.csv
  - <cluster>.normalized_counts.csv / <cluster>.ranked_genes.csv
- results/plots/
  - volcano/<cmp_id>/  ：all_cells.volcano.png 与 <cluster>.volcano.png
  - ma/<cmp_id>/        ：all_cells.MA.png 与 <cluster>.MA.png
  - heatmap/<cmp_id>/   ：all_cells.topGenes.heatmap.png 与 <cluster>.topGenes.heatmap.png
  - enrichment/<cmp_id>/：all_cells.fgsea.hallmark.png 与 <cluster>.fgsea.hallmark.png（启用时）
- results/enrichment/<cmp_id>/
  - all_cells.fgsea.results.csv 与 <cluster>.fgsea.results.csv（启用时）
- data/
  - metadata_snapshot.tsv（一次性快照，方便核查 timepoint 与簇）
  - cell_counts_threshold_failures.tsv（簇分析时低于最小细胞数的跳过项）
- logs/
  - deg/1_run_deg_analysis.log
  - viz/2_run_viz.log

注意事项
- comparisons 的标签必须与 group_label_format 构造出的标签一致（例如 timepoint 为 0W/1W/2W/6W）。若不一致，某比较可能被全部跳过。
- 整体分析不受 min_cells_per_cluster 限制，但仍需要每组至少有一个伪样本列（至少一个时间点对应 control 与 case）。
- 可视化脚本会在缺失 baseMean 时，以 |log2FC| 的秩替代横轴，以确保 MA 图仍可绘制。
- FGSEA 需要 ranked_genes.csv 的评分列（默认 "stat"），脚本会在缺失时以 sign(log2FC) * -log10(pvalue) 作为回退。
- 若编辑工具（lintr）提示 “unexpected end of input” 之类的告警，可先按上述命令直接运行脚本并查看日志；脚本结尾已包含会话信息输出与 sink 关闭语句，正常运行不受影响。
