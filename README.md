# scRNA-seq NASH NK细胞分析项目

本仓库用于在小鼠 NASH（非酒精性脂肪性肝炎）疾病模型中，围绕 NK 细胞进行单细胞 RNA 测序（scRNA-seq）分析的全流程管理与实现。包含数据组织、主数据处理管线（清理、调参、UMAP/聚类）、以及下游分析（簇比例、差异表达、基因集与特异基因分析）等。

## 项目概览

- 研究对象：NK 细胞在 NASH 进程中（0W、1W、2W、6W）的时序变化
- 核心对象：若干经过清理、注释、调参与重命名的 Seurat RDS 对象（例如 `nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds` 等）
- 统一的项目管理文档体系：每个模块均包含 `reports/guidedoc.md`、`reports/process.md`、`reports/summary.md`，确保可追溯、可复现

## 目录结构（摘要）

```
scRNA-seq/
├─ README.md
├─ .gitignore
├─ 1_Files/                          # 原始/预处理数据（分组）
│  ├─ Bam/
│  │  ├─ MCD-1w/ MCD-2w/ MCD-6w/ NCD/
│  ├─ Matrix/
│  │  ├─ CD45.2/                     # CD45.2 分选数据（含各时间点）
│  │  └─ NK1.1/                      # NK1.1 分选数据（含各时间点、RDS）
│  └─ RDS/                           # 原始或中间 RDS 对象（如有）
│
├─ 2_DataProcessing/                 # 主数据处理管线（脚本与产物）
│  ├─ 1_Samples_Merging/
│  │  ├─ guidedoc.md / guidedoc_*.md
│  │  └─ Results/
│  │     ├─ data/ plots/ rds/
│  ├─ 2_Doublet_Removed/
│  │  ├─ plots/                      # 清理/注释后 UMAP 等图件
│  │  └─ reports/                    # 清理与注释修复等报告
│  ├─ 3_UMAP-Tuning/
│  │  ├─ data/ logs/ plots/          # 调参与指标、会话信息与图件
│  ├─ reports/
│  │  ├─ analysis_files_update_report.md
│  │  ├─ cluster_renaming_*.md       # 簇重命名策略与执行报告
│  │  └─ visualization_report_*.md   # 可视化报告（含隐藏簇策略）
│  └─ scripts/
│     ├─ OVERVIEW.md
│     ├─ Scripts_New/                # 新版统一管线
│     │  ├─ run_pipeline.R
│     │  ├─ config/                  # cluster_mapping.yaml, parameters.yaml, paths.yaml
│     │  ├─ stages/                  # 01~06 分阶段脚本
│     │  │  ├─ 01_data_import.R
│     │  │  ├─ 02_basic_qc.R
│     │  │  ├─ 03_cell_annotation.R
│     │  │  ├─ 04_cell_filtering.R
│     │  │  ├─ 05_parameter_tuning.R
│     │  │  └─ 06_visualization.R
│     │  └─ utils/                   # 复用函数（plotting/seurat/validation）
│     ├─ Scripts_Backup/
│     │  ├─ original/                # 历史脚本归档（remove_* / tune_* 等）
│     │  └─ deprecated/
│     └─ Tailor_scripts/             # 定制化脚本（如过滤特定簇等）
│
└─ 3_Analysis/                       # 下游分析模块
   ├─ 1_ClusterAnalysis/
   │  └─ Enrichment/
   │     ├─ scripts/                 # 比例导出、标记基因等
   │     ├─ results/
   │     │  ├─ enrichment/ enrichment_gsea_gobp/
   │     │  ├─ plots/ markers/
   │     │  └─ tables/
   │     ├─ logs/
   │     └─ reports/                 # guidedoc/process/summary
   ├─ 2_DifferetialAnalysis/
   │  ├─ scripts/ config.yaml
   │  ├─ results/ plots/ tables/
   │  ├─ logs/
   │  └─ reports/                    # run_summary.md 等
   └─ 3_GeneAnalysis/
      ├─ 1_SpecificGene/             # 单基因分布与变化
      │  ├─ scripts/ config.yaml
      │  ├─ results/ plots/
      │  └─ reports/                 # guidedoc/process/summary
      └─ 2_GeneSetAnalysis/          # 基因集打分
         ├─ scripts/
         │  ├─ stages/               # 01_score_genesets.R
         │  ├─ utils/ config/ Scripts_Backup/
         ├─ results/ data/ plots/ tables/
         ├─ logs/
         └─ reports/                 # guidedoc/process/summary
```

> 注：项目结构清单为当前仓库中可见目录的聚合摘要，若有新模块（如 `3_Analysis/2_TimepointAnalysis/...`）在后续执行中将继续补充。

## 快速开始

1) 准备运行环境  
- R（建议 4.x）、Seurat v5、SingleR、scDblFinder、UCell、tidyverse 等。  
- 依赖以脚本为准，若缺包请按报错安装或在 `scripts/Scripts_New/utils` 中查看调用。

2) 配置参数  
- 编辑 `2_DataProcessing/scripts/Scripts_New/config/paths.yaml`、`parameters.yaml`、`cluster_mapping.yaml`。
- 路径与参数需与本机数据位置一致（如 `1_Files/Matrix/...`）。

3) 运行主管线（示例）  
```bash
Rscript 2_DataProcessing/scripts/Scripts_New/run_pipeline.R
```
- 阶段性产物与日志将输出到对应的 `2_Doublet_Removed/`、`3_UMAP-Tuning/`、以及 `reports/` 目录。
- 完整执行记录请参考各模块的 `reports/process.md`。

## 下游分析入口（示例）

- 簇比例与标记基因（ClusterAnalysis/Enrichment）
  - `3_Analysis/1_ClusterAnalysis/Enrichment/scripts/export_cluster_proportions.R`
  - `3_Analysis/1_ClusterAnalysis/Enrichment/scripts/find_markers_simple.R`
  - 产出：`results/tables/`、`results/plots/` 等

- 差异表达分析（DifferetialAnalysis）
  - `3_Analysis/2_DifferetialAnalysis/scripts/1_run_deg_analysis.R`
  - 可视化：`scripts/viz/2_run_viz.R` 与 `scripts/viz/run_module_viz.R`
  - 产出：`results/plots/`（volcano/heatmap/enrichment 等）、`results/tables/`

- 基因集分析（GeneSetAnalysis）
  - `3_Analysis/3_GeneAnalysis/2_GeneSetAnalysis/scripts/stages/01_score_genesets.R`
  - 产出：`results/plots/`（小提琴图等）、`results/data/`、`results/tables/`

- 特异基因分析（SpecificGene）
  - `3_Analysis/3_GeneAnalysis/1_SpecificGene/scripts/plot_specific_gene_variation.R`
  - 产出：`results/plots/` 与汇总 `results/tables/`

> 所有模块均配套 `reports/guidedoc.md`（目标与方案）、`reports/process.md`（执行记录与 Todo）、`reports/summary.md`（结果总结）。

## 关键对象与报告

- 清理与注释修复报告：
  - `2_DataProcessing/2_Doublet_Removed/reports/cleaning_report.md`
  - `2_DataProcessing/2_Doublet_Removed/reports/singleR_fix_report.md`
- 簇重命名与隐藏策略：
  - `2_DataProcessing/reports/cluster_renaming_final_report.md`
  - `2_DataProcessing/reports/visualization_report_hidden_cluster6.md`
- 最终/阶段性对象：
  - `2_DataProcessing/2_Doublet_Removed/RDS/*`（如有）
  - `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds`（如有）

## 近期更新（摘要）

- 2025-11-13 时间点特异基因表达分布（Lgals1）
  - 模块：`3_Analysis/2_TimepointAnalysis/2_SpecificGeneVariation/`（包含两版小提琴图与汇总表）
  - 文档：`reports/guidedoc.md`、`reports/process.md`（样式更新：叠加黑色小点）、`reports/summary.md`

- 2025-11-12 基因集打分（iNK, TR-NK, CD56bright-like, CD56dim-like）
  - 模块：`3_Analysis/3_GeneAnalysis/2_GeneSetAnalysis/`
  - 脚本：`scripts/stages/01_score_genesets.R`
  - 产出：`results/plots/` 各基因集小提琴图与完整文档集

- 2025-10-23 簇重命名与隐藏策略 + 重新调优
  - 移除污染簇、最佳参数（示例：dims=10、res=0.3）、获得具生物学意义的 NK 子群
  - 报告与图件已系统化存档至 `2_DataProcessing/reports/` 与 `3_Analysis/1_ClusterAnalysis/Enrichment/results/*`

- 2025-10-20~21 清理、SingleR 修复与全流程重跑
  - 集成 scDblFinder 去双胞、SingleR 自动注释、UCell 评分等
  - 规范化 SCE/logcounts 流程，稳健处理元数据对齐与日志/报告目录生成
  - 历史脚本已归档至 `2_DataProcessing/scripts/Scripts_Backup/original/`

> 更详细的逐日变更与产出示例请参考各子模块内的 `reports/*` 与 `results/*`。

## 贡献与协作

- 提交规范：建议使用约定式提交（如 `feat:`、`fix:`、`docs:`、`refactor:` 等），并在相关模块的 `reports/process.md` 中同步更新执行记录。
- 问题与建议：请在对应模块的 `reports/process.md` 中记录，或发起 Issue。

## 远程仓库

- origin: `git@github.com:PanCodeInventory/scRNA-seq-NASH-NK-analysis.git`
