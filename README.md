# NK1.1 Pipeline

三阶段单细胞 RNA-seq 分析流程：
- Stage 1 (`1_Matrix2Markers`): 10x 读取、QC、整合、聚类、marker 分析
- Stage 2 (`2_ClusterEditor`): 聚类裁剪、重命名、重处理
- Stage 3 (`3_BasicViz_v2`): 结果可视化（UMAP、dotplot、river、比例图）

## 项目总览

```text
NK1.1/
├── 1_Matrix2Markers/
├── 2_ClusterEditor/
├── 3_BasicViz_v2/
└── Files/
```

核心代码集中在 `scripts/` / `workflow/`；`results/`、`logs/`、`Files/` 主要是数据与生成产物。

## 主数据文件概况

目标文件：`Files/merged.processed.pruned.renamed.h5ad`

| 项目 | 信息 |
|------|------|
| 维度 | `17550 x 21500` (`n_obs x n_vars`) |
| `obs` 关键列 | `sample`, `group`, `time`, `leiden`, `cluster_annot` |
| `obsm` | `X_pca (17550x50)`, `X_pca_harmony (17550x50)`, `X_umap (17550x2)` |
| `layers` | `counts (int64, 17550x21500)` |

`cluster_annot` 共 5 类：Cytotoxic NK、ILC1、Circulating NK、High-Ribosomal NK、Homeostatic NK。

## 各板块基础信息

| 板块 | 输入 | 输出 | 关键文件 | 运行目录 |
|------|------|------|----------|----------|
| `1_Matrix2Markers` | 10x 原始矩阵目录 | `merged.processed.h5ad` + QC/marker 结果 | `1_Matrix2Markers/scripts/Snakefile`, `1_Matrix2Markers/scripts/config.yaml` | `1_Matrix2Markers/scripts` |
| `2_ClusterEditor` | Stage 1 输出 h5ad | `merged.processed.pruned.renamed.h5ad` | `2_ClusterEditor/workflow/Snakefile`, `2_ClusterEditor/config/workflow_config.yaml` | `2_ClusterEditor/workflow` |
| `3_BasicViz_v2` | Stage 2 输出 h5ad | PDF 图件与配色文件 | `3_BasicViz_v2/main.R`, `3_BasicViz_v2/config.yaml` | `3_BasicViz_v2` |

## 快速运行

```bash
# Stage 1
cd "1_Matrix2Markers/scripts" && snakemake --cores 4

# Stage 2
cd "2_ClusterEditor/workflow" && snakemake --cores 4

# Stage 3
cd "3_BasicViz_v2" && Rscript main.R
```

## 使用约定

- 每个模块在自己的工作目录执行，不在仓库根目录直接跑 workflow。
- 参数优先由 YAML 配置驱动，避免在脚本内硬编码。
- `results/`、`logs/`、`.snakemake/`、`Files/` 属于生成/数据目录，不作为源码编辑目标。
