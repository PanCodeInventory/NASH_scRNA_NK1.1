# 4_GenesetScore 技术讨论稿（先讨论，不执行）

## 1) 目标与边界

- 目标：基于基因集对每个 `leiden` 簇进行打分，描述簇的潜在 NK 亚型特征。
- 当前阶段：仅确定技术路线与参数，不跑全流程。
- 输入数据：`/home/user/Pan Chongshi/Projects/NASH/NK1.1/Files/merged.processed.pruned.renamed.h5ad`
- 输出目录：`4_GenesetScore/geneset_scoring/`

## 2) 现有数据与代码约束（已核对）

- `adata.obs` 已有 `leiden` 和 `cluster_annot`。
- `adata.X` 为 normalized + log1p 表达矩阵。
- `adata.layers['counts']` 为原始计数。
- 注意：本仓库里 `adata.raw` 不是原始 counts，而是 normalize/log1p 后设置的快照（不建议把 raw 当 counts 用）。
- 仓库目前没有现成 geneset score 实现，`4_GenesetScore/` 仅有规划文档。

## 3) 方法选择建议（默认方案）

### 默认方案：Scanpy 双通道打分 + 差值

对每个亚型签名分别计算：

1. `pos_score`：阳性 marker 基因集打分
2. `neg_score`：阴性 marker 基因集打分（若存在）
3. `final_score = Z(pos_score) - Z(neg_score)`（若无阴性集则 `final_score = Z(pos_score)`）

说明：

- 优点：与当前 Scanpy 管线兼容，不新增依赖，易维护。
- 与你的设想一致（`Z(阳性)-Z(阴性)`）。
- 细胞级评分后，再做簇级聚合，能保留簇内异质性信息。

### 备选方案：UCell/pyUCell

- 优点：对稀疏矩阵更稳，天然支持正负签名（rank-based）。
- 缺点：引入新依赖，和现有流水线耦合稍弱。
- 建议：若后续发现小基因集或掉零导致 Scanpy 分数不稳，再切换。

## 4) 表达矩阵与分组策略

- 评分矩阵默认使用 `adata.X`（log1p normalized），与当前项目可视化习惯一致。
- 分组键使用 `leiden`（用户目标是“每个 leiden 簇”）。
- `cluster_annot` 作为结果对照列，用于 sanity check，不作为主分组键。

## 5) 基因集定义与清洗规则

建议新增一个结构化配置文件（如 `genesets.yaml`），每个亚型包含：

- `name`
- `positive` 基因列表
- `negative` 基因列表（可空）

清洗规则：

1. 去重、去空值。
2. 基因名匹配采用多级回退：`exact -> upper -> title`。
3. 记录每个基因集的匹配率（匹配到的基因数 / 输入基因数）。
4. 低匹配率告警：例如 `< 30%` 标记为低可信。

## 6) 打分与标准化细节

对每个 subtype：

1. 用 `scanpy.tl.score_genes` 计算 `pos_score`。
2. 若存在阴性集，计算 `neg_score`。
3. 对 `pos_score`、`neg_score` 分别在全体细胞做 Z-score。
4. 计算 `final_score = z_pos - z_neg`。

簇级统计（按 `leiden`）：

- `mean_final`
- `median_final`（推荐作为主排序依据，抗离群点）
- `std_final`
- `n_cells`

## 7) 簇注释判定逻辑（建议）

对每个 `leiden` 簇：

1. 取 `median_final` 最高的 subtype 作为候选注释。
2. 计算置信度：`delta = top1 - top2`。
3. 若 `delta` 过小（阈值可配，如 `< 0.2`），标记为 `Ambiguous`。
4. 输出“候选注释 + 置信度 + 次优候选”。

## 8) 计划输出文件（落在 4_GenesetScore/）

建议后续实现时生成：

- `config/genesets.yaml`：基因集配置
- `scripts/run_geneset_score.py`：主脚本
- `results/cell_scores.tsv.gz`：细胞级分数（可压缩）
- `results/cluster_scores_by_leiden.tsv`：簇级统计表
- `results/cluster_annotation_candidates.tsv`：簇注释候选与置信度
- `results/adata.geneset_scored.h5ad`：带分数字段的 AnnData
- `results/plots/`：热图、UMAP score overlay

## 9) 关键技术分歧点（建议我们先确认）

以下是会显著影响结果解释的参数点：

1. **主方法**：先用 Scanpy（默认）还是直接上 UCell。
2. **主排序统计量**：`median_final`（默认）还是 `mean_final`。
3. **Ambiguous 阈值**：`top1-top2` 的阈值设为多少（默认建议 0.2）。
4. **最小匹配基因数**：例如每个基因集至少匹配 3 个基因才纳入判定。

## 10) 推荐的下一步（讨论后执行）

1. 固定参数：方法、阈值、最小匹配基因数。
2. 将 `geneset_info.md` 中的 CIMA 表转成结构化 `genesets.yaml`。
3. 实现 `run_geneset_score.py` 并先做小规模 dry-run（不改现有主流程）。
4. 产出簇级评分表和候选注释表，再与你一起复核。
