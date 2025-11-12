# 任务指南文档

## 基本信息
- **任务名称**: 单细胞聚类分析与差异基因富集分析
- **创建时间**: 2025-11-12 15:00:00
- **负责人**: Harry

## 任务目标
对 NASH 疾病模型单细胞数据进行聚类分析，识别不同细胞簇的标记基因，并进行 GSEA 富集分析，以理解不同细胞类型在疾病进展中的功能变化。

### 实现步骤

#### 阶段 1: 数据准备与标记基因鉴定 ✅ 已完成

1. **步骤1: 数据加载与预处理**
   - 输入: `2_DataProcessing/3_UMAP-Tuning/data/nk_noCluster6_final.rds`
   - 输出: 清洗后的 Seurat 对象
   - 工具/方法: Seurat v4/v5
   - 状态: ✅ 完成

2. **步骤2: 聚类分析与可视化**
   - 输入: Seurat 对象
   - 输出: 
     - UMAP 聚类图
     - 簇组成堆叠图
     - 簇比例线图
   - 工具/方法: Seurat 聚类算法, UMAP 降维
   - 状态: ✅ 完成
   - 结果: 识别出 7 个簇（Cluster 0-6）

3. **步骤3: 标记基因鉴定 (FindMarkers)**
   - 输入: 聚类后的 Seurat 对象
   - 输出: 
     - `results/tables/markers_all_clusters.csv` (4,286 个标记基因)
     - `results/tables/markers_top10_per_cluster.csv` (每簇 Top10)
     - `results/plots/markers/dotplot_markers_top10.png`
   - 工具/方法: 
     - Seurat::FindAllMarkers()
     - 统计方法: Wilcoxon rank sum test
     - 筛选阈值: p_val_adj < 0.05, avg_log2FC > 0
   - 状态: ✅ 完成

#### 阶段 2: 富集分析 🔄 进行中

4. **步骤4a: GO/KEGG 富集分析 (ORA)**
   - 输入: `results/tables/markers_all_clusters.csv`
   - 输出: 
     - 每个簇的 GO BP/CC/MF 富集结果
     - 每个簇的 KEGG 通路富集结果
     - 汇总气泡图和表格
   - 工具/方法: 
     - clusterProfiler::enrichGO()
     - clusterProfiler::enrichKEGG()
     - 物种: Mus musculus (org.Mm.eg.db)
     - 筛选: p.adjust < 0.05, q-value < 0.2
     - 排除: Cluster 6, 线粒体/核糖体相关通路
   - 脚本: `scripts/run_enrichment_by_cluster.R`
   - 状态: 放弃

5. **步骤4b: GO BP GSEA 分析** ⭐ 调整
   - **目标**: 对每个簇进行基因集富集分析 (GSEA)，以识别在每个簇中显著上调或下调的生物学过程。
   - **输入**: 
     - `results/tables/markers_all_clusters.csv` (包含所有基因的 log2FC 值)
     - GO BP (Gene Ontology Biological Process) 基因集
   - **输出**:
     - 每个簇的 GSEA 富集结果表格 (包含 NES, p.adjust, q-value)
     - 每个簇的 GSEA enrichment plot (Top 5-10 上调/下调通路)
     - 跨簇比较的 GSEA 结果热图 (按 NES 值)
     - 跨簇比较的气泡图
   - **工具/方法**:
     - **方法**: `clusterProfiler::gseGO()`
     - **基因集**: `org.Mm.eg.db` (Mus musculus)
     - **排序指标**: `avg_log2FC` (从 FindMarkers 结果中提取，并按簇分组)
     - **统计参数**:
       - nPerm = 10000 (置换检验次数)
       - minGSSize = 10, maxGSSize = 500 (基因集大小)
       - pvalueCutoff = 0.05
       - pAdjustMethod = "BH"
     - **可视化**:
       - `enrichplot::gseaplot2()` - 单通路富集图
       - `pheatmap::pheatmap()` - NES 热图
       - `ggplot2` - 气泡图和条形图
   - **脚本**: `scripts/run_gsea_gobp_by_cluster.R` (待创建)
   - **状态**: ⏳ 待执行

6. **步骤4c: 富集结果整合与比较**
   - 输入: 
     - GO/KEGG ORA 结果
     - HALLMARK GSEA 结果
   - 输出:
     - 整合的富集结果表格
     - 跨方法比较图表
     - 关键通路汇总报告
   - 工具/方法: dplyr, ggplot2
   - 状态: 放弃

#### 阶段 3: 高级可视化与报告

7. **步骤5: 综合可视化**
   - 输入: 所有富集分析结果
   - 输出:
     - 多面板组合图
     - 通路-基因网络图
     - 交互式可视化（可选）
   - 工具/方法: 
     - cowplot/patchwork - 组合图
     - enrichplot::cnetplot() - 网络图
     - plotly (可选) - 交互图
   - 状态: ⏳ 待执行

8. **步骤6: 生成分析报告**
   - 输入: 所有分析结果
   - 输出: 
     - `reports/summary.md` - 完整分析报告
     - `results/README.md` - 结果说明文档
   - 状态: ⏳ 待执行

## 预期输出
- 聚类分析结果和可视化图表
- 各细胞簇的标记基因列表
- GSEA 富集分析结果
- 富集分析可视化图表
- 完整的分析报告

## 技术路线
1. 使用 Seurat 进行标准单细胞分析流程
2. 采用 HALLMARK 基因集进行 GSEA 分析
3. 使用 Mus musculus 物种注释
4. 设置 min_size=10, max_size=500 的富集阈值
