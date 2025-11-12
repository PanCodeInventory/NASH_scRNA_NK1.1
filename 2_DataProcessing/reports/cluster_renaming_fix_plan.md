# 簇重命名修复计划

## 问题诊断

### 当前状态
- **标记基因数据显示**：簇5仍然包含B细胞标记基因（Iglc3、Cd79a等）和增殖基因（H2afx、Mki67等）
- **用户反馈**：原cluster 5（B细胞污染）仍然以cluster 5显示在UMAP图中，没有变成cluster 6
- **根本问题**：簇重命名没有正确执行，或者分析文件没有正确更新

### 预期状态
- **簇5**：应该只包含增殖相关基因（H2afx、Mki67、Hist1h3c等）
- **簇6**：应该包含B细胞污染基因（Iglc3、Cd79a、Ebf1等），并在可视化中隐藏
- **UMAP显示**：原增殖细胞应该显示为簇5，原B细胞污染应该显示为簇6（或隐藏）

## 修复策略

### 第一步：验证当前RDS对象状态
需要检查以下RDS文件的实际簇分配：
- `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds`

### 第二步：正确的簇重命名映射
```r
# 正确的映射关系
cluster_mapping <- c(
  "0" = "0",
  "1" = "1", 
  "2" = "2",
  "3" = "3",
  "4" = "4",
  "5" = "6",  # 原簇5（B细胞污染）→ 簇6
  "6" = "5"   # 原簇6（增殖）→ 簇5
)
```

### 第三步：修复RDS对象
1. 读取原始RDS对象
2. 检查当前的seurat_clusters分配
3. 应用正确的簇重命名映射
4. 验证重命名后的簇分配
5. 保存修复后的RDS对象

### 第四步：修复分析文件
需要修复以下文件：
- `3_Analysis/1.ClusterAnalysis/data/markers_all_clusters.csv`
- `3_Analysis/1.ClusterAnalysis/data/markers_top10_per_cluster.csv`
- `3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv`

### 第五步：重新生成标记基因分析
基于修复后的RDS对象，重新运行标记基因分析：
```bash
Rscript 3_Analysis/Scripts/find_markers_simple.R
```

### 第六步：重新生成可视化
生成正确的UMAP图件，显示重命名后的簇分配：
- 原增殖细胞（现为簇5）应该在UMAP中正确显示
- 原B细胞污染（现为簇6）应该在UMAP中隐藏

## 技术实现细节

### R脚本结构
```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# 1. 读取RDS对象
rds_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds"
obj <- readRDS(rds_path)

# 2. 检查当前簇分配
current_clusters <- obj$seurat_clusters
table(current_clusters)

# 3. 应用簇重命名映射
cluster_mapping <- c("0"="0", "1"="1", "2"="2", "3"="3", "4"="4", "5"="6", "6"="5")
obj$seurat_clusters <- factor(cluster_mapping[as.character(obj$seurat_clusters)])

# 4. 验证重命名结果
table(obj$seurat_clusters)

# 5. 保存修复后的对象
output_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.fixed.rds"
saveRDS(obj, file = output_path)

# 6. 更新find_markers_simple.R脚本的RDS路径
# 将rds_path更新为修复后的文件路径
```

### 关键检查点
1. **簇数量验证**：确保重命名后仍然是7个簇（0-6）
2. **细胞数量验证**：确保总细胞数保持不变（18,530个）
3. **标记基因验证**：确保簇5和簇6的标记基因符合预期
4. **可视化验证**：确保UMAP图正确显示重命名后的簇

## 预期结果

### 修复后的标记基因分布
- **簇5**：H2afx、Cdca3、Hist1h3c、Mki67、Ezh2、Tuba1c（增殖标记）
- **簇6**：Iglc3、Cd79a、Iglc2、Ebf1、Ms4a1、Cd79b（B细胞标记，将被隐藏）

### 修复后的可视化
- UMAP图中增殖细胞显示为簇5
- B细胞污染显示为簇6（或在可视化中隐藏）
- 簇比例分析图表正确反映重命名后的分布

## 执行建议

### 立即行动
1. 切换到Code模式
2. 执行修复脚本
3. 验证修复结果
4. 更新所有相关分析文件

### 质量保证
1. 在修复前备份所有原始文件
2. 逐步验证每个修复步骤
3. 生成详细的修复报告
4. 更新Git版本控制

## 文件清单

### 需要创建的文件
- `2_DataProcessing/Scripts/fix_cluster_renaming.R` - 修复脚本
- `2_DataProcessing/reports/cluster_renaming_fix_report.md` - 修复报告

### 需要更新的文件
- `3_Analysis/Scripts/find_markers_simple.R` - 更新RDS路径
- 所有分析数据文件（markers_*.csv, cluster_proportions_*.csv）
- 所有可视化图件

### 需要验证的文件
- `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.fixed.rds`
- `3_Analysis/1.ClusterAnalysis/data/markers_top10_per_cluster.csv`

---

**创建时间**：2025-10-23  
**优先级**：高  
**状态**：待执行  
**建议执行模式**：Code