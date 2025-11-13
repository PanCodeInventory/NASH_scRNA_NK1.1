# 任务指南文档

## 基本信息
- **任务名称**: 基因集分析 (Gene Set Analysis)
- **创建时间**: 2025-11-12 18:39:35

## 任务目标
使用 `geneset.txt` 文件中定义的多个基因集，对 `nk.integrated.v4.rds` 文件中不同细胞簇进行打分，并通过小提琴图展示各簇的特征得分情况，以判断不同簇的功能状态。

### 实现步骤
1. **步骤1**: 加载数据和基因集
   - 输入: 
     - `1_Files/RDS/nk.v5.rds` (Seurat 对象)
     - `3_Analysis/6_GeneSetAnalysis/geneset.txt` (基因集文件)
   - 输出: 加载到 R 环境中的 Seurat 对象和基因集列表
   - 工具/方法: R, `Seurat`, `readr`

2. **步骤2**: 计算基因集模块得分
   - 输入: Seurat 对象, 基因集列表
   - 操作: 使用 Seurat 的 `AddModuleScore` 函数为每个细胞计算每个基因集的得分
   - 输出: 更新后的 Seurat 对象，metadata 中包含新的得分列
   - 工具/方法: R, `Seurat::AddModuleScore`

3. **步骤3**: 可视化与结果导出
   - 输入: 包含模块得分的 Seurat 对象
   - 操作: 针对每个基因集得分，绘制按细胞簇分组的小提琴图（不显示细胞点）
   - 输出: 
     - `results/plots/` 目录下每个基因集的小提琴图 (e.g., `vlnplot_iNK_Markers.png`)
     - `results/tables/` 目录下包含各细胞得分的表格 (可选)
   - 工具/方法: R, `Seurat::VlnPlot`, `ggplot2`

## 预期输出
- **可视化图表**: `results/plots/` 目录下，每个基因集对应一张按簇分组的小提琴图。
- **分析总结报告**: `reports/summary.md` 中对各个簇在不同基因集上的得分情况进行解读。
- **（可选）数据表格**: `results/tables/module_scores.csv`，包含每个细胞在所有基因集上的得分。
