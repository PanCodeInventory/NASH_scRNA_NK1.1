# 任务指南文档

## 基本信息
- **任务名称**: 特定基因时间点表达分布（Lgals1）
- **创建时间**: 2025-11-13 13:25:00

## 任务目标
基于 NK 单细胞 Seurat 对象，展示并对比 Lgals1 在不同时间点（0W_NCD、1W_MCD、2W_MCD、6W_MCD）的表达分布，分别输出：
- 小提琴图（包含所有细胞）
- 小提琴图（仅包含表达>0的细胞）
- 每时间点的表达汇总表（均值/中位数/表达比例/细胞数）

### 实现步骤
1. **数据加载与校验**
   - 输入: `1_Files/RDS/nk.v5.rds`
   - 输出: 已校验的 Seurat 对象（含 `timepoint` 列与基因 `Lgals1`）
   - 工具/方法: Seurat v4/v5，yaml（配置）

2. **时间点顺序与分组设置**
   - 输入: `timepoint` 列与顺序 [0W_NCD, 1W_MCD, 2W_MCD, 6W_MCD]
   - 输出: 按顺序分组好的元数据
   - 工具/方法: dplyr，因子化排序

3. **可视化 - 小提琴图（两版）**
   - 输入: Seurat 对象与基因表达
   - 输出: 
     - `results/plots/Violin/Violin_Lgals1_byTimepoint_allcells.png/svg`
     - `results/plots/Violin/Violin_Lgals1_byTimepoint_expr_gt0.png/svg`
   - 工具/方法: Seurat::VlnPlot，ggplot2

4. **汇总统计表**
   - 输入: 表达向量与时间点分组
   - 输出: `results/tables/Lgals1_byTimepoint_summary.csv`
   - 指标: mean、median、pct_expr(Lgals1>0)、n_cells
   - 工具/方法: dplyr，readr

5. **日志与报告**
   - 输出: `logs/run_YYYYMMDD_HHMMSS.log`、`reports/process.md`、`reports/summary.md`

## 预期输出
- PNG/SVG 小提琴图（全体细胞、表达>0两版）
- Lgals1 在不同时间点的统计表（CSV）
- 执行过程与总结文档

## 技术路线
1. 使用 Seurat 读取 RDS 并提取基因表达（FetchData）
2. 使用 dplyr 按时间点计算统计指标
3. 使用 VlnPlot 生成可视化，ggplot2 调整主题与注释
4. 使用 yaml 管理参数与输出路径

## 备注
- 不排除任何 cluster；仅按 `timepoint` 分组展示
- 表达阈值用于“表达>0”视图，当前阈值为 0.0（可在 `scripts/config.yaml` 中调整）
