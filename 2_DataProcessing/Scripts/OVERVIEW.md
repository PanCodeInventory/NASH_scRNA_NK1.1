# scRNA-seq 数据处理脚本系统综述

## 📋 项目概述

本项目将原有的15个分散、重复的R脚本重构为一个模块化、配置驱动的scRNA-seq数据处理流程。新系统解决了硬编码路径、代码重复、缺乏统一配置等问题，提供了更好的可维护性和用户体验。

## 🏗️ 代码组织架构

### 目录结构
```
2_DataProcessing/Scripts/
├── Scripts_New/                    # 新脚本系统
│   ├── config/                     # 配置文件目录
│   │   ├── paths.yaml             # 数据路径配置
│   │   ├── parameters.yaml        # 分析参数配置
│   │   └── cluster_mapping.yaml  # 簇映射配置
│   ├── utils/                      # 工具函数库
│   │   ├── seurat_utils.R         # Seurat操作工具函数
│   │   ├── plotting_utils.R       # 绘图工具函数
│   │   └── validation_utils.R     # 验证工具函数
│   ├── stages/                     # 分阶段脚本
│   │   ├── 01_data_import.R        # 数据导入与预处理
│   │   ├── 02_basic_qc.R          # 基础质量控制
│   │   ├── 03_cell_annotation.R    # 细胞标注
│   │   ├── 04_cell_filtering.R     # 细胞过滤
│   │   ├── 05_parameter_tuning.R   # 参数优化
│   │   ├── 06_visualization.R       # 可视化
│   │   └── run_pipeline.R          # 主控制脚本
│   ├── README.md                   # 项目说明
│   ├── DEVELOPMENT_PLAN.md         # 开发计划
│   └── SCRIPTS_SUMMARY.md         # 脚本功能总结
├── Scripts_Backup/                 # 原始脚本备份
└── OVERVIEW.md                    # 本文档
```

### 架构设计原则

#### 1. 模块化设计
- **阶段分离**：每个处理阶段独立为单个脚本
- **功能解耦**：工具函数与业务逻辑分离
- **接口标准化**：统一的输入输出格式

#### 2. 配置驱动
- **路径配置**：所有文件路径通过`paths.yaml`管理
- **参数配置**：分析参数通过`parameters.yaml`管理
- **映射配置**：簇重命名通过`cluster_mapping.yaml`管理

#### 3. 可复用性
- **工具函数库**：775行可复用工具代码
- **通用接口**：支持不同数据集和参数配置
- **扩展性**：易于添加新的分析功能

## 🚀 如何使用

### 快速开始

#### 1. 运行完整流程
```bash
cd 2_DataProcessing/Scripts_New
Rscript run_pipeline.R
```

#### 2. 运行特定阶段
```bash
# 只运行数据导入和基础QC
Rscript run_pipeline.R --stages data_import,basic_qc

# 从特定阶段恢复执行
Rscript run_pipeline.R --resume cell_annotation
```

#### 3. 试运行模式
```bash
# 检查配置和依赖，不执行实际分析
Rscript run_pipeline.R --dry-run
```

#### 4. 单独运行阶段脚本
```bash
# 使用默认配置运行数据导入
Rscript 01_data_import.R

# 使用自定义配置
Rscript 01_data_import.R custom_params.yaml custom_paths.yaml
```

### 配置文件使用

#### 路径配置 (`config/paths.yaml`)
```yaml
# 输入数据路径
input:
  base_dir: "/home/harry/NASH/scRNA-seq/Files/Filter Files"
  samples:
    - "NCD_NK1.1"
    - "MCD-1W_NK1.1"
    - "MCD-2W_NK1.1"
    - "MCD-6W_NK1.1"

# 输出路径
output:
  base_dir: "/home/harry/NASH/scRNA-seq/2_DataProcessing"
  rds_dir: "RDS"
  plots_dir: "plots"
  reports_dir: "reports"
```

#### 参数配置 (`config/parameters.yaml`)
```yaml
# 预处理参数
preprocessing:
  nfeatures: 3000
  normalization_method: "SCT"
  integration_method: "SCT"

# 聚类参数
clustering:
  resolution: 0.3
  dims: 1:10
  algorithm: "louvain"
```

#### 簇映射配置 (`config/cluster_mapping.yaml`)
```yaml
# 簇重命名映射
cluster_renaming:
  "5": "6"  # 原簇5(B细胞污染) → 簇6
  "6": "5"  # 原簇6(增殖) → 簇5

# 隐藏簇配置
hidden_clusters:
  - "6"  # 隐藏B细胞污染簇
```

## 🔧 特定板块调整指南

### 1. 数据导入阶段调整

#### 定位相关文件
- **主脚本**：`01_data_import.R`
- **工具函数**：`utils/seurat_utils.R` 中的 `read_10x_data_safe()`, `create_seurat_objects()`, `integrate_data()`
- **配置文件**：`config/paths.yaml`, `config/parameters.yaml`

#### 常见调整场景
```bash
# 修改输入数据路径
# 编辑 config/paths.yaml 中的 input.base_dir 和 samples

# 调整SCTransform参数
# 编辑 config/parameters.yaml 中的 preprocessing 部分

# 修改特征选择数量
# 编辑 config/parameters.yaml 中的 preprocessing.nfeatures
```

### 2. 质量控制阶段调整

#### 定位相关文件
- **主脚本**：`02_basic_qc.R`
- **工具函数**：`utils/seurat_utils.R` 中的 `detect_doublets()`, `remove_doublets()`
- **配置文件**：`config/parameters.yaml` 中的 quality_control 部分

#### 常见调整场景
```bash
# 调整双胞检测参数
# 编辑 config/parameters.yaml 中的 quality_control.doublet_detection

# 修改双胞阈值
# 编辑 config/parameters.yaml 中的 quality_control.doublet_threshold
```

### 3. 细胞注释阶段调整

#### 定位相关文件
- **主脚本**：`03_cell_annotation.R`
- **工具函数**：`utils/seurat_utils.R` 中的 `annotate_cells_singler()`
- **配置文件**：`config/parameters.yaml` 中的 annotation 部分

#### 常见调整场景
```bash
# 更改参考数据集
# 编辑 config/parameters.yaml 中的 annotation.reference_dataset

# 调整注释阈值
# 编辑 config/parameters.yaml 中的 annotation.pruning_threshold
```

### 4. 细胞过滤阶段调整

#### 定位相关文件
- **主脚本**：`04_cell_filtering.R`
- **工具函数**：`utils/seurat_utils.R` 中的 `filter_cells_by_type()`
- **配置文件**：`config/cluster_mapping.yaml`, `config/parameters.yaml`

#### 常见调整场景
```bash
# 修改要过滤的细胞类型
# 编辑 config/cluster_mapping.yaml 中的 cell_type_filtering

# 调整过滤规则
# 编辑 config/parameters.yaml 中的 cell_filtering 部分
```

### 5. 参数优化阶段调整

#### 定位相关文件
- **主脚本**：`05_parameter_tuning.R`
- **工具函数**：`utils/seurat_utils.R` 中的 `tune_parameters_grid()`
- **配置文件**：`config/parameters.yaml` 中的 parameter_tuning 部分

#### 常见调整场景
```bash
# 调整参数搜索范围
# 编辑 config/parameters.yaml 中的 parameter_tuning.dims_grid 和 resolution_grid

# 修改评估指标
# 编辑 config/parameters.yaml 中的 parameter_tuning.evaluation_metrics
```

### 6. 可视化阶段调整

#### 定位相关文件
- **主脚本**：`06_visualization.R`
- **工具函数**：`utils/plotting_utils.R` 中的所有绘图函数
- **配置文件**：`config/parameters.yaml` 中的 visualization 部分

#### 常见调整场景
```bash
# 修改图片格式和分辨率
# 编辑 config/parameters.yaml 中的 visualization.formats 和 dpi

# 调整颜色方案
# 编辑 config/parameters.yaml 中的 visualization.color_palettes

# 自定义图件布局
# 编辑 config/parameters.yaml 中的 visualization.layout
```

## 📊 系统特性

### 1. 错误处理机制
- **完善的try-catch**：每个主要操作都有异常处理
- **详细日志**：结构化的错误信息和进度跟踪
- **优雅失败**：错误时保存中间结果，支持断点恢复

### 2. 性能监控
- **内存管理**：大文件分块处理，避免内存溢出
- **进度显示**：实时显示处理进度和预计完成时间
- **资源监控**：记录CPU和内存使用情况

### 3. 可扩展性
- **插件架构**：易于添加新的分析模块
- **配置模板**：支持不同项目的配置模板
- **版本控制**：配置文件版本管理，确保兼容性

## 🔍 故障排除

### 常见问题及解决方案

#### 1. 路径问题
```bash
# 问题：找不到输入文件
# 解决：检查 config/paths.yaml 中的路径配置
# 确保使用绝对路径或相对于项目根目录的路径
```

#### 2. 内存问题
```bash
# 问题：内存不足
# 解决：调整 config/parameters.yaml 中的 computational.chunk_size
# 减少并行处理的线程数
```

#### 3. 依赖包问题
```bash
# 问题：缺少R包
# 解决：运行以下命令安装依赖
Rscript -e "install.packages(c('Seurat', 'SingleR', 'scDblFinder', 'yaml'))"
```

#### 4. 配置文件问题
```bash
# 问题：配置文件格式错误
# 解决：使用YAML验证工具检查语法
# 参考配置文件模板确保格式正确
```

## 📈 性能优化建议

### 1. 大数据集优化
- **分块处理**：启用数据分块处理模式
- **并行计算**：利用多核CPU加速计算
- **内存映射**：使用内存映射文件减少内存占用

### 2. 参数调优优化
- **智能搜索**：使用贝叶斯优化替代网格搜索
- **早停机制**：设置性能阈值提前终止搜索
- **缓存机制**：缓存中间计算结果避免重复计算

### 3. 可视化优化
- **矢量图形**：对大规模数据优先使用SVG格式
- **分层渲染**：复杂图件分层渲染提升性能
- **压缩输出**：使用压缩算法减少文件大小

## 📚 相关文档

- **[README.md](Scripts_New/README.md)**：详细的项目说明和快速开始指南
- **[DEVELOPMENT_PLAN.md](Scripts_New/DEVELOPMENT_PLAN.md)**：详细的开发计划和技术实现
- **[SCRIPTS_SUMMARY.md](Scripts_New/SCRIPTS_SUMMARY.md)**：所有脚本的详细功能说明
- **[Scripts_Backup/README.md](Scripts_Backup/README.md)**：原始脚本的说明和对比

## 🤝 贡献指南

### 代码贡献
1. 遵循现有的代码风格和命名规范
2. 为新功能添加相应的单元测试
3. 更新相关文档和配置文件
4. 提交前进行完整的测试验证

### 问题反馈
1. 提供详细的错误信息和复现步骤
2. 包含系统环境信息（R版本、包版本等）
3. 附上相关的配置文件（脱敏后）
4. 描述期望的行为和实际行为的差异

---

*本文档随项目更新而持续完善，最后更新时间：2025-10-27*