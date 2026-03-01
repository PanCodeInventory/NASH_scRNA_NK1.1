背景：

- 参考数据：CIMA PBMC_markers.xlsx（分层细胞类型注释系统）
- 目标：设计基因集打分系统，对每个簇进行 NK 亚型注释

**CIMA NK 细胞层级结构（L2-L4）：**

| L2 亚型    | L3 亚型     | L4 亚型                 | 阳性 Marker                         | 阴性 Marker |
| ---------- | ----------- | ----------------------- | ----------------------------------- | ----------- |
| ILC2       | ILC2        | ILC2-IL2RA              | IL2RA, KIT, ILIR1                   | -           |
| NKbright   | NKbright    | NKbright-XCL1           | XCL1, IGFBP4, GZMK, SELL, TNFRSF11A | -           |
| Transit.NK | Transit.NK  | Transitional NK-GZMK    | GZMK, CXCR4, ITGA6, PIK3R1, PDE4B   | -           |
| NKdim      | Mature NK   | Mature NKdim-FCGR3A     | FCGR3A, CX3CR1, FGFBP2, GZMB, SPON2 | GZMK        |
| NKdim      | Mature NK   | Inflamed NKdim-IFIT1    | IFIT1, OAS3, IFI44L                 | -           |
| NKdim      | Terminal NK | Terminal NKdim-CD160neg | LILRB1, ZEB2, ZNF700, PCSK7, CLEC2D | CD160       |
| Cycling NK | Cycling NK  | Cycling NK              | MKI67, TYMS, STMN1                  | -           |

**混合阴阳 Marker 处理策略：**

- 方案：分离打分 + 标准化差值法
- 阳性marker和阴性marker分别打分
- 最终得分 = Z(阳性分) - Z(阴性分)

**工作流程规划：**

1. 解析 Marker 文件 → 拆分阴阳marker
2. 构建基因集字典
3. 对每个簇/细胞进行打分
4. 注释分配（取最高分亚型）
5. 可视化验证