# Matrix2Markers 分析结果解读

基于对分析报告和原始数据的深入分析，以下是详细的生物学解读：

## 1. 质量控制 (QC) 总结
- **样本概况**：共分析了 4 个样本，总细胞数为 19,420 个。
- **数据质量**：
  - 所有样本的细胞数均充足（>3000 cells/sample），其中对照组 NCD_0W 细胞数最多（7519个）。
  - 双胞率（Doublet Rate）极低（<0.2%），说明数据纯净度高。
  - 过滤后的基因数量在 1.6万-2万之间，表明测序深度良好。
- **结论**：数据质量优异，无异常样本，适合进行下游比较分析。

## 2. 细胞类型注释 (Cell Type Annotation)
基于 Top Marker Genes，我们将 12 个聚类（Cluster 0-11）注释为以下细胞类型：

*   **Cluster 0 (Gzma, Ccl5, Nkg7, Prf1): NK Cells (Cytotoxic)**
    *   *依据*：高表达经典的 NK 细胞毒性颗粒基因（Gzma, Prf1）和趋化因子（Ccl5）。
*   **Cluster 1 (Cd7, Xcl1, Klrb1b, Cd160): NK Cells (Tissue-Resident/Memory-like)**
    *   *依据*：表达 NK 受体（Klrb1b/NK1.1）和组织驻留相关基因（Xcl1, Cd160）。
*   **Cluster 2 (Zbtb20, Ptprc, Itga4, Runx1): Lymphoid Progenitors / ILCs**
    *   *依据*：表达造血/淋巴系转录因子（Runx1, Zbtb20）和整合素（Itga4）。
*   **Cluster 3 (Rps/Rpl genes): Proliferating / Metabolically Active Cells**
    *   *依据*：主要特征是高表达核糖体蛋白基因，提示蛋白质合成活跃。
*   **Cluster 4 (Cdk8, Lars2, Malat1): Cycling / Stressed Cells**
    *   *依据*：表达细胞周期激酶（Cdk8）和长链非编码RNA（Malat1），可能处于应激或分裂状态。
*   **Cluster 5 (Tox, Cd226, Itga1): Exhausted / Activated NK Cells**
    *   *依据*：表达耗竭相关转录因子（Tox）和激活受体（Cd226）。
*   **Cluster 6 (Stmn1, Top2a, Mcm5, Mki67): Proliferating NK Cells**
    *   *依据*：高表达细胞周期标志物（Top2a, Mcm5），表明这是一群正在分裂的 NK 细胞。
*   **Cluster 7 (Cd74, Cd79a, Ms4a1): B Cells**
    *   *依据*：明确表达 B 细胞特征基因（Cd79a, Ms4a1/CD20）和 MHC-II 类分子（Cd74）。这可能混入的少量 B 细胞。
*   **Cluster 8 (Isg15, Stat1, Ifit1): Interferon-Stimulated NK Cells**
    *   *依据*：显著高表达干扰素刺激基因（ISGs），表明处于抗病毒或炎症激活状态。
*   **Cluster 9 (Cd3d, Cd3e, Trbc2): T Cells**
    *   *依据*：明确表达 T 细胞受体复合物基因（CD3D/E），这是 NK 细胞分选中常见的 T 细胞混杂。
*   **Cluster 10 (Ccl4, Ccl3, Egr1): Activated / Inflammatory NK Cells**
    *   *依据*：高表达炎症趋化因子（Ccl3, Ccl4）和早期激活基因（Egr1）。
*   **Cluster 11 (Lyz2, Cst3, Spi1): Myeloid Cells (Monocytes/Macrophages)**
    *   *依据*：表达典型的髓系标记（Lyz2, Spi1），为少量混杂细胞。

## 3. 细胞组成分析 (Composition)
- **主要群体**：
  - **Cluster 0 (Cytotoxic NK)** 是最主要的群体（6103个细胞），表明样本中绝大多数是具有杀伤功能的成熟 NK 细胞。
  - **Cluster 1 & 2** 也是主要群体，分别代表组织驻留型和前体样 NK 细胞。
- **混杂群体**：
  - 检测到少量的 **T 细胞 (Cluster 9)**、**B 细胞 (Cluster 7)** 和 **髓系细胞 (Cluster 11)**，这在流式分选后的 scRNA-seq 数据中是正常的背景噪音。
- **状态特异性群体**：
  - **Cluster 8 (IFN-response)** 和 **Cluster 6 (Proliferating)** 代表了 NK 细胞的特定功能状态（干扰素响应和增殖），值得比较其在 MCD 模型组与 NCD 对照组中的比例差异。

## 4. 生物学总结 (Executive Summary)
该数据集主要由 **NK 细胞** 组成，展现了从前体（Cluster 2）到成熟杀伤（Cluster 0）及组织驻留（Cluster 1）的完整分化谱系。我们识别出了具有特定功能状态的亚群，包括 **干扰素响应型 NK 细胞 (Cluster 8)** 和 **高炎症型 NK 细胞 (Cluster 10)**，这些亚群可能在 NASH 模型（MCD饮食）的进展中发挥关键作用。此外，数据中包含少量的 T/B/髓系细胞混杂，建议在后续针对 NK 细胞的深入分析中将其剔除。