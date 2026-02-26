"""
Shared utility functions for ClusterEditer workflow.
ClusterEditer 流程的共享工具函数。
"""

from __future__ import annotations

import os
import re
from typing import List, Tuple, Dict, Any, Optional

import yaml
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt


def load_config(config_path: str, section: Optional[str] = None) -> Dict[str, Any]:
    """
    Load YAML configuration file.
    加载 YAML 配置文件。

    参数:
        config_path (str): 配置文件路径
        section (str, optional): 指定读取的 section 名称，如 'prune' 或 'rename'

    返回:
        Dict[str, Any]: 配置字典
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    if section:
        if section not in cfg:
            raise KeyError(
                f"Section '{section}' not found in config file: {config_path}"
            )
        cfg = cfg[section]

    return cfg


def read_adata(input_h5ad: str) -> AnnData:
    """
    Read h5ad file and return AnnData object.
    读取 h5ad 文件并返回 AnnData 对象。
    """
    if not os.path.exists(input_h5ad):
        raise FileNotFoundError(f"h5ad file not found: {input_h5ad}")
    return sc.read_h5ad(input_h5ad)


def ensure_dir(path: str) -> None:
    """
    Ensure directory for a file exists.
    确保文件的父目录存在，如果不存在则创建。
    Handles both absolute and relative paths.
    """
    path = os.path.abspath(path)
    dirpath = os.path.dirname(path)
    if dirpath and not os.path.exists(dirpath):
        os.makedirs(dirpath, exist_ok=True)


def _remove_none_recursively(d: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively remove None values from nested dictionaries."""
    cleaned = {}
    for k, v in d.items():
        if v is None:
            print(f"[clean_adata] Removing None: uns['{k}']")
            continue
        elif isinstance(v, dict):
            cleaned[k] = _remove_none_recursively(v)
        else:
            cleaned[k] = v
    return cleaned


def clean_adata_for_export(adata: AnnData) -> None:
    """
    Clean AnnData object to ensure compatibility when saving.
    清理 AnnData 对象以确保保存时的兼容性。

    This function removes None values from uns that cause IORegistryError
    when reading with older anndata versions.
    此函数移除 uns 中的 None 值，避免旧版本 anndata 读取时出错。
    """
    # 1. Recursively clean uns: remove all None values at any depth
    adata.uns = _remove_none_recursively(dict(adata.uns))

    # 2. Clean obs: remove unused categories from categorical columns
    for col in adata.obs.columns:
        if hasattr(adata.obs[col], "cat"):
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()

    # 3. Clean var: remove unused categories from categorical columns
    for col in adata.var.columns:
        if hasattr(adata.var[col], "cat"):
            adata.var[col] = adata.var[col].cat.remove_unused_categories()

    # 4. Handle raw: remove it if it causes issues after subsetting
    # After subsetting, raw can have mismatched indices causing serialization issues
    if adata.raw is not None:
        try:
            for col in adata.raw.var.columns:
                if hasattr(adata.raw.var[col], "cat"):
                    adata.raw.var[col] = adata.raw.var[
                        col
                    ].cat.remove_unused_categories()
        except Exception:
            print("[clean_adata] Removing problematic raw attribute")
            adata.raw = None


def save_adata(
    adata: AnnData, output_path: str, compression: Optional[str] = None
) -> None:
    """
    Save AnnData to file, creating directory if needed.
    保存 AnnData 对象到文件，自动创建目录。

    Automatically cleans the AnnData to ensure cross-version compatibility.
    自动清理 AnnData 以确保跨版本兼容性。
    """
    ensure_dir(output_path)
    clean_adata_for_export(adata)
    adata.write(output_path, compression=compression)


def get_cluster_counts(adata: AnnData, cluster_col: str) -> pd.DataFrame:
    """
    Count cells per cluster.
    统计每个簇的细胞数量。

    参数:
        adata: AnnData 对象
        cluster_col: 聚类列名

    返回:
        DataFrame: 包含 'cluster' 和 'n_cells' 列，按细胞数降序排列
    """
    if cluster_col not in adata.obs.columns:
        raise KeyError(f"Column '{cluster_col}' not found in adata.obs")

    vc = adata.obs[cluster_col].value_counts()
    df_counts = pd.DataFrame({"cluster": vc.index.astype(str), "n_cells": vc.values})

    # Sort by n_cells descending
    df_counts = df_counts.sort_values("n_cells", ascending=False).reset_index(drop=True)
    return df_counts


def prune_clusters(
    adata: AnnData,
    cluster_col: str,
    clusters_to_remove: List[str],
) -> AnnData:
    """
    Remove cells belonging to specified clusters.
    删除属于指定簇的细胞（剪枝）。

    参数:
        adata: AnnData 对象
        cluster_col: 聚类列名
        clusters_to_remove: 需要删除的簇ID列表

    返回:
        AnnData: 剪枝后的新 AnnData 对象
    """
    if cluster_col not in adata.obs.columns:
        raise KeyError(f"Column '{cluster_col}' not found in adata.obs")

    cluster_series = adata.obs[cluster_col].astype(str)
    mask_keep = ~cluster_series.isin(clusters_to_remove)

    n_total = adata.n_obs
    n_keep = int(mask_keep.sum())
    n_removed = int(n_total - n_keep)

    print(f"[prune_clusters] total: {n_total}, kept: {n_keep}, removed: {n_removed}")

    return adata[mask_keep].copy()


def prepare_before_after_counts(
    df_before: pd.DataFrame,
    df_after: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compare cell counts before and after pruning.
    对比剪枝前后的细胞数量统计。
    """
    df_b = df_before.rename(columns={"n_cells": "n_cells_before"})
    df_a = df_after.rename(columns={"n_cells": "n_cells_after"})

    df_merged = pd.merge(df_b, df_a, on="cluster", how="outer")
    df_merged["n_cells_before"] = df_merged["n_cells_before"].fillna(0).astype(int)
    df_merged["n_cells_after"] = df_merged["n_cells_after"].fillna(0).astype(int)

    df_merged["n_cells_removed"] = (
        df_merged["n_cells_before"] - df_merged["n_cells_after"]
    )

    def _calc_pct_removed(row) -> float:
        before = row["n_cells_before"]
        if before <= 0:
            return 0.0
        return float(row["n_cells_removed"]) / float(before)

    df_merged["pct_removed"] = df_merged.apply(_calc_pct_removed, axis=1)
    df_merged = df_merged.sort_values("n_cells_before", ascending=False).reset_index(
        drop=True
    )

    return df_merged


def save_counts_tables(
    df_before: pd.DataFrame,
    df_after: pd.DataFrame,
    detail_path: str,
    summary_path: str,
) -> None:
    """
    Save detailed and summary count tables.
    保存详细和汇总的细胞计数表。
    """
    df_detail = prepare_before_after_counts(df_before, df_after)

    ensure_dir(detail_path)
    ensure_dir(summary_path)

    df_detail.to_csv(detail_path, sep="\t", index=False)

    total_before = int(df_detail["n_cells_before"].sum())
    total_after = int(df_detail["n_cells_after"].sum())
    total_removed = int(df_detail["n_cells_removed"].sum())
    pct_removed_total = (
        float(total_removed) / float(total_before) if total_before > 0 else 0.0
    )

    df_summary = pd.DataFrame(
        {
            "total_cells_before": [total_before],
            "total_cells_after": [total_after],
            "total_removed": [total_removed],
            "pct_removed_total": [pct_removed_total],
        }
    )

    df_summary.to_csv(summary_path, sep="\t", index=False)


def parse_merge_rules(merge_rules: Dict[str, List[Any]]) -> Dict[str, str]:
    """
    Parse merge rules into a one-to-one mapping dictionary.
    将合并规则解析为一对一的映射字典。

    参数:
        merge_rules: 字典，格式为 {"新名称": [旧ID1, 旧ID2, ...]}

    返回:
        Dict[str, str]: 格式为 {"旧ID1": "新名称", "旧ID2": "新名称"}
    """
    mapping = {}
    for new_name, old_ids in merge_rules.items():
        for old_id in old_ids:
            # Convert to string to ensure matching
            mapping[str(old_id)] = str(new_name)
    return mapping


def apply_cluster_mapping(
    adata: AnnData,
    cluster_col: str,
    new_col: str,
    mapping: Dict[str, str],
    keep_original: bool = True,
    error_if_unmapped: bool = True,
    default_label: str = "",
) -> None:
    """
    Apply cluster mapping to create new annotation.
    应用簇映射以创建新的注释（原地修改）。

    参数:
        adata: AnnData 对象
        cluster_col: 原始簇列名
        new_col: 新的注释列名
        mapping: 映射字典 {old_id: new_name}
        keep_original: 是否保留原列
        error_if_unmapped: 是否对未映射的簇报错
        default_label: 未映射簇的默认标签
    """
    if cluster_col not in adata.obs:
        raise KeyError(f"obs does not contain column '{cluster_col}'")

    original_values = adata.obs[cluster_col].astype(str)
    unique_clusters = sorted(original_values.unique())
    unmapped = [c for c in unique_clusters if c not in mapping]

    if unmapped:
        if error_if_unmapped:
            raise ValueError(f"Clusters not in mapping: {unmapped}")
        else:
            for c in unmapped:
                mapping[c] = default_label

    new_labels = original_values.map(mapping)

    if keep_original:
        adata.obs[new_col] = new_labels
    else:
        adata.obs[cluster_col] = new_labels


def summarize_mapping(adata: AnnData, cluster_col: str, new_col: str) -> pd.DataFrame:
    """
    Generate summary stats for mapping.
    生成映射后的统计汇总。
    """
    if cluster_col not in adata.obs:
        raise KeyError(f"'{cluster_col}' missing")
    if new_col not in adata.obs:
        raise KeyError(f"'{new_col}' missing")

    df_grouped = (
        adata.obs.groupby([cluster_col, new_col]).size().reset_index(name="n_cells")
    )
    total_cells = int(df_grouped["n_cells"].sum())
    df_grouped["pct_cells"] = df_grouped["n_cells"] / float(total_cells)

    return df_grouped.sort_values("n_cells", ascending=False).reset_index(drop=True)


def plot_umap(
    adata: AnnData,
    out_path: str,
    color: List[str],
    show: bool = False,
) -> None:
    """
    Plot UMAP with specified coloring.
    绘制 UMAP 并按指定列着色。
    """
    if "X_umap" not in adata.obsm_keys():
        print("No 'X_umap' found in adata.obsm, skipping plot.")
        return

    ensure_dir(out_path)

    # Use matplotlib backend explicitly to avoid issues
    sc.pl.umap(adata, color=color, show=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def refine_hvgs_by_gene_patterns(
    adata: AnnData,
    mt_pattern: str = "^MT-",
    rp_pattern: str = "^RP[SL]",
    ncrna_pattern: str = "^[A-Z][A-Z][0-9].*\\.[0-9]",
    linc_pattern: str = "(^LOC|LINC)[1-9]*",
) -> List[str]:
    """
    Refine HVGs by removing specific gene patterns (MT, RP, ncRNA, LINC, LOC).
    根据基因模式过滤高变基因（剔除线粒体、核糖体等干扰基因）。
    """
    if "highly_variable" not in adata.var.columns:
        print("[WARN] 'highly_variable' not found in adata.var. Skipping refinement.")
        return []

    hvgs = adata.var_names[adata.var["highly_variable"]].tolist()
    original_count = len(hvgs)

    patterns = [mt_pattern, rp_pattern, ncrna_pattern, linc_pattern]
    combined_pattern = "|".join([p for p in patterns if p])

    if not combined_pattern:
        return hvgs

    regex = re.compile(combined_pattern)
    filtered_hvgs = [gene for gene in hvgs if not regex.search(gene)]

    adata.var["highly_variable"] = False
    adata.var.loc[filtered_hvgs, "highly_variable"] = True

    print(f"[refine_hvgs] Filtered HVGs from {original_count} to {len(filtered_hvgs)}")
    return filtered_hvgs


def reprocess_adata(
    adata: AnnData,
    n_top_genes: int = 2000,
    target_sum: float = 1e4,
    scale_max_value: float = 10.0,
    n_neighbors: int = 10,
    gene_patterns: Dict[str, str] = {},
) -> AnnData:
    """
    Reprocess AnnData: Normalize -> HVG -> Filter HVG -> Scale -> PCA -> Neighbors -> UMAP.
    重新处理流程：归一化 -> 选高变 -> 过滤高变 -> Scale -> PCA -> Neighbors -> UMAP。

    Returns:
        AnnData: Processed AnnData with new embeddings.
    """
    # 1. Backup counts if not exists
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # 2. Normalize & Log1p
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    if "log1p" in adata.uns:
        del adata.uns["log1p"]

    # 3. HVG Selection
    print("Selecting HVGs...")
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=n_top_genes, layer="counts", subset=False
    )

    # 4. Refine HVG (Filter out MT/RP/etc)
    refine_hvgs_by_gene_patterns(
        adata,
        mt_pattern=gene_patterns.get("mt", ""),
        rp_pattern=gene_patterns.get("rp", ""),
        ncrna_pattern=gene_patterns.get("ncrna", ""),
        linc_pattern=gene_patterns.get("linc", ""),
    )

    # 5. Process on Subset (HVG) for PCA
    hv_mask = adata.var["highly_variable"].values
    adata_hvg = adata[:, hv_mask].copy()

    print("Scaling and PCA...")
    sc.pp.scale(adata_hvg, max_value=scale_max_value)
    sc.pp.pca(adata_hvg, svd_solver="arpack")

    print(f"Computing Neighbors (k={n_neighbors})...")
    sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors, use_rep="X_pca")

    print("Computing UMAP...")
    sc.tl.umap(adata_hvg)

    # Sync back to original adata
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]

    if "neighbors" in adata_hvg.uns:
        adata.uns["neighbors"] = adata_hvg.uns["neighbors"]

    return adata
