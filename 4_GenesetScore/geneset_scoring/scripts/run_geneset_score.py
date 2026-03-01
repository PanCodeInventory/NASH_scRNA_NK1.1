#!/usr/bin/env python
"""
run_geneset_score.py

功能：
- 基于基因集对单细胞进行打分（阳性/阴性分开）
- 计算 final_score = z(pos_score) - z(neg_score)
- 按 leiden 聚合得到簇级统计
- 输出簇注释候选（top1/top2 + delta + status）
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

sc = importlib.import_module("scanpy")


def log_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def ensure_dir(path: str) -> None:
    if not path:
        return
    os.makedirs(path, exist_ok=True)


def zscore_series(series: pd.Series) -> pd.Series:
    std = float(series.std(ddof=0))
    if std == 0.0 or np.isnan(std):
        return pd.Series(np.zeros(series.shape[0], dtype=float), index=series.index)
    mean = float(series.mean())
    return (series - mean) / std


def obs_series(adata: ad.AnnData, col: str) -> pd.Series:
    return pd.Series(adata.obs[col], index=adata.obs.index)


def unique_keep_order(items: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def safe_token(text: str) -> str:
    keep = []
    for ch in text.strip():
        if ch.isalnum() or ch in ["-", "_"]:
            keep.append(ch)
        else:
            keep.append("_")
    token = "".join(keep)
    while "__" in token:
        token = token.replace("__", "_")
    return token.strip("_") or "group"


def build_varname_lookup(var_names: List[str]) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for gene in var_names:
        lookup.setdefault(gene, gene)
        lookup.setdefault(gene.upper(), gene)
        lookup.setdefault(gene.title(), gene)
    return lookup


def resolve_genes(raw_genes: List[str], lookup: Dict[str, str]) -> Tuple[List[str], List[str]]:
    matched: List[str] = []
    missing: List[str] = []
    for g in raw_genes:
        key = g.strip()
        if not key:
            continue
        if key in lookup:
            matched.append(lookup[key])
            continue
        upper_key = key.upper()
        if upper_key in lookup:
            matched.append(lookup[upper_key])
            continue
        title_key = key.title()
        if title_key in lookup:
            matched.append(lookup[title_key])
            continue
        missing.append(g)
    return unique_keep_order(matched), missing


def load_genesets(path: str) -> List[Dict[str, List[str]]]:
    with open(path, "r", encoding="utf-8") as f:
        payload = yaml.safe_load(f)
    genesets = payload.get("genesets", [])
    if not isinstance(genesets, list) or not genesets:
        raise ValueError("No valid genesets found in YAML.")
    return genesets


def score_one_signature(
    adata: ad.AnnData,
    sig_name: str,
    pos_genes: List[str],
    neg_genes: List[str],
    gene_pool: List[str],
    random_state: int,
    use_raw: bool,
) -> Dict[str, str]:
    pos_raw_col = f"gs_{sig_name}_pos_raw"
    neg_raw_col = f"gs_{sig_name}_neg_raw"
    raw_col = f"gs_{sig_name}_raw"
    z_col = f"gs_{sig_name}_z"

    ctrl_size_pos = 50
    sc.tl.score_genes(
        adata,
        gene_list=pos_genes,
        score_name=pos_raw_col,
        ctrl_size=ctrl_size_pos,
        gene_pool=gene_pool,
        random_state=random_state,
        use_raw=use_raw,
    )

    if neg_genes:
        ctrl_size_neg = 50
        sc.tl.score_genes(
            adata,
            gene_list=neg_genes,
            score_name=neg_raw_col,
            ctrl_size=ctrl_size_neg,
            gene_pool=gene_pool,
            random_state=random_state,
            use_raw=use_raw,
        )
    else:
        adata.obs[neg_raw_col] = 0.0

    adata.obs[raw_col] = obs_series(adata, pos_raw_col).astype(float) - obs_series(adata, neg_raw_col).astype(float)
    adata.obs[z_col] = zscore_series(obs_series(adata, raw_col).astype(float))

    return {
        "pos_raw_col": pos_raw_col,
        "neg_raw_col": neg_raw_col,
        "raw_col": raw_col,
        "z_col": z_col,
    }


def parse_args() -> argparse.Namespace:
    module_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Gene set scoring for NK subtype characterization.")
    parser.add_argument(
        "--input-h5ad",
        type=str,
        default="/home/user/Pan Chongshi/Projects/NASH/NK1.1/Files/merged.processed.pruned.renamed.h5ad",
        help="Path to input h5ad",
    )
    parser.add_argument(
        "--geneset-yaml",
        type=str,
        default=str(module_root / "config" / "genesets.yaml"),
        help="Path to geneset YAML",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(module_root / "results"),
        help="Output directory",
    )
    parser.add_argument("--groupby", type=str, default="leiden", help="Cluster column in adata.obs")
    parser.add_argument(
        "--output-group-name",
        type=str,
        default="",
        help="Output column name for grouping key. Default: same as --groupby",
    )
    parser.add_argument("--ambiguous-delta", type=float, default=0.3, help="Top1-Top2 threshold for Ambiguous")
    parser.add_argument("--assign-min-z", type=float, default=0.7, help="Minimum top1 mean z for assignment")
    parser.add_argument("--assign-min-pct", type=float, default=0.25, help="Minimum top1 pct(z>=1) for assignment")
    parser.add_argument("--min-matched-genes", type=int, default=3, help="Minimum matched positive genes per signature")
    parser.add_argument("--random-state", type=int, default=0, help="Random state for scanpy score_genes")
    parser.add_argument(
        "--violin-dir",
        type=str,
        default="",
        help="Output directory for per-leiden violin plots. Default: <output-dir>/plots/leiden_violin",
    )
    return parser.parse_args()


def save_group_violin_plots(
    adata: ad.AnnData,
    groupby: str,
    signature_cols: Dict[str, Dict[str, str]],
    out_dir: str,
) -> None:
    ensure_dir(out_dir)
    grouped = adata.obs.groupby(groupby, observed=False)
    for cluster_id, cluster_obs in grouped:
        names: List[str] = []
        data: List[np.ndarray] = []
        for sig_name, cols in signature_cols.items():
            vals = cluster_obs[cols["z_col"]].astype(float).to_numpy()
            if vals.size == 0:
                continue
            names.append(sig_name)
            data.append(vals)
        if not data:
            continue

        fig_w = max(10.0, float(len(names)) * 1.25)
        fig, ax = plt.subplots(figsize=(fig_w, 5.5))
        violin = ax.violinplot(data, showmeans=True, showmedians=True)
        for body in violin["bodies"]:
            body.set_alpha(0.5)
        ax.set_xticks(np.arange(1, len(names) + 1))
        ax.set_xticklabels(names, rotation=30, ha="right")
        ax.set_ylabel("z score")
        ax.set_title(f"{groupby} = {cluster_id} signature score distribution")
        ax.axhline(0.0, linestyle="--", linewidth=1.0)
        fig.tight_layout()
        cluster_token = safe_token(str(cluster_id))
        group_token = safe_token(groupby)
        out_path = os.path.join(out_dir, f"{group_token}_{cluster_token}_violin.png")
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def save_cluster_signature_heatmap(cluster_df: pd.DataFrame, groupby: str, out_path: str) -> None:
    plot_df = cluster_df.pivot(index=groupby, columns="signature", values="mean_z")
    plot_df = plot_df.sort_index(axis=0)
    plot_df = plot_df.reindex(sorted(plot_df.columns), axis=1)

    fig_w = max(8.0, float(plot_df.shape[1]) * 1.1)
    fig_h = max(4.0, float(plot_df.shape[0]) * 0.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = float(np.nanmax(np.abs(plot_df.to_numpy())))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0
    im = ax.imshow(plot_df.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(plot_df.shape[1]))
    ax.set_xticklabels(plot_df.columns.tolist(), rotation=30, ha="right")
    ax.set_yticks(np.arange(plot_df.shape[0]))
    ax.set_yticklabels(plot_df.index.astype(str).tolist())
    ax.set_xlabel("signature")
    ax.set_ylabel(groupby)
    ax.set_title("Cluster signature mean z heatmap")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("mean z")
    fig.tight_layout()

    ensure_dir(os.path.dirname(out_path))
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_group_name = args.output_group_name.strip() if args.output_group_name else args.groupby
    group_token = safe_token(output_group_name)

    ensure_dir(args.output_dir)

    try:
        print(f"[{log_timestamp()}] [INFO] Loading AnnData: {args.input_h5ad}")
        adata = ad.read_h5ad(args.input_h5ad)
        print(f"[{log_timestamp()}] [INFO] Shape: cells={adata.n_obs}, genes={adata.n_vars}")

        if args.groupby not in adata.obs.columns:
            raise KeyError(f"groupby column '{args.groupby}' not found in adata.obs")

        print(f"[{log_timestamp()}] [INFO] Loading genesets: {args.geneset_yaml}")
        genesets = load_genesets(args.geneset_yaml)

        use_raw = adata.raw is not None
        if use_raw:
            score_var_names = list(adata.raw.var_names)
        else:
            score_var_names = list(adata.var_names)

        var_lookup = build_varname_lookup(score_var_names)
        gene_pool = list(score_var_names)

        print(f"[{log_timestamp()}] [INFO] score_genes source: {'adata.raw' if use_raw else 'adata.X'}")

        report_rows: List[Dict[str, object]] = []
        signature_cols: Dict[str, Dict[str, str]] = {}
        resolved_meta: List[Dict[str, Any]] = []

        for gs in genesets:
            sig_name = str(gs.get("name", "")).strip()
            if not sig_name:
                continue
            sig_key = sig_name.replace(" ", "_").replace("-", "_")

            pos_input = [str(x) for x in gs.get("positive", [])]
            neg_input = [str(x) for x in gs.get("negative", [])]

            pos_genes, pos_missing = resolve_genes(pos_input, var_lookup)
            neg_genes, neg_missing = resolve_genes(neg_input, var_lookup)

            pos_match_rate = (len(pos_genes) / len(pos_input)) if pos_input else 0.0
            neg_match_rate = (len(neg_genes) / len(neg_input)) if neg_input else 1.0

            skipped = False
            skip_reason = ""
            if len(pos_genes) < args.min_matched_genes:
                skipped = True
                skip_reason = f"matched positive genes {len(pos_genes)} < min {args.min_matched_genes}"

            if not skipped:
                cols = score_one_signature(
                    adata=adata,
                    sig_name=sig_key,
                    pos_genes=pos_genes,
                    neg_genes=neg_genes,
                    gene_pool=gene_pool,
                    random_state=args.random_state,
                    use_raw=use_raw,
                )
                signature_cols[sig_name] = cols
                print(
                    f"[{log_timestamp()}] [INFO] Scored signature '{sig_name}': "
                    f"pos={len(pos_genes)}, neg={len(neg_genes)}"
                )
            else:
                print(f"[{log_timestamp()}] [WARN] Skip signature '{sig_name}': {skip_reason}")

            report_rows.append(
                {
                    "signature": sig_name,
                    "n_pos_input": len(pos_input),
                    "n_pos_matched": len(pos_genes),
                    "n_pos_missing": len(pos_missing),
                    "pos_missing_genes": ",".join(pos_missing),
                    "pos_match_rate": round(pos_match_rate, 4),
                    "n_neg_input": len(neg_input),
                    "n_neg_matched": len(neg_genes),
                    "n_neg_missing": len(neg_missing),
                    "neg_missing_genes": ",".join(neg_missing),
                    "neg_match_rate": round(neg_match_rate, 4),
                    "skipped": skipped,
                    "skip_reason": skip_reason,
                }
            )
            resolved_meta.append(
                {
                    "signature": sig_name,
                    "positive_resolved": pos_genes,
                    "negative_resolved": neg_genes,
                    "positive_missing": pos_missing,
                    "negative_missing": neg_missing,
                    "skipped": skipped,
                    "skip_reason": skip_reason,
                }
            )

        report_df = pd.DataFrame(report_rows)
        report_path = os.path.join(args.output_dir, "geneset_match_report.tsv")
        report_df.to_csv(report_path, sep="\t", index=False)
        print(f"[{log_timestamp()}] [INFO] Saved: {report_path}")

        if not signature_cols:
            raise RuntimeError("No signatures scored successfully. Check geneset matching report.")

        cell_base_cols = [args.groupby]
        for extra in ["sample", "group", "time", "cluster_annot"]:
            if extra in adata.obs.columns:
                cell_base_cols.append(extra)
        cell_df = adata.obs[cell_base_cols].copy()
        if output_group_name != args.groupby:
            cell_df[output_group_name] = cell_df[args.groupby].astype(str)
            cell_df = cell_df.drop(columns=[args.groupby])
        cell_df.insert(0, "cell_id", adata.obs.index.astype(str))

        for sig_name, cols in signature_cols.items():
            cell_df[f"{sig_name}__pos_raw"] = obs_series(adata, cols["pos_raw_col"]).astype(float)
            cell_df[f"{sig_name}__neg_raw"] = obs_series(adata, cols["neg_raw_col"]).astype(float)
            cell_df[f"{sig_name}__raw"] = obs_series(adata, cols["raw_col"]).astype(float)
            cell_df[f"{sig_name}__z"] = obs_series(adata, cols["z_col"]).astype(float)

        cell_scores_path = os.path.join(args.output_dir, "cell_scores.tsv.gz")
        cell_df.to_csv(cell_scores_path, sep="\t", index=False, compression="gzip")
        print(f"[{log_timestamp()}] [INFO] Saved: {cell_scores_path}")

        cluster_rows: List[Dict[str, object]] = []
        grouped = adata.obs.groupby(args.groupby, observed=False)
        for cluster_id, cluster_obs in grouped:
            for sig_name, cols in signature_cols.items():
                raw_vals = cluster_obs[cols["raw_col"]].astype(float)
                z_vals = cluster_obs[cols["z_col"]].astype(float)
                cluster_rows.append(
                    {
                        output_group_name: str(cluster_id),
                        "signature": sig_name,
                        "mean_raw": float(raw_vals.mean()),
                        "mean_z": float(z_vals.mean()),
                        "median_z": float(z_vals.median()),
                        "std_z": float(z_vals.std(ddof=0)),
                        "pct_z_ge_1": float((z_vals >= 1.0).mean()),
                        "n_cells": int(z_vals.shape[0]),
                    }
                )

        cluster_df = pd.DataFrame(cluster_rows)
        cluster_scores_path = os.path.join(args.output_dir, f"cluster_scores_by_{group_token}.tsv")
        cluster_df.to_csv(cluster_scores_path, sep="\t", index=False)
        print(f"[{log_timestamp()}] [INFO] Saved: {cluster_scores_path}")

        candidate_rows: List[Dict[str, object]] = []
        for cluster_id, one_cluster_df in cluster_df.groupby(output_group_name):
            ranked = one_cluster_df.sort_values("mean_z", ascending=False).reset_index(drop=True)
            for idx in range(ranked.shape[0]):
                row = ranked.iloc[idx]
                next_mean_z = float(ranked.iloc[idx + 1]["mean_z"]) if idx + 1 < ranked.shape[0] else float("nan")
                delta_to_next = float(row["mean_z"] - next_mean_z) if idx + 1 < ranked.shape[0] else float("inf")
                candidate_rows.append(
                    {
                        output_group_name: str(cluster_id),
                        "rank": int(idx + 1),
                        "signature": str(row["signature"]),
                        "mean_z": float(row["mean_z"]),
                        "pct_z_ge_1": float(row["pct_z_ge_1"]),
                        "delta_to_next": delta_to_next,
                        "n_cells": int(row["n_cells"]),
                    }
                )

        candidates_df = pd.DataFrame(candidate_rows)

        recommended_labels: Dict[str, str] = {}
        for cluster_id, one_cluster_df in candidates_df.groupby(output_group_name):
            first = one_cluster_df.sort_values("rank", ascending=True).iloc[0]
            top1_sig = str(first["signature"])
            top1_mean_z = float(first["mean_z"])
            top1_pct = float(first["pct_z_ge_1"])
            delta = float(first["delta_to_next"])
            if top1_mean_z < args.assign_min_z or top1_pct < args.assign_min_pct:
                recommended = "unassigned"
            elif delta < args.ambiguous_delta:
                recommended = "ambiguous"
            else:
                recommended = top1_sig
            recommended_labels[str(cluster_id)] = recommended

        candidates_df["recommended_label"] = candidates_df[output_group_name].astype(str).apply(
            lambda x: recommended_labels.get(x, "unassigned")
        )
        candidates_path = os.path.join(args.output_dir, f"cluster_annotation_candidates_by_{group_token}.tsv")
        candidates_df.to_csv(candidates_path, sep="\t", index=False)
        print(f"[{log_timestamp()}] [INFO] Saved: {candidates_path}")

        violin_dir = args.violin_dir if args.violin_dir else os.path.join(args.output_dir, "plots", "leiden_violin")
        save_group_violin_plots(
            adata=adata,
            groupby=args.groupby,
            signature_cols=signature_cols,
            out_dir=violin_dir,
        )
        print(f"[{log_timestamp()}] [INFO] Saved violin plots: {violin_dir}")

        heatmap_path = os.path.join(args.output_dir, "plots", f"cluster_signature_mean_z_heatmap_by_{group_token}.png")
        save_cluster_signature_heatmap(cluster_df=cluster_df, groupby=output_group_name, out_path=heatmap_path)
        print(f"[{log_timestamp()}] [INFO] Saved heatmap: {heatmap_path}")

        adata.uns["genesets_scoring"] = {
            "input_h5ad": args.input_h5ad,
            "geneset_yaml": args.geneset_yaml,
            "groupby": args.groupby,
            "output_group_name": output_group_name,
            "use_raw": use_raw,
            "ctrl_size": 50,
            "n_bins": 25,
            "gene_pool_rule": "all_genes_in_scoring_matrix",
            "z_scope": "global_per_signature",
            "random_state": args.random_state,
            "ambiguous_delta": args.ambiguous_delta,
            "assign_min_z": args.assign_min_z,
            "assign_min_pct": args.assign_min_pct,
            "min_matched_genes": args.min_matched_genes,
            "resolved_genesets_json": json.dumps(resolved_meta, ensure_ascii=False),
        }

        scored_h5ad_path = os.path.join(args.output_dir, "adata.geneset_scored.h5ad")
        adata.write(scored_h5ad_path)
        print(f"[{log_timestamp()}] [INFO] Saved: {scored_h5ad_path}")

        print(f"[{log_timestamp()}] [INFO] Done.")

    except Exception as exc:
        print(f"[ERROR] Gene set scoring failed: {exc}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
