#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Xcl1 violin across cluster_annot groups."
    )
    parser.add_argument(
        "--input-h5ad",
        type=str,
        default="/home/user/Pan Chongshi/Projects/NASH/NK1.1/Files/merged.processed.pruned.renamed.h5ad",
        help="Path to input h5ad file.",
    )
    parser.add_argument(
        "--gene",
        type=str,
        default="Xcl1",
        help="Gene symbol to plot.",
    )
    parser.add_argument(
        "--groupby",
        type=str,
        default="cluster_annot",
        help="Column in adata.obs for grouping. Prefer cluster_annot.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "results"),
        help="Output directory.",
    )
    parser.add_argument(
        "--use-raw",
        action="store_true",
        help="Plot from adata.raw if available.",
    )
    return parser.parse_args()


def resolve_groupby(adata: Any, groupby: str) -> str:
    if groupby in adata.obs.columns:
        return groupby
    if groupby == "cluster_anno" and "cluster_annot" in adata.obs.columns:
        return "cluster_annot"
    raise KeyError(f"Grouping column '{groupby}' not found in adata.obs")


def resolve_gene(adata: Any, gene: str, use_raw: bool) -> str:
    var_names = adata.raw.var_names if use_raw and adata.raw is not None else adata.var_names
    if gene in var_names:
        return gene
    lookup = {name.upper(): name for name in var_names}
    upper_gene = gene.upper()
    if upper_gene in lookup:
        return lookup[upper_gene]
    raise KeyError(f"Gene '{gene}' not found in expression matrix")


def main() -> None:
    args = parse_args()

    try:
        ad = importlib.import_module("anndata")
        sc = importlib.import_module("scanpy")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Missing dependency. Please install scanpy/anndata in the current environment before running this script."
        ) from exc

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    groupby = resolve_groupby(adata, args.groupby)
    gene = resolve_gene(adata, args.gene, args.use_raw)
    use_raw = bool(args.use_raw and adata.raw is not None)

    cluster_counts = (
        adata.obs[groupby]
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis(groupby)
        .reset_index(name="n_cells")
    )
    summary_path = output_dir / f"{groupby}_summary.tsv"
    cluster_counts.to_csv(summary_path, sep="\t", index=False)

    order = cluster_counts[groupby].tolist()

    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.violin(
        adata,
        keys=gene,
        groupby=groupby,
        order=order,
        stripplot=False,
        rotation=30,
        use_raw=use_raw,
        ax=ax,
        show=False,
    )
    ax.set_title(f"{gene} expression by {groupby}")
    ax.set_xlabel(groupby)
    ax.set_ylabel("expression")
    fig.tight_layout()

    png_path = output_dir / f"{gene}_{groupby}_violin.png"
    pdf_path = output_dir / f"{gene}_{groupby}_violin.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    print(f"groupby used: {groupby}")
    print(f"gene used: {gene}")
    print(f"saved: {summary_path}")
    print(f"saved: {png_path}")
    print(f"saved: {pdf_path}")


if __name__ == "__main__":
    main()
