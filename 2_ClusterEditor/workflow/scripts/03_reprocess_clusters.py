"""
Script: Reprocess Clusters (03_reprocess_clusters.py)
功能：在修剪和重命名后，基于纯化后的高变基因（去除线粒体/核糖体等）重新进行降维（PCA）、邻管图构建、UMAP 和聚类。
这样可以去除已删除簇的影响，使 UMAP 可视化更加干净、美观。
"""

import os
import sys
import argparse
import logging

# Add current directory to path to allow importing utils
sys.path.append(os.path.dirname(__file__))
import utils


def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="Cluster Reprocessing")
    parser.add_argument("--config", required=True, help="Path to config file")
    parser.add_argument(
        "--section", required=True, help="Config section (e.g. 'reprocess')"
    )
    return parser.parse_args()


def setup_logger(log_path: str):
    """配置日志"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(__name__)


def main():
    args = parse_args()

    # Load Config
    cfg = utils.load_config(args.config, args.section)

    # Setup Paths
    os.makedirs("logs", exist_ok=True)
    logger = setup_logger("logs/reprocess_clusters.log")

    logger.info("=== Starting Cluster Reprocessing (开始重新降维聚类) ===")

    input_path = cfg["input_h5ad"]
    output_path = cfg["output_h5ad"]

    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_path}")

    # Read Data
    adata = utils.read_adata(input_path)
    logger.info(f"Loaded adata: {adata.n_obs} cells, {adata.n_vars} genes")

    # Reprocess
    logger.info(
        "Running reprocessing pipeline (Normalize -> refined HVG -> Scale -> PCA -> UMAP)..."
    )

    # Extract params with defaults
    n_top_genes = cfg.get("n_top_genes", 2000)
    target_sum = cfg.get("target_sum", 1e4)
    scale_max = cfg.get("scale_max_value", 10.0)
    n_neighbors = cfg.get("n_neighbors", 10)
    gene_patterns = cfg.get("gene_patterns", {})

    adata_processed = utils.reprocess_adata(
        adata,
        n_top_genes=n_top_genes,
        target_sum=target_sum,
        scale_max_value=scale_max,
        n_neighbors=n_neighbors,
        gene_patterns=gene_patterns,
    )

    logger.info("Reprocessing complete.")

    # Plotting
    if cfg.get("plots", {}).get("do_umap", False):
        plot_cfg = cfg["plots"]
        out_dir = plot_cfg.get("umap_output_dir", "results/plots/reprocess")
        color_by = plot_cfg.get("umap_color_by", ["cluster_annot"])

        utils.ensure_dir(os.path.join(out_dir, "placeholder.txt"))

        for col in color_by:
            if col in adata_processed.obs.columns:
                out_fig = os.path.join(out_dir, f"umap_reprocessed_by_{col}.png")
                logger.info(f"Saving plot to {out_fig}")
                utils.plot_umap(adata_processed, out_fig, [col])
            else:
                logger.warning(f"Column '{col}' not found in obs, skipping plot.")

    # Save
    utils.save_adata(adata_processed, output_path)
    logger.info(f"Saved processed adata to {output_path}")
    logger.info("Done.")


if __name__ == "__main__":
    main()
