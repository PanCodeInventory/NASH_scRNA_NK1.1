"""
Script: Prune Clusters (01_prune_clusters.py)
功能：读取 h5ad 文件，根据配置文件删除指定的簇（Cluster），并保存处理后的数据。
"""

import os
import sys
import argparse
import logging
from datetime import datetime

# Add current directory to path to allow importing utils
# 将当前目录加入 path 以便导入 utils
sys.path.append(os.path.dirname(__file__))

import utils


def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="Cluster Pruning")
    parser.add_argument(
        "--config", required=True, help="Path to config file (配置文件路径)"
    )
    parser.add_argument(
        "--section", required=True, help="Config section (配置部分，如 'prune')"
    )
    return parser.parse_args()


def setup_logger(log_path: str):
    """配置日志记录"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(__name__)


def main():
    """主流程函数"""
    args = parse_args()

    # Load config / 加载配置
    cfg = utils.load_config(args.config, args.section)

    # Setup paths / 设置路径
    os.makedirs("logs", exist_ok=True)
    logger = setup_logger("logs/prune_clusters.log")

    logger.info("=== Starting Cluster Pruning (开始簇剪枝) ===")

    # Extract params / 提取参数
    input_path = cfg["input_h5ad"]
    output_path = cfg["output"]["pruned_h5ad"]
    cluster_col = cfg["cluster_column"]
    clusters_to_remove = [str(c) for c in cfg["clusters_to_remove"]]

    logger.info(f"Input file (输入): {input_path}")
    logger.info(f"Output file (输出): {output_path}")
    logger.info(f"Removing clusters from column '{cluster_col}': {clusters_to_remove}")

    # Read Data / 读取数据
    adata = utils.read_adata(input_path)
    logger.info(f"Loaded adata: {adata.n_obs} cells")

    # Stats before / 剪枝前统计
    df_before = utils.get_cluster_counts(adata, cluster_col)

    # Prune / 执行剪枝
    adata_pruned = utils.prune_clusters(adata, cluster_col, clusters_to_remove)
    logger.info(f"Pruned adata: {adata_pruned.n_obs} cells (after pruning)")

    # Stats after / 剪枝后统计
    df_after = utils.get_cluster_counts(adata_pruned, cluster_col)

    # Save stats / 保存统计表
    if cfg.get("options", {}).get("save_cell_counts_table", True):
        utils.save_counts_tables(
            df_before,
            df_after,
            "results/tables/prune_detail.tsv",
            "results/tables/prune_summary.tsv",
        )
        logger.info("Saved statistics tables (统计表已保存)")

    # Plot UMAP / 绘制 UMAP
    if cfg.get("plots", {}).get("do_umap_after_prune", False):
        umap_cols = cfg["plots"].get("umap_color_by", [cluster_col])
        utils.plot_umap(adata_pruned, "results/plots/umap_after_prune.png", umap_cols)
        logger.info("Saved UMAP plot (UMAP图已保存)")

    # Save Data / 保存结果数据
    utils.save_adata(adata_pruned, output_path)
    logger.info("Done (完成).")


if __name__ == "__main__":
    main()
