"""
Script: Rename Clusters (02_rename_clusters.py)
功能：读取 h5ad 文件，根据配置文件对簇进行重命名（支持合并），并保存结果。
"""

import os
import sys
import argparse
import logging

# Add current directory to path to allow importing utils
# 将当前目录加入 path 以便导入 utils
sys.path.append(os.path.dirname(__file__))
import utils


def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="Cluster Renaming")
    parser.add_argument(
        "--config", required=True, help="Path to config file (配置文件路径)"
    )
    parser.add_argument(
        "--section", required=True, help="Config section (配置部分，如 'rename')"
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
    logger = setup_logger("logs/rename_clusters.log")

    logger.info("=== Starting Cluster Renaming (开始簇重命名) ===")

    # Extract params / 提取参数
    input_path = cfg["input_h5ad"]
    output_path = cfg["output_h5ad"]
    cluster_col = cfg["cluster_column"]
    new_col = cfg["new_column"]

    # Process mapping and merge rules / 处理映射和合并规则
    mapping = {}

    # 1. Load direct mapping if exists / 加载直接映射
    if "cluster_mapping" in cfg and cfg["cluster_mapping"]:
        mapping.update({str(k): str(v) for k, v in cfg["cluster_mapping"].items()})

    # 2. Load merge rules if exists / 加载合并规则
    # merge_rules format: {"NewName": ["old1", "old2"]}
    if "merge_rules" in cfg and cfg["merge_rules"]:
        merge_mapping = utils.parse_merge_rules(cfg["merge_rules"])
        # Merge rules override direct mapping if conflict
        mapping.update(merge_mapping)
        logger.info(f"Loaded merge rules for {len(merge_mapping)} clusters")

    logger.info(f"Input file (输入): {input_path}")
    logger.info(f"Output file (输出): {output_path}")
    logger.info(
        f"Mapping from '{cluster_col}' to '{new_col}' with {len(mapping)} rules"
    )

    # Read Data / 读取数据
    adata = utils.read_adata(input_path)
    logger.info(f"Loaded adata: {adata.n_obs} cells")

    # Rename / 执行重命名
    opts = cfg.get("options", {})
    utils.apply_cluster_mapping(
        adata,
        cluster_col,
        new_col,
        mapping,
        keep_original=opts.get("keep_original_column", True),
        error_if_unmapped=opts.get("error_if_unmapped", True),
        default_label=opts.get("default_label_for_unmapped", ""),
    )
    logger.info(f"Created new column (已创建新列): {new_col}")

    # Summarize / 生成汇总表
    df_summary = utils.summarize_mapping(adata, cluster_col, new_col)

    # TODO: Make this path configurable in YAML if needed
    summary_path = "results/tables/rename_mapping.tsv"
    utils.ensure_dir(summary_path)
    df_summary.to_csv(summary_path, sep="\t", index=False)
    logger.info(f"Saved mapping summary to {summary_path}")

    # Plot UMAP / 绘制 UMAP
    if cfg.get("plots", {}).get("do_umap", False):
        out_plot = cfg["plots"].get("umap_output", "results/plots/umap_renamed.png")
        color_by = cfg["plots"].get("umap_color_by", new_col)
        utils.plot_umap(adata, out_plot, color_by)
        logger.info(f"Saved UMAP plot to {out_plot}")

    # Save Data / 保存结果数据
    utils.save_adata(adata, output_path)
    logger.info(f"Saved renamed adata to {output_path}")
    logger.info("Done (完成).")


if __name__ == "__main__":
    main()
