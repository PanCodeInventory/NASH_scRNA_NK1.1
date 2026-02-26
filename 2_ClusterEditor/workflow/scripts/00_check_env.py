"""
Script: Check Environment (00_check_env.py)
功能：检查当前运行环境是否满足项目依赖要求。
"""

import sys
import pkg_resources
import logging


def setup_logger():
    """配置日志输出到控制台"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(__name__)


def check_requirements(requirements_path: str):
    """
    检查 requirements.txt 中的依赖是否已安装。

    参数:
        requirements_path: requirements.txt 文件的路径
    """
    logger = setup_logger()
    logger.info("=== 开始环境检查 ===")

    try:
        with open(requirements_path, "r") as f:
            requirements = [
                str(requirement) for requirement in pkg_resources.parse_requirements(f)
            ]
    except FileNotFoundError:
        logger.error(f"找不到文件: {requirements_path}")
        sys.exit(1)

    missing = []
    for requirement in requirements:
        try:
            pkg_resources.require(requirement)
            logger.info(f"✓ {requirement} 已安装")
        except pkg_resources.DistributionNotFound:
            missing.append(requirement)
            logger.error(f"✗ {requirement} 未安装")
        except pkg_resources.VersionConflict as e:
            missing.append(requirement)
            logger.error(f"✗ {requirement} 版本冲突: {e}")

    if missing:
        logger.error(f"环境检查失败。缺少的包: {', '.join(missing)}")
        sys.exit(1)

    logger.info("=== 环境检查通过 ===")


if __name__ == "__main__":
    # 默认检查根目录下的 requirements.txt
    check_requirements("requirements.txt")
