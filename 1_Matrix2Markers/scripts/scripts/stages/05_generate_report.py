import argparse
import os
import sys
import yaml
import pandas as pd
import base64
from datetime import datetime

# 获取项目根目录 (假设当前脚本在 stages/ 目录下，上一级是根目录)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PROJECT_ROOT)

def load_config():
    config_path = os.path.join(PROJECT_ROOT, "config.yaml")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def img_to_base64(img_path):
    """将图片转换为 Base64 字符串以嵌入 HTML"""
    if not os.path.exists(img_path):
        return ""
    with open(img_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')

def generate_html_report(config, output_file):
    dirs = config['directories']
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # 1. 准备数据路径
    # 使用 config['directories'] 中的路径，假设所有 artifacts 都在那里
    qc_summary_path = os.path.join(dirs['tables'], "qc_single_summary.tsv")
    cell_counts_path = os.path.join(dirs['tables'], "cell_counts_by_leiden.tsv")
    markers_path = os.path.join(dirs['tables'], "all_markers.csv")
    umap_path = os.path.join(dirs['plots'], "umap", "umap_leiden.png")
    dotplot_path = os.path.join(dirs['plots'], "markers", "dotplot_leiden_top5.png")

    # 2. 读取表格数据
    qc_html = "<p>No QC data found.</p>"
    if os.path.exists(qc_summary_path):
        df = pd.read_csv(qc_summary_path, sep='\t')
        qc_html = df.to_html(classes="table table-striped", index=False, border=0)

    counts_html = "<p>No cell count data found.</p>"
    if os.path.exists(cell_counts_path):
        df = pd.read_csv(cell_counts_path, sep='\t', index_col=0)
        counts_html = df.to_html(classes="table table-striped", border=0)

    markers_html = "<p>No marker data found.</p>"
    top_markers_summary = ""
    if os.path.exists(markers_path):
        df = pd.read_csv(markers_path)
        # 简化 Marker 表：只取每个 Cluster 的前 5 个
        if 'cluster' in df.columns and 'names' in df.columns:
            # 假设 df 已经是排好序的，或者我们需要根据 score/pvals 排序
            # 这里简单取头部
            summary_data = []
            for cl in sorted(df['cluster'].unique()):
                top_genes = df[df['cluster'] == cl].head(5)['names'].tolist()
                summary_data.append({"Cluster": cl, "Top Markers": ", ".join(top_genes)})
            
            markers_summary_df = pd.DataFrame(summary_data)
            markers_html = markers_summary_df.to_html(classes="table table-striped", index=False, border=0)
            
            # 为 Agent 准备的纯文本摘要
            top_markers_summary = markers_summary_df.to_string(index=False)

    # 3. 准备图片 (Base64 嵌入)
    umap_b64 = img_to_base64(umap_path)
    dotplot_b64 = img_to_base64(dotplot_path)

    # --- 准备分板块的 AI 上下文数据 ---

    # A. QC Context
    qc_context_txt = f"Total Samples: {len(config.get('samples', {}))}. "
    if os.path.exists(qc_summary_path):
         qcdf = pd.read_csv(qc_summary_path, sep='\t')
         if 'n_cells_after_qc' in qcdf.columns:
            low_cell_samples = qcdf[qcdf['n_cells_after_qc'] < 500]['sample'].tolist()
            if low_cell_samples:
                qc_context_txt += f"WARNING: Samples with <500 cells: {', '.join(low_cell_samples)}. "
            else:
                qc_context_txt += "All samples have >500 cells. "
    
    qc_section = f"""
        <div class="section">
            <h2>1. Quality Control Summary</h2>
            <div class="table-responsive">
                {qc_html}
            </div>
            <div class="hidden-context" data-section="qc">
                [AI_DATA_QC]
                {qc_context_txt}
                Instructions: Check for samples with low cell counts or abnormal metrics.
            </div>
        </div>
    """

    # B. UMAP Context (No AI Context as requested)
    umap_section = f"""
        <div class="section">
            <h2>2. Visualization (UMAP)</h2>
            <div class="img-container">
                {'<img src="data:image/png;base64,' + umap_b64 + '" alt="UMAP Plot">' if umap_b64 else '<p>UMAP plot not found</p>'}
            </div>
        </div>
    """

    # C. Composition Context
    comp_context_txt = ""
    if os.path.exists(cell_counts_path):
        cdf = pd.read_csv(cell_counts_path, sep='\t', index_col=0)
        top_clusters = cdf.sum().sort_values(ascending=False).head(3).index.tolist()
        comp_context_txt += f"Dominant Clusters: {', '.join(map(str, top_clusters))}. "
        
        col_norm = cdf / cdf.sum(axis=0)
        specific_clusters = []
        for cl in col_norm.columns:
            max_sample = col_norm[cl].idxmax()
            max_prop = col_norm[cl].max()
            if max_prop > 0.8:
                specific_clusters.append(f"{cl} (dominated by {max_sample})")
        
        if specific_clusters:
            comp_context_txt += f"Sample-Specific Clusters: {', '.join(specific_clusters)}."
    
    comp_section = f"""
        <div class="section">
            <h2>4. Cell Composition</h2>
            <div class="table-responsive">
                {counts_html}
            </div>
            <div class="hidden-context" data-section="composition">
                [AI_DATA_COMPOSITION]
                {comp_context_txt}
                Instructions: Analyze the abundance of cell types. Note any clusters that are specific to certain samples.
                CRITICAL: Use the cell types you identified in the 'Markers' section to explain these composition changes (e.g., 'T cells are expanded in Sample A' instead of 'Cluster 1').
            </div>
        </div>
    """

    # D. Markers Context
    marker_section = f"""
        <div class="section">
            <h2>3. Marker Genes & Annotation</h2>
            <div class="img-container">
                {'<img src="data:image/png;base64,' + dotplot_b64 + '" alt="Dot Plot">' if dotplot_b64 else '<p>Dot plot not found</p>'}
            </div>
            <h3>Top Markers per Cluster</h3>
            {markers_html}
            <div class="hidden-context" data-section="markers">
                [AI_DATA_MARKERS]
                Here is the list of top marker genes for annotation:
                {top_markers_summary}
                Instructions: Propose cell types for each cluster based on these genes.
            </div>
        </div>
    """

    # 4. 构建 HTML 内容 (顺序: QC -> UMAP -> Markers -> Composition)
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Matrix2Markers Analysis Report</title>
        <style>
            body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; line-height: 1.6; color: #333; max-width: 1200px; margin: 0 auto; padding: 20px; }}
            h1, h2, h3 {{ color: #2c3e50; }}
            .section {{ margin-bottom: 40px; background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            .table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
            .table th, .table td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
            .table th {{ background-color: #f8f9fa; }}
            .img-container {{ text-align: center; margin: 20px 0; }}
            img {{ max-width: 100%; height: auto; border-radius: 4px; box-shadow: 0 2px 8px rgba(0,0,0,0.15); }}
            .meta {{ color: #666; font-size: 0.9em; margin-bottom: 30px; }}
            .hidden-context {{ display: none; }}
        </style>
    </head>
    <body>
        <h1>Matrix2Markers Analysis Report</h1>
        <div class="meta">
            Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} <br>
            Samples: {len(config.get('samples', {}))}
        </div>

        {qc_section}
        {umap_section}
        {marker_section}
        {comp_section}

    </body>
    </html>
    """

    with open(output_file, "w") as f:
        f.write(html_content)
    
    print(f"HTML Report generated: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate HTML Analysis Report")
    parser.add_argument("--output-file", required=True, help="Path to the output HTML report file")
    args = parser.parse_args()

    cfg = load_config()
    generate_html_report(cfg, args.output_file)