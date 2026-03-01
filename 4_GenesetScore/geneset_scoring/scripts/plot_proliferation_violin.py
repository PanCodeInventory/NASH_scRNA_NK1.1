#!/usr/bin/env python
"""
plot_proliferation_violin.py

Generate violin plots comparing proliferation scores across clusters.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    # Paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.dirname(script_dir) + "/results"
    cell_scores_path = results_dir + "/cell_scores.tsv.gz"
    output_dir = results_dir + "/plots/proliferation"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load cell-level scores
    print(f"Loading: {cell_scores_path}")
    df = pd.read_csv(cell_scores_path, sep="\t", compression="gzip")
    
    # Proliferation signatures
    sig_short = "Cycling NK"
    sig_extended = "Proliferation-Extended"
    
    z_col_short = f"{sig_short}__z"
    z_col_extended = f"{sig_extended}__z"
    
    # Cluster order
    cluster_order = ["Cytotoxic NK", "ILC1", "Circulating NK", "High-Ribosomal NK", "Homeostatic NK"]
    
    # =========== Plot 1: Side-by-side comparison ===========
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Cycling NK (3 genes)
    ax1 = axes[0]
    data1 = [df[df["cluster_annot"] == c][z_col_short].values for c in cluster_order]
    parts1 = ax1.violinplot(data1, positions=range(len(cluster_order)), showmeans=True, showmedians=True)
    for pc in parts1['bodies']:
        pc.set_facecolor('#3498db')
        pc.set_alpha(0.6)
    ax1.set_xticks(range(len(cluster_order)))
    ax1.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax1.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax1.set_ylabel("Z-score")
    ax1.set_title(f"{sig_short}\n(MKI67, TYMS, STMN1)")
    ax1.set_ylim(-3, 5)
    
    # Proliferation-Extended (16 genes)
    ax2 = axes[1]
    data2 = [df[df["cluster_annot"] == c][z_col_extended].values for c in cluster_order]
    parts2 = ax2.violinplot(data2, positions=range(len(cluster_order)), showmeans=True, showmedians=True)
    for pc in parts2['bodies']:
        pc.set_facecolor('#e74c3c')
        pc.set_alpha(0.6)
    ax2.set_xticks(range(len(cluster_order)))
    ax2.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax2.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax2.set_ylabel("Z-score")
    ax2.set_title(f"{sig_extended}\n(16 proliferation markers)")
    ax2.set_ylim(-3, 5)
    
    fig.suptitle("Proliferation Score Distribution by Cluster", fontsize=14, fontweight='bold')
    plt.tight_layout()
    out_path1 = output_dir + "/proliferation_violin_comparison.png"
    fig.savefig(out_path1, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path1}")
    
    # =========== Plot 2: Overlay comparison ===========
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(cluster_order))
    width = 0.35
    
    means_short = [df[df["cluster_annot"] == c][z_col_short].mean() for c in cluster_order]
    means_extended = [df[df["cluster_annot"] == c][z_col_extended].mean() for c in cluster_order]
    
    bars1 = ax.bar(x - width/2, means_short, width, label=f'{sig_short} (3 genes)', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, means_extended, width, label=f'{sig_extended} (16 genes)', color='#e74c3c', alpha=0.8)
    
    ax.set_ylabel('Mean Z-score')
    ax.set_title('Proliferation Score Comparison by Cluster', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax.legend()
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    out_path2 = output_dir + "/proliferation_barplot_comparison.png"
    fig.savefig(out_path2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path2}")
    
    # =========== Plot 3: Single violin plot for Extended signature ===========
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Prepare data for seaborn-style violin
    plot_data = []
    for c in cluster_order:
        cluster_vals = df[df["cluster_annot"] == c][z_col_extended].values
        for v in cluster_vals:
            plot_data.append({"cluster": c, "z_score": v})
    plot_df = pd.DataFrame(plot_data)
    
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6', '#f39c12']
    
    violin_parts = ax.violinplot(
        [df[df["cluster_annot"] == c][z_col_extended].values for c in cluster_order],
        positions=range(len(cluster_order)),
        showmeans=True,
        showmedians=True
    )
    
    for i, pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.6)
    
    ax.set_xticks(range(len(cluster_order)))
    ax.set_xticklabels(cluster_order, rotation=30, ha='right', fontsize=11)
    ax.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax.axhline(1, linestyle=':', color='orange', linewidth=0.8, label='z=1 threshold')
    ax.set_ylabel("Proliferation Z-score", fontsize=12)
    ax.set_xlabel("")
    ax.set_title(f"Proliferation-Extended Score by Cluster\n(16 genes: MKI67, PCNA, TOP2A, MCM2/4/5/6, CCNB1, CDK1, AURKB, BIRC5, UBE2C, TYMS, STMN1, HMGB2, HIST1H4C)", 
                 fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_ylim(-3, 5)
    
    plt.tight_layout()
    out_path3 = output_dir + "/proliferation_extended_violin.png"
    fig.savefig(out_path3, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path3}")
    
    # Print summary
    print("\n" + "="*60)
    print("PROLIFERATION SCORE SUMMARY BY CLUSTER")
    print("="*60)
    print(f"\n{'Cluster':<20} {'Cycling NK (z)':<15} {'Extended (z)':<15} {'pct_z>=1':<10}")
    print("-"*60)
    for c in cluster_order:
        subset = df[df["cluster_annot"] == c]
        mean_short = subset[z_col_short].mean()
        mean_ext = subset[z_col_extended].mean()
        pct_ge1 = (subset[z_col_extended] >= 1).mean() * 100
        print(f"{c:<20} {mean_short:>12.3f} {mean_ext:>15.3f} {pct_ge1:>9.1f}%")
    print("="*60)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
