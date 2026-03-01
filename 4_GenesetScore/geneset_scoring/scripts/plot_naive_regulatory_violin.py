#!/usr/bin/env python
"""
plot_naive_regulatory_violin.py

Generate violin plots comparing Naive-NK and Regulatory-NK scores across clusters.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.dirname(script_dir) + "/results"
    cell_scores_path = results_dir + "/cell_scores.tsv.gz"
    output_dir = results_dir + "/plots/naive_regulatory"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load cell-level scores
    print(f"Loading: {cell_scores_path}")
    df = pd.read_csv(cell_scores_path, sep="\t", compression="gzip")
    
    # Signatures
    sig_naive = "Naive-NK"
    sig_regulatory = "Regulatory-NK"
    
    z_col_naive = f"{sig_naive}__z"
    z_col_regulatory = f"{sig_regulatory}__z"
    
    # Cluster order
    cluster_order = ["Cytotoxic NK", "ILC1", "Circulating NK", "High-Ribosomal NK", "Homeostatic NK"]
    colors = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6', '#f39c12']
    
    # =========== Plot 1: Side-by-side comparison ===========
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Naive-NK
    ax1 = axes[0]
    data1 = [df[df["cluster_annot"] == c][z_col_naive].values for c in cluster_order]
    parts1 = ax1.violinplot(data1, positions=range(len(cluster_order)), showmeans=True, showmedians=True)
    for i, pc in enumerate(parts1['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.6)
    ax1.set_xticks(range(len(cluster_order)))
    ax1.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax1.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax1.set_ylabel("Z-score")
    ax1.set_title(f"{sig_naive}\n(Sell, Il7r, Cd27, Tcf7, Lef1, Ccr7)")
    ax1.set_ylim(-3, 5)
    
    # Regulatory-NK
    ax2 = axes[1]
    data2 = [df[df["cluster_annot"] == c][z_col_regulatory].values for c in cluster_order]
    parts2 = ax2.violinplot(data2, positions=range(len(cluster_order)), showmeans=True, showmedians=True)
    for i, pc in enumerate(parts2['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.6)
    ax2.set_xticks(range(len(cluster_order)))
    ax2.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax2.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax2.set_ylabel("Z-score")
    ax2.set_title(f"{sig_regulatory}\n(Il10, Klrc1, Tgfb1, Il2ra, Ctla4, Pdcd1, Havcr2, Lag3, Tigit)")
    ax2.set_ylim(-3, 5)
    
    fig.suptitle("Naive & Regulatory NK Score Distribution by Cluster", fontsize=14, fontweight='bold')
    plt.tight_layout()
    out_path1 = output_dir + "/naive_regulatory_violin_comparison.png"
    fig.savefig(out_path1, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path1}")
    
    # =========== Plot 2: Barplot comparison ===========
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(cluster_order))
    width = 0.35
    
    means_naive = [df[df["cluster_annot"] == c][z_col_naive].mean() for c in cluster_order]
    means_regulatory = [df[df["cluster_annot"] == c][z_col_regulatory].mean() for c in cluster_order]
    
    bars1 = ax.bar(x - width/2, means_naive, width, label=f'{sig_naive}', color='#27ae60', alpha=0.8)
    bars2 = ax.bar(x + width/2, means_regulatory, width, label=f'{sig_regulatory}', color='#8e44ad', alpha=0.8)
    
    ax.set_ylabel('Mean Z-score')
    ax.set_title('Naive-NK vs Regulatory-NK Score Comparison by Cluster', fontweight='bold')
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
                    ha='center', va='bottom', fontsize=9)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    out_path2 = output_dir + "/naive_regulatory_barplot_comparison.png"
    fig.savefig(out_path2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path2}")
    
    # =========== Plot 3: Combined heatmap-style visualization ===========
    # Create a summary heatmap
    fig, ax = plt.subplots(figsize=(10, 6))
    
    summary_data = []
    for c in cluster_order:
        subset = df[df["cluster_annot"] == c]
        summary_data.append([
            subset[z_col_naive].mean(),
            subset[z_col_regulatory].mean()
        ])
    
    summary_df = pd.DataFrame(summary_data, columns=['Naive-NK', 'Regulatory-NK'], index=cluster_order)
    
    im = ax.imshow(summary_df.values, cmap='RdYlBu_r', aspect='auto', vmin=-0.5, vmax=0.5)
    ax.set_xticks(range(len(summary_df.columns)))
    ax.set_xticklabels(summary_df.columns, fontsize=11)
    ax.set_yticks(range(len(summary_df.index)))
    ax.set_yticklabels(summary_df.index, fontsize=11)
    
    # Add value annotations
    for i in range(len(summary_df.index)):
        for j in range(len(summary_df.columns)):
            val = summary_df.values[i, j]
            color = 'white' if abs(val) > 0.3 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', color=color, fontsize=12, fontweight='bold')
    
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Mean Z-score')
    ax.set_title('Naive & Regulatory NK Scores Heatmap', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    out_path3 = output_dir + "/naive_regulatory_heatmap.png"
    fig.savefig(out_path3, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path3}")
    
    # Print summary
    print("\n" + "="*75)
    print("NAIVE & REGULATORY NK SCORE SUMMARY BY CLUSTER")
    print("="*75)
    print(f"\n{'Cluster':<20} {'Naive-NK (z)':<15} {'Regulatory (z)':<15} {'Naive%':<10} {'Reg%':<10}")
    print("-"*75)
    for c in cluster_order:
        subset = df[df["cluster_annot"] == c]
        mean_naive = subset[z_col_naive].mean()
        mean_reg = subset[z_col_regulatory].mean()
        pct_naive = (subset[z_col_naive] >= 1).mean() * 100
        pct_reg = (subset[z_col_regulatory] >= 1).mean() * 100
        print(f"{c:<20} {mean_naive:>12.3f} {mean_reg:>15.3f} {pct_naive:>9.1f}% {pct_reg:>9.1f}%")
    print("="*75)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
