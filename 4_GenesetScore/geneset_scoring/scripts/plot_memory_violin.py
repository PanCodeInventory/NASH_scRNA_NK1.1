#!/usr/bin/env python
"""
plot_memory_violin.py

Generate violin plots comparing memory NK scores across clusters.
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
    output_dir = results_dir + "/plots/memory"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load cell-level scores
    print(f"Loading: {cell_scores_path}")
    df = pd.read_csv(cell_scores_path, sep="\t", compression="gzip")
    
    # Memory signatures
    sig_memory = "Memory-NK"
    sig_adaptive = "Adaptive-NK-CMV"
    
    z_col_memory = f"{sig_memory}__z"
    z_col_adaptive = f"{sig_adaptive}__z"
    
    # Cluster order
    cluster_order = ["Cytotoxic NK", "ILC1", "Circulating NK", "High-Ribosomal NK", "Homeostatic NK"]
    colors = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6', '#f39c12']
    
    # =========== Plot 1: Side-by-side comparison ===========
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Memory-NK
    ax1 = axes[0]
    data1 = [df[df["cluster_annot"] == c][z_col_memory].values for c in cluster_order]
    parts1 = ax1.violinplot(data1, positions=range(len(cluster_order)), showmeans=True, showmedians=True)
    for i, pc in enumerate(parts1['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.6)
    ax1.set_xticks(range(len(cluster_order)))
    ax1.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax1.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax1.set_ylabel("Z-score")
    ax1.set_title(f"{sig_memory}\n(Sell, Il7r, Cxcr3, Tcf7, Cd69, Klrg1)")
    ax1.set_ylim(-3, 5)
    
    # Adaptive-NK-CMV
    ax2 = axes[1]
    data2 = [df[df["cluster_annot"] == c][z_col_adaptive].values for c in cluster_order]
    parts2 = ax2.violinplot(data2, positions=range(len(cluster_order)), showmeans=True, showmedians=True)
    for i, pc in enumerate(parts2['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.6)
    ax2.set_xticks(range(len(cluster_order)))
    ax2.set_xticklabels(cluster_order, rotation=30, ha='right')
    ax2.axhline(0, linestyle='--', color='gray', linewidth=0.8)
    ax2.set_ylabel("Z-score")
    ax2.set_title(f"{sig_adaptive}\n(Klrc2, Fcgr3, Lilrb4a, B3gat1, Zeb2, Fasl)")
    ax2.set_ylim(-3, 5)
    
    fig.suptitle("Memory NK Score Distribution by Cluster", fontsize=14, fontweight='bold')
    plt.tight_layout()
    out_path1 = output_dir + "/memory_violin_comparison.png"
    fig.savefig(out_path1, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path1}")
    
    # =========== Plot 2: Barplot comparison ===========
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(cluster_order))
    width = 0.35
    
    means_memory = [df[df["cluster_annot"] == c][z_col_memory].mean() for c in cluster_order]
    means_adaptive = [df[df["cluster_annot"] == c][z_col_adaptive].mean() for c in cluster_order]
    
    bars1 = ax.bar(x - width/2, means_memory, width, label=f'{sig_memory}', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, means_adaptive, width, label=f'{sig_adaptive}', color='#e74c3c', alpha=0.8)
    
    ax.set_ylabel('Mean Z-score')
    ax.set_title('Memory NK Score Comparison by Cluster', fontweight='bold')
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
    out_path2 = output_dir + "/memory_barplot_comparison.png"
    fig.savefig(out_path2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path2}")
    
    # =========== Plot 3: Single violin plot for Memory-NK ===========
    fig, ax = plt.subplots(figsize=(12, 7))
    
    violin_parts = ax.violinplot(
        [df[df["cluster_annot"] == c][z_col_memory].values for c in cluster_order],
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
    ax.set_ylabel("Memory-NK Z-score", fontsize=12)
    ax.set_xlabel("")
    ax.set_title(f"Memory-NK Score by Cluster\n(Sell, Il7r, Cxcr3, Tcf7, Cd69, Klrg1)", 
                 fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_ylim(-3, 5)
    
    plt.tight_layout()
    out_path3 = output_dir + "/memory_nk_violin.png"
    fig.savefig(out_path3, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path3}")
    
    # =========== Plot 4: Single violin plot for Adaptive-NK-CMV ===========
    fig, ax = plt.subplots(figsize=(12, 7))
    
    violin_parts = ax.violinplot(
        [df[df["cluster_annot"] == c][z_col_adaptive].values for c in cluster_order],
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
    ax.set_ylabel("Adaptive-NK-CMV Z-score", fontsize=12)
    ax.set_xlabel("")
    ax.set_title(f"Adaptive-NK-CMV Score by Cluster\n(Klrc2, Fcgr3, Lilrb4a, B3gat1, Zeb2, Fasl)", 
                 fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_ylim(-3, 5)
    
    plt.tight_layout()
    out_path4 = output_dir + "/adaptive_nk_cmv_violin.png"
    fig.savefig(out_path4, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path4}")
    
    # Print summary
    print("\n" + "="*70)
    print("MEMORY NK SCORE SUMMARY BY CLUSTER")
    print("="*70)
    print(f"\n{'Cluster':<20} {'Memory-NK (z)':<15} {'Adaptive (z)':<15} {'Mem%':<10} {'Adp%':<10}")
    print("-"*70)
    for c in cluster_order:
        subset = df[df["cluster_annot"] == c]
        mean_mem = subset[z_col_memory].mean()
        mean_adp = subset[z_col_adaptive].mean()
        pct_mem = (subset[z_col_memory] >= 1).mean() * 100
        pct_adp = (subset[z_col_adaptive] >= 1).mean() * 100
        print(f"{c:<20} {mean_mem:>12.3f} {mean_adp:>15.3f} {pct_mem:>9.1f}% {pct_adp:>9.1f}%")
    print("="*70)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
