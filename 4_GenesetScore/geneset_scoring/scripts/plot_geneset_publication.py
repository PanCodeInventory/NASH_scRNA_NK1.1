#!/usr/bin/env python
"""
plot_geneset_publication.py

Publication-quality visualization for NK cell gene sets.
Each gene set gets its own independent figure.

Following scientific-visualization skill guidelines:
- Colorblind-safe Okabe-Ito palette
- Sans-serif fonts (Arial)
- Proper figure dimensions for journals
- PDF + high-DPI PNG export
- Remove unnecessary spines
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# ============================================================================
# Publication Style Configuration
# ============================================================================

def apply_publication_style():
    """Apply publication-quality matplotlib settings."""
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    mpl.rcParams['font.size'] = 8
    mpl.rcParams['axes.labelsize'] = 9
    mpl.rcParams['axes.titlesize'] = 10
    mpl.rcParams['xtick.labelsize'] = 7
    mpl.rcParams['ytick.labelsize'] = 7
    mpl.rcParams['legend.fontsize'] = 7
    
    mpl.rcParams['figure.dpi'] = 100
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams['savefig.bbox'] = 'tight'
    mpl.rcParams['savefig.pad_inches'] = 0.05
    
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['xtick.major.width'] = 0.5
    mpl.rcParams['ytick.major.width'] = 0.5
    
    mpl.rcParams['legend.frameon'] = False
    mpl.rcParams['legend.borderpad'] = 0.3

# Okabe-Ito colorblind-safe palette
OKABE_ITO = {
    'orange': '#E69F00',
    'sky_blue': '#56B4E9',
    'bluish_green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermillion': '#D55E00',
    'reddish_purple': '#CC79A7',
    'black': '#000000'
}

CLUSTER_COLORS = [
    OKABE_ITO['blue'],          # Cytotoxic NK
    OKABE_ITO['vermillion'],    # ILC1
    OKABE_ITO['sky_blue'],      # Circulating NK
    OKABE_ITO['reddish_purple'], # High-Ribosomal NK
    OKABE_ITO['bluish_green'],  # Homeostatic NK
]

CLUSTER_ORDER = ["Cytotoxic NK", "ILC1", "Circulating NK", "High-Ribosomal NK", "Homeostatic NK"]

# Gene set definitions (based on literature)
GENESETS = {
    "Proliferation-Extended": {
        "genes": "Mki67, Pcna, Top2a, Mcm2/4/5/6, Ccnb1, Cdk1, Aurkb, Birc5, Ube2c, Tyms, Stmn1, Hmgb2, Hist1h4c",
        "short": "Proliferation",
        "description": "Cell cycle and proliferation markers"
    },
    "Memory-NK": {
        "genes": "Sell, Il7r, Cxcr3, Tcf7, Cd69, Klrg1",
        "short": "Memory",
        "description": "Memory NK cell signature (tissue-resident memory)"
    },
    "Immature-NK": {
        "genes": "Cd27, Tcf7, Lef1, Sell, Zbtb16, Cxcr3",
        "short": "Immature",
        "description": "Immature NK cell markers (CD11b-CD27+ stage)"
    },
    "Regulatory-NK": {
        "genes": "Il10, Klrc1, Tgfb1, Il2ra, Ctla4, Pdcd1, Havcr2, Lag3, Tigit",
        "short": "Regulatory",
        "description": "Regulatory/exhaustion markers"
    }
}

GENESET_ORDER = ["Proliferation-Extended", "Memory-NK", "Immature-NK", "Regulatory-NK"]


def save_publication_figure(fig, name, output_dir, formats=['pdf', 'png']):
    """Save figure in publication formats."""
    for fmt in formats:
        path = os.path.join(output_dir, f"{name}.{fmt}")
        if fmt == 'pdf':
            fig.savefig(path, format='pdf', bbox_inches='tight')
        else:
            fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {path}")


def create_single_violin(df, geneset_name, output_dir):
    """Create publication-quality violin plot for a single gene set."""
    z_col = f"{geneset_name}__z"
    geneset_info = GENESETS[geneset_name]
    
    fig_width = 3.5
    fig_height = 2.8
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    data = [df[df["cluster_annot"] == c][z_col].values for c in CLUSTER_ORDER]
    positions = np.arange(len(CLUSTER_ORDER))
    
    violin_parts = ax.violinplot(data, positions=positions, showmeans=True, showmedians=False, widths=0.7)
    
    for i, (pc, color) in enumerate(zip(violin_parts['bodies'], CLUSTER_COLORS)):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_linewidth(0.5)
        pc.set_alpha(0.7)
    
    violin_parts['cmeans'].set_color('white')
    violin_parts['cmeans'].set_linewidth(1.5)
    
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, zorder=0)
    ax.axhline(1, color='orange', linestyle=':', linewidth=0.5, zorder=0, alpha=0.7)
    ax.axhline(-1, color='orange', linestyle=':', linewidth=0.5, zorder=0, alpha=0.7)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(CLUSTER_ORDER, rotation=30, ha='right')
    ax.set_ylabel('Z-score')
    ax.set_xlabel('')
    ax.set_ylim(-3, 5)
    ax.set_title(geneset_info['short'], fontweight='bold', fontsize=10)
    
    ax.yaxis.grid(True, linestyle='-', linewidth=0.3, alpha=0.3)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    save_publication_figure(fig, f"{geneset_name.replace('-', '_')}_violin", output_dir)
    plt.close(fig)


def create_barplot(df, geneset_name, output_dir):
    """Create publication-quality bar plot for a single gene set."""
    z_col = f"{geneset_name}__z"
    geneset_info = GENESETS[geneset_name]
    
    fig_width = 3.5
    fig_height = 2.5
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    means = []
    sems = []
    for c in CLUSTER_ORDER:
        vals = df[df["cluster_annot"] == c][z_col].values
        means.append(np.mean(vals))
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)))
    
    x = np.arange(len(CLUSTER_ORDER))
    bars = ax.bar(x, means, yerr=sems, color=CLUSTER_COLORS, edgecolor='black', 
                  linewidth=0.5, capsize=2, error_kw={'linewidth': 0.5})
    
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(CLUSTER_ORDER, rotation=30, ha='right')
    ax.set_ylabel('Mean Z-score ± SEM')
    ax.set_xlabel('')
    ax.set_ylim(-1, 1)
    ax.set_title(geneset_info['short'], fontweight='bold', fontsize=10)
    
    ax.yaxis.grid(True, linestyle='-', linewidth=0.3, alpha=0.3)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    save_publication_figure(fig, f"{geneset_name.replace('-', '_')}_barplot", output_dir)
    plt.close(fig)


def create_combined_figure(df, output_dir):
    """Create multi-panel figure with all gene sets."""
    fig_width = 7.2
    fig_height = 5.0
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.35,
                          left=0.10, right=0.98, top=0.92, bottom=0.15)
    
    axes = [
        fig.add_subplot(gs[0, 0]),  # Proliferation
        fig.add_subplot(gs[0, 1]),  # Memory
        fig.add_subplot(gs[1, 0]),  # Immature
        fig.add_subplot(gs[1, 1]),  # Regulatory
    ]
    
    panel_labels = ['A', 'B', 'C', 'D']
    
    for ax, label, geneset_name in zip(axes, panel_labels, GENESET_ORDER):
        z_col = f"{geneset_name}__z"
        geneset_info = GENESETS[geneset_name]
        
        data = [df[df["cluster_annot"] == c][z_col].values for c in CLUSTER_ORDER]
        positions = np.arange(len(CLUSTER_ORDER))
        
        violin_parts = ax.violinplot(data, positions=positions, showmeans=True, showmedians=False, widths=0.7)
        
        for i, (pc, color) in enumerate(zip(violin_parts['bodies'], CLUSTER_COLORS)):
            pc.set_facecolor(color)
            pc.set_edgecolor('black')
            pc.set_linewidth(0.5)
            pc.set_alpha(0.7)
        
        violin_parts['cmeans'].set_color('white')
        violin_parts['cmeans'].set_linewidth(1.5)
        
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.5, zorder=0)
        
        ax.set_xticks(positions)
        ax.set_xticklabels(CLUSTER_ORDER, rotation=30, ha='right', fontsize=6)
        ax.set_ylabel('Z-score', fontsize=8)
        ax.set_ylim(-3, 5)
        ax.set_title(geneset_info['short'], fontweight='bold', fontsize=9)
        
        ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=11, fontweight='bold', va='top')
        
        ax.yaxis.grid(True, linestyle='-', linewidth=0.3, alpha=0.3)
        ax.set_axisbelow(True)
    
    fig.suptitle('NK Cell Gene Set Scores by Cluster', fontsize=11, fontweight='bold', y=0.98)
    
    save_publication_figure(fig, "geneset_combined_figure", output_dir)
    plt.close(fig)


def create_heatmap(df, output_dir):
    """Create publication-quality heatmap of all gene set scores."""
    heatmap_data = []
    for geneset in GENESET_ORDER:
        z_col = f"{geneset}__z"
        row = []
        for cluster in CLUSTER_ORDER:
            mean_z = df[df["cluster_annot"] == cluster][z_col].mean()
            row.append(mean_z)
        heatmap_data.append(row)
    
    heatmap_df = pd.DataFrame(
        heatmap_data,
        index=[GENESETS[g]['short'] for g in GENESET_ORDER],
        columns=CLUSTER_ORDER
    )
    
    fig_width = 4.5
    fig_height = 3.2
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    im = ax.imshow(heatmap_df.values, cmap='RdBu_r', aspect='auto', vmin=-0.6, vmax=0.6)
    
    ax.set_xticks(np.arange(len(CLUSTER_ORDER)))
    ax.set_xticklabels(CLUSTER_ORDER, rotation=30, ha='right')
    ax.set_yticks(np.arange(len(GENESET_ORDER)))
    ax.set_yticklabels([GENESETS[g]['short'] for g in GENESET_ORDER])
    
    for i in range(len(GENESET_ORDER)):
        for j in range(len(CLUSTER_ORDER)):
            val = heatmap_df.values[i, j]
            color = 'white' if abs(val) > 0.3 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', color=color, fontsize=7, fontweight='bold')
    
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Mean Z-score', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    
    ax.set_title('Gene Set Scores Heatmap', fontweight='bold', fontsize=10)
    
    plt.tight_layout()
    save_publication_figure(fig, "geneset_heatmap", output_dir)
    plt.close(fig)


def print_summary(df):
    """Print summary statistics for all gene sets."""
    print("\n" + "="*80)
    print("GENE SET SCORES SUMMARY")
    print("="*80)
    
    for geneset in GENESET_ORDER:
        z_col = f"{geneset}__z"
        info = GENESETS[geneset]
        
        print(f"\n{info['short']} ({geneset})")
        print(f"  Genes: {info['genes']}")
        print(f"  Description: {info['description']}")
        print(f"  {'Cluster':<20} {'Mean Z':>10} {'% z≥1':>10} {'n':>8}")
        print("  " + "-"*50)
        
        for cluster in CLUSTER_ORDER:
            vals = df[df["cluster_annot"] == cluster][z_col].values
            mean_z = np.mean(vals)
            pct_ge1 = (vals >= 1).mean() * 100
            n = len(vals)
            print(f"  {cluster:<20} {mean_z:>10.3f} {pct_ge1:>9.1f}% {n:>8}")
    
    print("="*80)


def main():
    apply_publication_style()
    
    script_dir = Path(__file__).parent
    results_dir = script_dir.parent / "results"
    cell_scores_path = results_dir / "cell_scores.tsv.gz"
    output_dir = results_dir / "plots" / "publication"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading: {cell_scores_path}")
    df = pd.read_csv(cell_scores_path, sep="\t", compression="gzip")
    
    print(f"\nOutput directory: {output_dir}")
    print("="*60)
    
    # Create individual plots
    print("\n[1/4] Creating individual violin and bar plots...")
    for geneset in GENESET_ORDER:
        print(f"  Processing: {geneset}")
        create_single_violin(df, geneset, output_dir)
        create_barplot(df, geneset, output_dir)
    
    # Create combined figure
    print("\n[2/4] Creating combined multi-panel figure...")
    create_combined_figure(df, output_dir)
    
    # Create heatmap
    print("\n[3/4] Creating heatmap...")
    create_heatmap(df, output_dir)
    
    # Print summary
    print("\n[4/4] Summary statistics...")
    print_summary(df)
    
    print(f"\n✓ All figures saved to: {output_dir}")


if __name__ == "__main__":
    main()
