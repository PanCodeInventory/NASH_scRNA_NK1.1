# Visualization Report (Hidden Cluster 6)

## Strategy
- **Hidden Cluster**: 6 (B cell contamination)
- **Visible Clusters**: 0, 1, 2, 3, 4, 5
- **Original Assignments**: 
  - Cluster 5 (B cell contamination) → Cluster 6 (hidden)
  - Cluster 6 (proliferation) → Cluster 5 (visible)

## Generated Plots

### UMAP Visualizations
1. **UMAP_timepoint_hidden_cluster6**: UMAP colored by timepoint
2. **UMAP_cluster_hidden_cluster6**: UMAP colored by cluster
3. **UMAP_singler_hidden_cluster6**: UMAP colored by SingleR labels

### Proportion Analysis
1. **cluster_proportion_lineplot_hidden_cluster6**: Line plot of cluster proportions over time
2. **cluster_composition_stacked_percent_hidden_cluster6**: Stacked percentage bar plot
3. **cluster_composition_stacked_count_hidden_cluster6**: Stacked count bar plot

### Combined Visualizations
1. **combined_UMAP_proportions_hidden_cluster6**: Combined UMAP and proportion plot

## File Locations
- **Plots Directory**: `3_Analysis/1.ClusterAnalysis/plots`
- **Backup Directory**: `3_Analysis/1.ClusterAnalysis/plots/backup`
- **Data Source**: `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.fixed.rds`

## Format Options
- **PNG**: High-resolution raster images (300 DPI)
- **PDF**: Vector graphics for publication quality

## Quality Notes
- All plots exclude cluster 6 (B cell contamination)
- Cluster numbering reflects corrected assignments
- Timepoints are ordered: 0W_NCD → 1W_MCD → 2W_MCD → 6W_MCD
- Color schemes use ColorBrewer palettes for accessibility

---
*Generated: 2025-10-23 20:03:52*
