#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(stringr)
})

# Paths
rds_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.fixed.rds"
proportions_path <- "3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv"
plots_dir <- "3_Analysis/1.ClusterAnalysis/plots"
backup_dir <- file.path(plots_dir, "backup")

# Create directories
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters
hidden_clusters <- "6"  # Hide B cell contamination cluster
formats <- c("png", "pdf")
width <- 10
height <- 8
dpi <- 300

cat("Starting visualization with hidden clusters...\n")

# Load the renamed Seurat object
cat(sprintf("Loading Seurat object: %s\n", rds_path))
obj <- readRDS(rds_path)
cat(sprintf("Loaded object with %d cells\n", ncol(obj)))

# Function to backup existing plots
backup_existing_plots <- function() {
  existing_plots <- list.files(plots_dir, pattern = "\\.(png|pdf)$", full.names = TRUE)
  if (length(existing_plots) > 0) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    for (plot_file in existing_plots) {
      backup_name <- file.path(backup_dir, basename(plot_file) %>% str_replace("\\.(png|pdf)$", sprintf("_%s.\\1", timestamp)))
      file.copy(plot_file, backup_name)
    }
    cat(sprintf("Backed up %d existing plots\n", length(existing_plots)))
  }
}

# Function to generate UMAP plots
generate_umap_plots <- function(seurat_obj, hidden_clusters) {
  cat("Generating UMAP plots...\n")
  
  # Filter out hidden clusters
  cells_to_keep <- !seurat_obj$seurat_clusters %in% hidden_clusters
  filtered_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
  
  # UMAP by timepoint
  p1 <- DimPlot(filtered_obj, reduction = "umap", group.by = "timepoint", 
                label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("UMAP by Timepoint (Cluster 6 Hidden)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # UMAP by cluster
  p2 <- DimPlot(filtered_obj, reduction = "umap", group.by = "seurat_clusters", 
                label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("UMAP by Cluster (Cluster 6 Hidden)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # UMAP by SingleR labels
  if ("SingleR_label" %in% colnames(filtered_obj@meta.data)) {
    p3 <- DimPlot(filtered_obj, reduction = "umap", group.by = "SingleR_label", 
                  label = TRUE, repel = TRUE, label.size = 3) +
      ggtitle("UMAP by SingleR Label (Cluster 6 Hidden)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14)) +
      theme(legend.position = "bottom")
  } else {
    p3 <- DimPlot(filtered_obj, reduction = "umap", group.by = "ident", 
                  label = TRUE, repel = TRUE, label.size = 3) +
      ggtitle("UMAP by Identity (Cluster 6 Hidden)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14))
  }
  
  return(list(by_timepoint = p1, by_cluster = p2, by_singler = p3))
}

# Function to generate cluster proportion plots
generate_proportion_plots <- function(proportions_df, hidden_clusters) {
  cat("Generating cluster proportion plots...\n")
  
  # Filter out hidden clusters
  filtered_proportions <- proportions_df %>% 
    filter(!cluster %in% hidden_clusters)
  
  # Convert cluster to factor for proper ordering
  filtered_proportions$cluster <- factor(filtered_proportions$cluster)
  filtered_proportions$timepoint <- factor(filtered_proportions$timepoint, 
                                          levels = c("0W_NCD", "1W_MCD", "2W_MCD", "6W_MCD"))
  
  # Line plot of cluster proportions over time
  p1 <- ggplot(filtered_proportions, aes(x = timepoint, y = proportion, 
                                         color = factor(cluster), group = factor(cluster))) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    labs(title = "Cluster Proportions Over Time (Cluster 6 Hidden)",
         x = "Timepoint", y = "Proportion", color = "Cluster") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_color_brewer(palette = "Set2")
  
  # Stacked bar plot
  p2 <- ggplot(filtered_proportions, aes(x = timepoint, y = proportion, 
                                         fill = factor(cluster))) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Cluster Composition Over Time (Cluster 6 Hidden)",
         x = "Timepoint", y = "Proportion", fill = "Cluster") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(labels = scales::percent)
  
  # Absolute count stacked bar plot
  p3 <- ggplot(filtered_proportions, aes(x = timepoint, y = count, 
                                         fill = factor(cluster))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Cluster Cell Counts Over Time (Cluster 6 Hidden)",
         x = "Timepoint", y = "Cell Count", fill = "Cluster") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_brewer(palette = "Set2")
  
  return(list(lineplot = p1, stacked_percent = p2, stacked_count = p3))
}

# Function to save plots in multiple formats
save_plots <- function(plot_list, plot_names, formats, width, height, dpi) {
  for (i in seq_along(plot_list)) {
    plot_name <- plot_names[i]
    plot_obj <- plot_list[[i]]
    
    for (format in formats) {
      filename <- file.path(plots_dir, sprintf("%s_hidden_cluster6.%s", plot_name, format))
      
      if (format == "png") {
        ggsave(filename, plot = plot_obj, width = width, height = height, 
               dpi = dpi, bg = "white")
      } else if (format == "pdf") {
        ggsave(filename, plot = plot_obj, width = width, height = height, 
               device = "pdf")
      }
      
      cat(sprintf("Saved: %s\n", filename))
    }
  }
}

# Backup existing plots
backup_existing_plots()

# Generate UMAP plots
umap_plots <- generate_umap_plots(obj, hidden_clusters)

# Generate proportion plots
if (file.exists(proportions_path)) {
  proportions_df <- read_csv(proportions_path, show_col_types = FALSE)
  proportion_plots <- generate_proportion_plots(proportions_df, hidden_clusters)
} else {
  cat("Warning: Proportions file not found, skipping proportion plots\n")
  proportion_plots <- list()
}

# Save all plots
cat("Saving plots...\n")

# Save UMAP plots
umap_plot_names <- c("UMAP_timepoint", "UMAP_cluster", "UMAP_singler")
save_plots(umap_plots, umap_plot_names, formats, width, height, dpi)

# Save proportion plots
if (length(proportion_plots) > 0) {
  prop_plot_names <- c("cluster_proportion_lineplot", "cluster_composition_stacked_percent", 
                       "cluster_composition_stacked_count")
  save_plots(proportion_plots, prop_plot_names, formats, width, height, dpi)
}

# Create combined plots
if (length(proportion_plots) > 0) {
  combined_plot <- (umap_plots$by_cluster / proportion_plots$lineplot) + 
    plot_layout(heights = c(2, 1))
  
  for (format in formats) {
    filename <- file.path(plots_dir, sprintf("combined_UMAP_proportions_hidden_cluster6.%s", format))
    if (format == "png") {
      ggsave(filename, plot = combined_plot, width = width, height = height * 1.5, 
             dpi = dpi, bg = "white")
    } else if (format == "pdf") {
      ggsave(filename, plot = combined_plot, width = width, height = height * 1.5, 
             device = "pdf")
    }
    cat(sprintf("Saved combined plot: %s\n", filename))
  }
}

# Generate visualization report
viz_report <- sprintf(
  "# Visualization Report (Hidden Cluster 6)\n\n## Strategy\n- **Hidden Cluster**: 6 (B cell contamination)\n- **Visible Clusters**: 0, 1, 2, 3, 4, 5\n- **Original Assignments**: \n  - Cluster 5 (B cell contamination) → Cluster 6 (hidden)\n  - Cluster 6 (proliferation) → Cluster 5 (visible)\n\n## Generated Plots\n\n### UMAP Visualizations\n1. **UMAP_timepoint_hidden_cluster6**: UMAP colored by timepoint\n2. **UMAP_cluster_hidden_cluster6**: UMAP colored by cluster\n3. **UMAP_singler_hidden_cluster6**: UMAP colored by SingleR labels\n\n### Proportion Analysis\n1. **cluster_proportion_lineplot_hidden_cluster6**: Line plot of cluster proportions over time\n2. **cluster_composition_stacked_percent_hidden_cluster6**: Stacked percentage bar plot\n3. **cluster_composition_stacked_count_hidden_cluster6**: Stacked count bar plot\n\n### Combined Visualizations\n1. **combined_UMAP_proportions_hidden_cluster6**: Combined UMAP and proportion plot\n\n## File Locations\n- **Plots Directory**: `%s`\n- **Backup Directory**: `%s`\n- **Data Source**: `%s`\n\n## Format Options\n- **PNG**: High-resolution raster images (300 DPI)\n- **PDF**: Vector graphics for publication quality\n\n## Quality Notes\n- All plots exclude cluster 6 (B cell contamination)\n- Cluster numbering reflects corrected assignments\n- Timepoints are ordered: 0W_NCD → 1W_MCD → 2W_MCD → 6W_MCD\n- Color schemes use ColorBrewer palettes for accessibility\n\n---\n*Generated: %s*",
  plots_dir,
  backup_dir,
  rds_path,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(viz_report, file.path("2_DataProcessing/reports", "visualization_report_hidden_cluster6.md"))
cat("Visualization generation completed successfully!\n")
cat(sprintf("Visualization report: %s\n", file.path("2_DataProcessing/reports", "visualization_report_hidden_cluster6.md")))
cat(sprintf("All plots saved to: %s\n", plots_dir))