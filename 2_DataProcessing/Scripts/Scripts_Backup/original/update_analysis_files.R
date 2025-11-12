#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Paths
markers_all_path <- "3_Analysis/1.ClusterAnalysis/data/markers_all_clusters.csv"
markers_top10_path <- "3_Analysis/1.ClusterAnalysis/data/markers_top10_per_cluster.csv"
proportions_path <- "3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv"
reports_dir <- "2_DataProcessing/reports"
backup_dir <- file.path(reports_dir, "backup")

# Create backup directory
dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)

# Cluster renaming mapping: 5→6 (B cell contamination), 6→5 (proliferation)
cluster_mapping <- c(
  "0" = "0",
  "1" = "1", 
  "2" = "2",
  "3" = "3",
  "4" = "4",
  "5" = "6",  # Original cluster 5 (B cell contamination) → cluster 6
  "6" = "5"   # Original cluster 6 (proliferation) → cluster 5
)

# Function to backup file
backup_file <- function(path) {
  if (file.exists(path)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    backup_name <- file.path(backup_dir, basename(path) %>% str_replace("\\.csv$", sprintf("_%s.csv", timestamp)))
    file.copy(path, backup_name)
    cat(sprintf("Backed up: %s -> %s\n", path, backup_name))
    return(backup_name)
  }
  return(NULL)
}

# Function to apply cluster mapping to a dataframe
apply_cluster_mapping <- function(df, cluster_col) {
  df[[cluster_col]] <- as.character(df[[cluster_col]])
  
  # Apply mapping in correct order to avoid conflicts
  # First rename cluster 6 to 5, then cluster 5 to 6
  df[[cluster_col]][df[[cluster_col]] == "6"] <- "temp_cluster_6"
  df[[cluster_col]][df[[cluster_col]] == "5"] <- "6"
  df[[cluster_col]][df[[cluster_col]] == "temp_cluster_6"] <- "5"
  
  df[[cluster_col]] <- as.integer(df[[cluster_col]])
  return(df)
}

cat("Starting to update analysis files with cluster renaming...\n")

# Backup original files
cat("Creating backups...\n")
backup_file(markers_all_path)
backup_file(markers_top10_path)
backup_file(proportions_path)

# Update markers_all_clusters.csv
cat("Updating markers_all_clusters.csv...\n")
if (file.exists(markers_all_path)) {
  markers_all <- read_csv(markers_all_path, show_col_types = FALSE)
  markers_all_updated <- apply_cluster_mapping(markers_all, "cluster")
  write_csv(markers_all_updated, markers_all_path)
  cat(sprintf("Updated %d marker entries\n", nrow(markers_all_updated)))
}

# Update markers_top10_per_cluster.csv
cat("Updating markers_top10_per_cluster.csv...\n")
if (file.exists(markers_top10_path)) {
  markers_top10 <- read_csv(markers_top10_path, show_col_types = FALSE)
  markers_top10_updated <- apply_cluster_mapping(markers_top10, "cluster")
  write_csv(markers_top10_updated, markers_top10_path)
  cat(sprintf("Updated %d top10 marker entries\n", nrow(markers_top10_updated)))
}

# Update cluster_proportions_by_timepoint.csv
cat("Updating cluster_proportions_by_timepoint.csv...\n")
if (file.exists(proportions_path)) {
  proportions <- read_csv(proportions_path, show_col_types = FALSE)
  proportions_updated <- apply_cluster_mapping(proportions, "cluster")
  write_csv(proportions_updated, proportions_path)
  cat(sprintf("Updated %d proportion entries\n", nrow(proportions_updated)))
}

# Generate update report
update_report <- sprintf(
  "# Analysis Files Update Report\n\n## Files Updated\n1. **markers_all_clusters.csv**: All cluster markers with renamed cluster numbers\n2. **markers_top10_per_cluster.csv**: Top 10 markers per cluster with renamed cluster numbers\n3. **cluster_proportions_by_timepoint.csv**: Cluster proportions by timepoint with renamed cluster numbers\n\n## Cluster Mapping Applied\n- **Cluster 5 → Cluster 6**: B cell contamination cluster (will be hidden in visualizations)\n- **Cluster 6 → Cluster 5**: Proliferation cluster\n- **Clusters 0-4**: No change\n\n## Backup Location\nAll original files have been backed up to: `%s`\n\n## Next Steps\n1. Generate visualizations with hidden cluster 6\n2. Validate that all downstream analyses use the updated cluster numbers\n3. Update documentation to reflect the new cluster numbering\n\n---\n*Generated: %s*",
  backup_dir,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(update_report, file.path(reports_dir, "analysis_files_update_report.md"))
cat("Analysis files update completed successfully!\n")
cat(sprintf("Update report: %s\n", file.path(reports_dir, "analysis_files_update_report.md")))