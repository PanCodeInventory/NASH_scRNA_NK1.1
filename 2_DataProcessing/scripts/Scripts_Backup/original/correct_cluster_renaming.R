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

cat("Starting correct cluster renaming...\n")

# Backup current files
cat("Creating backups...\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
file.copy(markers_all_path, file.path(backup_dir, sprintf("markers_all_clusters_%s.csv", timestamp)))
file.copy(markers_top10_path, file.path(backup_dir, sprintf("markers_top10_per_cluster_%s.csv", timestamp)))
file.copy(proportions_path, file.path(backup_dir, sprintf("cluster_proportions_%s.csv", timestamp)))

# Function to correctly apply cluster mapping
apply_correct_cluster_mapping <- function(df, cluster_col) {
  df[[cluster_col]] <- as.character(df[[cluster_col]])
  
  # Step 1: Identify B cell contamination markers (original cluster 5)
  b_cell_genes <- c("Iglc3", "Cd79a", "Iglc2", "Ebf1", "Ms4a1", "Cd79b", "Ighm", "Iglc1", "Ly6d", "Pax5")
  
  # Step 2: Identify proliferation markers (original cluster 6)  
  prolif_genes <- c("H2afx", "Cdca3", "Hist1h3c", "H2afv", "Hmgn2", "Knl1", "Hist1h2ab", "Mki67", "Ezh2", "Tuba1c")
  
  # Step 3: Apply renaming based on gene patterns
  # Original cluster 5 (B cell) → cluster 6
  # Original cluster 6 (proliferation) → cluster 5
  
  # For markers files, we need to check gene names
  if ("gene" %in% colnames(df)) {
    # Cluster 5 entries with B cell genes → cluster 6
    b_cell_mask <- df[[cluster_col]] == "5" & df$gene %in% b_cell_genes
    df[[cluster_col]][b_cell_mask] <- "6"
    
    # Cluster 6 entries with proliferation genes → cluster 5  
    prolif_mask <- df[[cluster_col]] == "6" & df$gene %in% prolif_genes
    df[[cluster_col]][prolif_mask] <- "5"
    
    # Any remaining cluster 5 entries (shouldn't be many) → cluster 6
    remaining_cluster5 <- df[[cluster_col]] == "5"
    df[[cluster_col]][remaining_cluster5] <- "6"
    
    # Any remaining cluster 6 entries (shouldn't be many) → cluster 5
    remaining_cluster6 <- df[[cluster_col]] == "6"
    df[[cluster_col]][remaining_cluster6] <- "5"
  } else {
    # For proportions file, just swap the numbers
    df[[cluster_col]][df[[cluster_col]] == "5"] <- "temp_5"
    df[[cluster_col]][df[[cluster_col]] == "6"] <- "5"
    df[[cluster_col]][df[[cluster_col]] == "temp_5"] <- "6"
  }
  
  df[[cluster_col]] <- as.integer(df[[cluster_col]])
  return(df)
}

# Update markers_all_clusters.csv
cat("Updating markers_all_clusters.csv...\n")
if (file.exists(markers_all_path)) {
  markers_all <- read_csv(markers_all_path, show_col_types = FALSE)
  markers_all_updated <- apply_correct_cluster_mapping(markers_all, "cluster")
  write_csv(markers_all_updated, markers_all_path)
  cat(sprintf("Updated %d marker entries\n", nrow(markers_all_updated)))
  
  # Verify the changes
  cat("Verification:\n")
  cluster5_genes <- markers_all_updated %>% filter(cluster == 5) %>% pull(gene) %>% head(5)
  cluster6_genes <- markers_all_updated %>% filter(cluster == 6) %>% pull(gene) %>% head(5)
  cat(sprintf("Cluster 5 sample genes: %s\n", paste(cluster5_genes, collapse = ", ")))
  cat(sprintf("Cluster 6 sample genes: %s\n", paste(cluster6_genes, collapse = ", ")))
}

# Update markers_top10_per_cluster.csv
cat("Updating markers_top10_per_cluster.csv...\n")
if (file.exists(markers_top10_path)) {
  markers_top10 <- read_csv(markers_top10_path, show_col_types = FALSE)
  markers_top10_updated <- apply_correct_cluster_mapping(markers_top10, "cluster")
  write_csv(markers_top10_updated, markers_top10_path)
  cat(sprintf("Updated %d top10 marker entries\n", nrow(markers_top10_updated)))
  
  # Verify the changes
  cat("Verification:\n")
  cluster5_top10 <- markers_top10_updated %>% filter(cluster == 5) %>% pull(gene)
  cluster6_top10 <- markers_top10_updated %>% filter(cluster == 6) %>% pull(gene)
  cat(sprintf("Cluster 5 top10 genes: %s\n", paste(cluster5_top10, collapse = ", ")))
  cat(sprintf("Cluster 6 top10 genes: %s\n", paste(cluster6_top10, collapse = ", ")))
}

# Update cluster_proportions_by_timepoint.csv
cat("Updating cluster_proportions_by_timepoint.csv...\n")
if (file.exists(proportions_path)) {
  proportions <- read_csv(proportions_path, show_col_types = FALSE)
  proportions_updated <- apply_correct_cluster_mapping(proportions, "cluster")
  write_csv(proportions_updated, proportions_path)
  cat(sprintf("Updated %d proportion entries\n", nrow(proportions_updated)))
  
  # Verify the changes
  cat("Verification:\n")
  cluster_counts <- proportions_updated %>% count(cluster) %>% arrange(cluster)
  for (i in 1:nrow(cluster_counts)) {
    cat(sprintf("Cluster %d: %d entries\n", cluster_counts$cluster[i], cluster_counts$n[i]))
  }
}

# Generate correction report
correction_report <- sprintf(
  "# Cluster Renaming Correction Report\n\n## Correction Strategy\nBased on gene expression patterns:\n- **Original Cluster 5 (B cell contamination)** → **Cluster 6**\n  - Key markers: Iglc3, Cd79a, Iglc2, Ebf1, Ms4a1, Cd79b, Ighm, etc.\n- **Original Cluster 6 (Proliferation)** → **Cluster 5**\n  - Key markers: H2afx, Cdca3, Hist1h3c, Mki67, Ezh2, Tuba1c, etc.\n\n## Files Updated\n1. **markers_all_clusters.csv**: All cluster markers with corrected cluster numbers\n2. **markers_top10_per_cluster.csv**: Top 10 markers per cluster with corrected cluster numbers\n3. **cluster_proportions_by_timepoint.csv**: Cluster proportions by timepoint with corrected cluster numbers\n\n## Backup Location\nAll files have been backed up to: `%s`\n\n## Next Steps\n1. Generate visualizations with hidden cluster 6 (B cell contamination)\n2. Proceed with functional enrichment analysis using corrected cluster numbering\n3. Update documentation to reflect the corrected cluster assignments\n\n---\n*Generated: %s*",
  backup_dir,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(correction_report, file.path(reports_dir, "cluster_renaming_correction_report.md"))
cat("Cluster renaming correction completed successfully!\n")
cat(sprintf("Correction report: %s\n", file.path(reports_dir, "cluster_renaming_correction_report.md")))