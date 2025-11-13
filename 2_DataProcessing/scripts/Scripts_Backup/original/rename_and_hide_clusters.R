#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Paths
rds_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds"
out_rds_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds"
reports_dir <- "2_DataProcessing/reports"
backup_dir <- file.path(reports_dir, "backup")

# Create directories
dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)

# Logging
log_file <- file.path(reports_dir, sprintf("cluster_renaming_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
con <- file(log_file, open = "wt")
wlog <- function(...) writeLines(paste0(...), con = con)

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

tryCatch({
  wlog("[INFO] Starting cluster renaming process")
  wlog(sprintf("[INFO] Loading RDS: %s", rds_path))
  
  # Load original object
  obj <- readRDS(rds_path)
  wlog(sprintf("[INFO] Loaded object with %d cells", ncol(obj)))
  
  # Backup original cluster assignments
  original_clusters <- obj@meta.data$seurat_clusters
  if (is.null(original_clusters)) {
    original_clusters <- Idents(obj)
  }
  wlog("[INFO] Backed up original cluster assignments")
  
  # Apply cluster renaming
  wlog("[INFO] Applying cluster renaming mapping")
  renamed_clusters <- as.character(original_clusters)
  
  # Apply mapping
  for (old_name in names(cluster_mapping)) {
    new_name <- cluster_mapping[old_name]
    renamed_clusters[renamed_clusters == old_name] <- as.character(new_name)
    wlog(sprintf("[INFO] Cluster %s → %s", old_name, new_name))
  }
  
  # Update Seurat object
  obj@meta.data$seurat_clusters <- factor(renamed_clusters)
  Idents(obj) <- factor(renamed_clusters)
  
  # Add original cluster as metadata for reference
  obj@meta.data$original_cluster <- original_clusters
  
  wlog("[INFO] Updated Seurat object with new cluster assignments")
  
  # Save renamed object
  saveRDS(obj, out_rds_path)
  wlog(sprintf("[INFO] Saved renamed RDS: %s", out_rds_path))
  
  # Generate cluster summary
  cluster_summary <- table(
    original = as.character(original_clusters),
    renamed = renamed_clusters
  )
  wlog("[INFO] Generated cluster summary:")
  for (i in 1:nrow(cluster_summary)) {
    wlog(sprintf("[INFO]   Original %s → Renamed %s: %d cells", 
                 rownames(cluster_summary)[i], 
                 colnames(cluster_summary)[i], 
                 cluster_summary[i]))
  }
  
  # Save cluster summary
  summary_df <- as.data.frame(cluster_summary)
  summary_df <- summary_df[summary_df$Freq > 0, ]
  write_csv(summary_df, file.path(reports_dir, "cluster_renaming_summary.csv"))
  wlog(sprintf("[INFO] Saved cluster summary: %s", file.path(reports_dir, "cluster_renaming_summary.csv")))
  
  # Create renaming report
  summary_table <- paste(sprintf("| %s | %s | %d |",
                                rownames(cluster_summary),
                                colnames(cluster_summary),
                                cluster_summary), collapse = "\n")
  
  report_content <- sprintf(
    "# Cluster Renaming Report\n\n## Renaming Strategy\n- **Cluster 5 → Cluster 6**: B cell contamination cluster (Iglc3, Cd79a, Iglc2, Ebf1, Ms4a1, etc.)\n- **Cluster 6 → Cluster 5**: Proliferation cluster (H2afx, Cdca3, Hist1h3c, Mki67, etc.)\n\n## Cluster Summary\n\n| Original Cluster | Renamed Cluster | Cell Count |\n|------------------|-----------------|------------|\n%s\n\n## Files Generated\n- **Renamed RDS**: `%s`\n- **Original RDS**: `%s`\n- **Cluster Summary**: `cluster_renaming_summary.csv`\n\n## Next Steps\n1. Update analysis result files (markers, proportions)\n2. Generate visualizations with hidden cluster 6\n3. Validate downstream analysis compatibility\n\n---\n*Generated: %s*",
    summary_table,
    out_rds_path,
    rds_path,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  writeLines(report_content, file.path(reports_dir, "cluster_renaming_report.md"))
  wlog(sprintf("[INFO] Created renaming report: %s", file.path(reports_dir, "cluster_renaming_report.md")))
  
  wlog("[INFO] Cluster renaming completed successfully")
  
}, error = function(e) {
  wlog(paste0("[ERROR] ", e$message))
  close(con)
  quit(status = 1)
})

close(con)
cat("Cluster renaming completed successfully!\n")
cat(sprintf("Check the report: %s\n", file.path(reports_dir, "cluster_renaming_report.md")))