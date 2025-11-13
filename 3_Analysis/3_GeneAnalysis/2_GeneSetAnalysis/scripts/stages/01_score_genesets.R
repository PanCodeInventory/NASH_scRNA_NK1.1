# 01_score_genesets.R

# 1. Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

# 2. Define file paths
rds_file <- "1_Files/RDS/nk.v5.rds"
gene_set_file <- "3_Analysis/6_GeneSetAnalysis/geneset.txt"
output_plot_dir <- "3_Analysis/6_GeneSetAnalysis/results/plots/"

# Create output directory if it doesn't exist
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Read and parse the gene set file
message("Parsing gene sets from: ", gene_set_file)
lines <- readLines(gene_set_file)
gene_sets <- list()
current_set_name <- NULL

for (line in lines) {
  line <- trimws(line)
  if (startsWith(line, "#")) {
    current_set_name <- sub("# ", "", line)
    gene_sets[[current_set_name]] <- c()
  } else if (nchar(line) > 0 && !is.null(current_set_name)) {
    # Extract genes, remove quotes, and split by comma and whitespace
    genes <- str_extract_all(line, '"([^"]+)"')[[1]]
    genes <- gsub('"', '', genes)
    gene_sets[[current_set_name]] <- c(gene_sets[[current_set_name]], genes)
  }
}
message("Successfully parsed ", length(gene_sets), " gene sets.")

# 4. Load the Seurat object
message("Loading Seurat object from: ", rds_file)
seurat_obj <- readRDS(rds_file)
message("Seurat object loaded successfully.")

# 5. Calculate module scores for each gene set
message("Calculating module scores...")
for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  # Ensure genes are present in the object to avoid errors
  genes_in_object <- genes[genes %in% rownames(seurat_obj)]
  
  if (length(genes_in_object) == 0) {
    warning("No genes from '", set_name, "' found in the Seurat object. Skipping.")
    next
  }
  
  # Sanitize set_name for use as a metadata column name
  feature_name <- str_replace_all(set_name, " ", "_")
  
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(genes_in_object),
    name = feature_name
  )
  message("Calculated score for '", set_name, "'")
}

# 6. Generate and save violin plots for each score
message("Generating violin plots...")
for (set_name in names(gene_sets)) {
  feature_name <- str_replace_all(set_name, " ", "_")
  # The actual feature name created by AddModuleScore is name + number (e.g., Mouse_iNK_Markers1)
  score_col <- paste0(feature_name, "1") 
  
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    warning("Score column '", score_col, "' not found in metadata. Skipping plot.")
    next
  }

  p <- VlnPlot(
    object = seurat_obj,
    features = score_col,
    group.by = "seurat_clusters", # Assuming clusters are stored here
    pt.size = 0
  ) + 
  labs(title = paste("Expression of", set_name), y = "Module Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Sanitize file name
  plot_filename <- paste0("vlnplot_", str_replace_all(set_name, "[^a-zA-Z0-9_]", "-"), ".png")
  output_path <- file.path(output_plot_dir, plot_filename)
  
  ggsave(
    filename = output_path,
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
  message("Saved plot to: ", output_path)
}

message("Analysis complete. All plots saved.")
