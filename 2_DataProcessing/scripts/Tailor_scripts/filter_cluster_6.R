# filter_cluster_6.R
# Purpose: Remove all cells from cluster 6 in a Seurat object and save a robust, readable RDS
# Input:  1_Files/RDS/nk.integrated.v4.rds
# Output: 1_Files/RDS/nk.v5.rds (atomic write + verification)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
})

# ------------------ Paths ------------------
input_rds_path  <- "1_Files/RDS/nk.integrated.v4.rds"
output_rds_path <- "1_Files/RDS/nk.v5.rds"

message("[1/6] Loading Seurat object from ", input_rds_path, " ...")
nk_integrated <- tryCatch({
  readRDS(input_rds_path)
}, error = function(e) {
  stop("Failed to read input RDS: ", conditionMessage(e))
})

# Update object for compatibility (safe no-op if already updated)
message("[2/6] Updating Seurat object (if needed)...")
nk_integrated <- tryCatch({
  UpdateSeuratObject(nk_integrated)
}, error = function(e) {
  warning("UpdateSeuratObject failed: ", conditionMessage(e), ". Proceeding with original object.")
  nk_integrated
})

# ------------------ Filtering logic ------------------
cluster_column <- "seurat_clusters"
if (!cluster_column %in% colnames(nk_integrated@meta.data)) {
  stop("Cluster column '", cluster_column, "' not found in metadata.")
}

message("[3/6] Original cluster distribution:")
print(table(nk_integrated@meta.data[[cluster_column]], useNA = "ifany"))

# Ensure we compare as character to avoid factor/numeric pitfalls
clusters_vec <- as.character(nk_integrated@meta.data[[cluster_column]])
keep_cells <- rownames(nk_integrated@meta.data)[!is.na(clusters_vec) & clusters_vec != "6"]

message("[4/6] Subsetting to remove cluster '6' ...")
nk_subset <- subset(nk_integrated, cells = keep_cells)

message("[5/6] New cluster distribution after filtering:")
print(table(nk_subset@meta.data[[cluster_column]], useNA = "ifany"))
message("Dimensions (features x cells): ", paste(dim(nk_subset), collapse = " x "))

# ------------------ Robust save (atomic write + verification) ------------------
# Write to a temporary file first, verify readability, then rename atomically

.tmp_path <- paste0(output_rds_path, ".tmp")
if (file.exists(.tmp_path)) {
  unlink(.tmp_path)
}

message("[6/6] Saving to temporary file: ", .tmp_path)
# Use RDS v3 for better compatibility; xz gives good compression and integrity
save_ok <- tryCatch({
  saveRDS(nk_subset, file = .tmp_path, compress = "xz", version = 3)
  TRUE
}, error = function(e) {
  message("Save failed: ", conditionMessage(e))
  FALSE
})

if (!save_ok || !file.exists(.tmp_path)) {
  stop("Temporary save failed; aborting to avoid producing a corrupted file.")
}

# Verify by reloading
message("Verifying temporary file readability...")
verify_ok <- tryCatch({
  reloaded <- readRDS(.tmp_path)
  if (!inherits(reloaded, "Seurat")) {
    stop("Reloaded object is not a Seurat object.")
  }
  message("Verification OK. Reloaded dims: ", paste(dim(reloaded), collapse = " x "))
  TRUE
}, error = function(e) {
  message("Verification failed: ", conditionMessage(e))
  FALSE
})

if (!verify_ok) {
  unlink(.tmp_path)
  stop("Verification failed; temporary file removed. No output produced.")
}

# Atomic replace
if (file.exists(output_rds_path)) {
  unlink(output_rds_path)
}
file.rename(.tmp_path, output_rds_path)

message("Done. Final RDS written to: ", output_rds_path)
