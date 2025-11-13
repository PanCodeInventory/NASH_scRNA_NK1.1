# check_rds.R
# Script to check if an RDS file is readable

# Define file path
rds_path <- "1_Files/RDS/nk.v5.rds"

# Try to read the RDS file
tryCatch({
    seurat_object <- readRDS(rds_path)
    print("Successfully read the RDS file.")
    print("Object summary:")
    print(seurat_object)
}, error = function(e) {
    print("Failed to read the RDS file.")
    print("Error message:")
    print(e$message)
})
