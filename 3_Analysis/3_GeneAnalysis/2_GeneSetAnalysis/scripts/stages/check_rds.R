# check_rds.R
message("Attempting to read RDS file: 1_Files/RDS/nk.v5.rds")
tryCatch({
  # Attempt to read the file
  obj <- readRDS("1_Files/RDS/nk.v5.rds")
  
  # If successful, print basic info
  message("Successfully read the RDS file.")
  message("Object class: ", class(obj))
  
  # If it's a Seurat object, print more details
  if (inherits(obj, "Seurat")) {
    message("Object is a Seurat object.")
    message("Dimensions: ", paste(dim(obj), collapse = " x "))
  }
  
}, error = function(e) {
  # If an error occurs, print the error message
  message("Failed to read the RDS file.")
  message("Error message: ", e$message)
})
