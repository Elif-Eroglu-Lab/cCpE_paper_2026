# README
# Author: Inika Prasad
# Aim: To merge together samples into the same Seurat object
# retaining unique columns for per-sample clustering
# and merging columns with metrics like nFeature, nCount, etc.
# This also where doublet removal would happen

# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

## Load libraries
library(Seurat)
library(dplyr)

set.seed(100)

# Evaluate Snakemake variables
filepaths <- snakemake@input[["rds"]]

mode <- snakemake@params$filtering_mode
exclude_categories <- snakemake@params$exclude_categories
metadata_field <- snakemake@params$metadata_field
metadata_fields2merge <- snakemake@params$metadata_fields2merge

output_txt <- snakemake@output[["txt"]]
output_rds <- snakemake@output[["rds"]]
metadata_backup <- snakemake@output[["metadata_backup"]]

# # Read in all the objects
seurat_list <- lapply(filepaths, readRDS)

# For each sample, all metadata should get the sample ID as prefix (exceptions in config file)
for (i in seq_along(seurat_list)){
    seu_obj <- seurat_list[[i]] # extract Seurat object from the list
    prefix <- unique(seu_obj$orig.ident)
    print(paste0("prefix is ", prefix))
    all_cols <- colnames(seu_obj@meta.data) # extract metadata columns
    rename_bool <- !(all_cols %in% metadata_fields2merge) #the columns to rename will be TRUE
    all_cols[rename_bool] <- paste0(prefix, "_", all_cols[rename_bool]) #rename the columns
    all_cols #verify that the columns are correctly renamed
    print(paste0("renamed_columns are ", paste(all_cols, collapse = ", ")))
    colnames(seu_obj@meta.data) <- all_cols
    colnames(seu_obj@meta.data)
    seurat_list[[i]] <- seu_obj # reassign modified seurat object to list
}

# # Merge all the objects
combined <- merge(seurat_list[[1]], y = seurat_list[-1])
combined$orig.ident %>% unique()

# Subset cells based on Doublet Confidence
cells2keep <- !combined[[metadata_field]][, 1] %in% exclude_categories
seurat_obj <- subset(combined, cells = colnames(combined)[cells2keep])

# Save info about the merger
sink(output_txt)
print("The files merged together can be found in:")
print(filepaths)
print("The chosen filtering mode is:")
print(mode)
print("The categories to be excluded are:")
print(exclude_categories)
print("The categories are defined in this metadata column")
print(metadata_field)
print("This metadata column has these levels:")
print(seurat_obj[[metadata_field]] %>% unique())
sink()

# Save metadata for the object before we remove a bunch of columns
write.csv(seurat_obj@meta.data, metadata_backup, row.names = TRUE)

# Remove metadata that is sample-specific to make saving easier/output object lighter
keep <- intersect(metadata_fields2merge,
                  colnames(seurat_obj@meta.data))

# drop = FALSE ensures result is a df
seurat_obj@meta.data <- seurat_obj@meta.data[, keep, drop = FALSE]

# Save the merged object
saveRDS(seurat_obj, output_rds)

# Close log file
sink(type = "message")
sink()
close(logfile)