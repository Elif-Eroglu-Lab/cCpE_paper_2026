############################## README
# Author: Inika Prasad
# Aim: To merge together spatial transcriptomics samples
# into the same Seurat object

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

###############################  Load libraries
library(Seurat)
library(dplyr)
library(readxl)

source("workflow/functions.R")

set.seed(100)
############################### Load Snakemake variables
metadata_file <- snakemake@input[["metadata"]]
output_rds <- snakemake@output[["rds"]]

############################### Tester variables
# metadata_file <- "../data/metadata.xlsx"
# output_rds <- "../results/merged.rds"

############################### Read in metadata
metadata <- read_xlsx(metadata_file)

# Add to metadata
metadata$orig.ident <- metadata$sample # create a column called orig.ident
metadata$slidenum <- sub("-[^-]+$", "", metadata$slide) #remove A1, B1 suffixes

################################ Read in all the objects
# Identify the filepaths for objects to be read in
filepaths <- metadata$path
names(filepaths) <- metadata$sample

filepaths

# Read in all the objects
seurat_list <- mapply(read_spatial,
                      path  = filepaths,
                      label = names(filepaths),
                      slice = names(filepaths),
                      SIMPLIFY = FALSE)

# Merge the objects
obj <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

################################ Add metadata to object
# Extract Seurat object's metadata to populate it
metadata2 <- obj[[]]
metadata2$spot <- rownames(metadata2)
head(metadata2)

# Combine existing and new metadata
metadata3 <- merge(metadata2, metadata, by = "orig.ident", all.x = TRUE)
rownames(metadata3) <- metadata3$spot

head(metadata3)

# Add metadata to Seurat object
obj <- AddMetaData(obj, metadata3)

# Make sure variables are factors as needed
obj$path <- NULL
obj$condition <- factor(obj$condition, levels = c("naive", "vehicle", "c-CpE"))
obj$replicate <- factor(obj$replicate)
obj$batch <- factor(obj$batch, levels = c("1", "2"))
obj$animal_readable <- factor(obj$animal_readable)
obj$slide <- factor(obj$slidenum)

# Make sample order more intuitive
obj$sample <- factor(obj$sample,
                     levels = c("naive-1", "vehicle-1", "vehicle-2", "vehicle-3",
                                "vehicle-4", "vehicle-5", "vehicle-6",
                                "vehicle-7", "vehicle-8",
                                "c-CpE-1", "c-CpE-2", "c-CpE-3",
                                "c-CpE-4", "c-CpE-5", "c-CpE-6"))

obj$orig.ident <- factor(obj$orig.ident,
                         levels = c("naive-1", "vehicle-1", "vehicle-2", "vehicle-3",
                                    "vehicle-4", "vehicle-5", "vehicle-6",
                                    "vehicle-7", "vehicle-8",
                                    "c-CpE-1", "c-CpE-2", "c-CpE-3",
                                    "c-CpE-4", "c-CpE-5", "c-CpE-6"))


levels(obj$orig.ident)
Idents(obj) <- "orig.ident"

################################ Make count layers nicely renamed
# Joining & resplitting gives the count layers sensible names instead of counts.1, counts.2
obj <- JoinLayers(obj)
obj[['Spatial']] <-  split(x = obj[['Spatial']], f = obj$sample)

# Check how the naming went
obj

################################ Rotate the sections
# Aim: to rotate the tissue section anticlockwise 90 degrees
# for the images that have the injury site facing left

# All the coordinates have positive x and y values
# i.e. they're in quadrant I
# (i.e. the origin is on the bottom left)
print("rotating sections now...")
obj$x %>% range()
obj$y %>% range()

# Which samples is this going to apply to?
obj$sample %>% unique()
samples2rotate <- c("naive-1", "vehicle-1", "vehicle-2", "vehicle-3", "vehicle-4",
                    "c-CpE-1", "c-CpE-2", "c-CpE-3", "c-CpE-4", "c-CpE-5", "c-CpE-6")

# Rotate and transpose the coordinates so that they return to quadrant I
subset_metadata <- obj@meta.data %>%
  filter(sample %in% samples2rotate) %>%
  mutate(
    spun_x = -y,
    spun_y = x
  ) %>%
  mutate(
    rot_x = spun_x - min(spun_x),
    rot_y = spun_y - min(spun_y)
  )

# Rename columns so that rotated coordinates are in x and y columns
renaming <- c(x = "rot_x", y = "rot_y", original_x = "x", original_y = "y")
subset_metadata_rot <- rename(subset_metadata, all_of(renaming))

# Replace old metadata for samples2rotate
obj@meta.data[rownames(subset_metadata_rot), c("x", "y")] <- subset_metadata_rot[, c("x", "y")]


################################## Rotate images themselves
for (sample in samples2rotate) {
  sample <- gsub("-", ".", sample)
  if (sample %in% names(obj@images)) # make sure the section name is present
  obj <- rotate_visiumv2(obj, slice = sample)

}

################################# Save merged object
# Save the merged object
saveRDS(obj, output_rds)

# Output which samples were in the object
print("Samples merged in this object are:")
obj$orig.ident %>% unique()

################################# Close log file
sink(type = "message")
sink()
close(logfile)