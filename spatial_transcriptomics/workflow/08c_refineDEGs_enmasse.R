# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aims:

# Since there is a batch effect AND
# differences in degree of tissue regenerated in vehicle batch 1
# vs. vehicle batch2,
# we did a veh vs. ccpe, and a vehbatch1 vs. vehbatch2
# and now we will select the DEGs that are ONLY present in
# veh vs. ccpe.
# These are the refined DEGs.

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

################################ Load Libraries
library(dplyr)
library(writexl)
library(readxl)

############################ README
# Author: Inika Prasad
# Aim: We want to take the intersection of differentially expressed genes in
# 10dpci ccpe vs 10dpci veh (genes disrupted by epi TJ modulation)
# 10dpci veh vs 14dpsham (genes important for regeneration)

# to identy genes that are important for regeneration AND
# messed up by TJ modulation

############################ Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(writexl)

set.seed(100)

############################ Specify the unique codes
code1 <- "ZANIA"
code2 <- "ZAURAK"

# Unique codes
code1 <- snakemake@params[["code1"]]
code2 <- snakemake@params[["code2"]]

compname <- paste0(code1, "_", code2)
print(compname)

# Read in the two folders
folder1 <- paste0("../results/08_DGE/", code1)
folder2  <- paste0("../results/08_DGE/", code2)

# intersection_file <- paste0("../results/DGE/intersections_", compname, ".xlsx")

intersection_file <- snakemake@output[["refdegs"]]

print("Folders to be used")
print(folder1)
print(folder2)

############################ Identify pairs of files

# List Excel files
files1 <- list.files(folder1, pattern = "\\.xlsx$", full.names = TRUE)
files2 <- list.files(folder2, pattern = "\\.xlsx$", full.names = TRUE)

print("Available files inside folders")
print(files1)
print(files2)

# Extract cluster names from filenames
clust1 <- sub(paste0(".*", code1, "(.*)", "_DGE_results.xlsx", "$"), "\\1", files1)
clust2 <- sub(paste0(".*", code2, "(.*)", "_DGE_results.xlsx", "$"), "\\1", files2)

# Create data frames with filenames and cluster numbers
df1 <- data.frame(file1 = files1, clust = clust1, stringsAsFactors = FALSE)
df2 <- data.frame(file2 = files2, clust = clust2, stringsAsFactors = FALSE)

# Match file pairs by number
matched_pairs <- merge(df1, df2, by = "clust")

# Check the matched pairs
print(matched_pairs)

# Remove any unmatched rows, and the summary excel
matched_pairs <- matched_pairs %>% filter(!is.na(clust) & clust != "/summary" & clust != "")

############################ Intersection per pair

# Empty list to store intersection genes
filtered_deg_dfs <- list()

# Identify intersection genes
for (i in seq_len(nrow(matched_pairs))) {

  # Define the filenames of interest
  degfile1 <- matched_pairs[i, "file1"]
  degfile2 <- matched_pairs[i, "file2"]
  clustnum <- matched_pairs[i, "clust"]

  # Name of the comparison:
  sheetname <- paste0(compname, "_", clustnum)

  # Read in the pairs and find the intersections
  degs1 <- read_xlsx(degfile1, sheet = "sig_degs")
  degs2 <- read_xlsx(degfile2, sheet = "sig_degs")

  # Find genes only changing in vehccpe (only in degs1)
  refined_genes <- setdiff(degs1$gene, degs2$gene)

  # if there are refined_genes, save their stats
  if (length(refined_genes) > 0) {
    filt_degs1 <- degs1 %>% filter(gene %in% refined_genes)
    # Save the joined dataframe of common genes
    filtered_deg_dfs[[sheetname]] <- filt_degs1

    # QC: this should be 0!
    filt_degs2 <- degs2 %>% filter(gene %in% refined_genes)
    print(filt_degs2)
  }

}

# Save the intersection of the two dataset DEGs
write_xlsx(filtered_deg_dfs, intersection_file)
