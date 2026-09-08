############################ README
# Author: Inika Prasad
# Aim: We want to make heatmaps for DEGs from
# 10dpci ccpe vs 10dpci veh (genes disrupted by epi TJ modulation)
# 10dpci veh vs 14dpsham (genes important for regeneration)
# 14dpsham vs PMnaive (genes disrupted by non-injurious surgery, 
# i.e. generic epi activation)

############################ Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(writexl)
library(RColorBrewer)
library(Seurat)

set.seed(100)

############################ Specify the unique codes
# Unique codes
code1 <- "ALGEDImod"
code2 <- "ALGENIB"
code3 <- "ALKAID"

# Unique codes for Epi peelout
code1 <- "ALKALUROPS" # naive-sham
code2 <- "ALGORAB" # sham-veh
code3 <- "ALGOL" # veh-ccpe

# WE WILL REMOVE GENES THAT ARE NOT DE IN CODE3
# SO THE ORDER MATTERS
# See line 121 for exact usage

############################ Identify pairs of files
# Read in the two folders
folder1 <- paste0("results/DGE/", code1)
folder2  <- paste0("results/DGE/", code2)
folder3 <- paste0("results/DGE/", code3)

compname <- paste0(code1, "_", code2, "_", code3)
intersection_file <- paste0("results/DGE/intersections_", compname, ".xlsx")

# List Excel files
files1 <- list.files(folder1, pattern = "\\.xlsx$", full.names = TRUE)
files2 <- list.files(folder2, pattern = "\\.xlsx$", full.names = TRUE)
files3 <- list.files(folder3, pattern = "\\.xlsx$", full.names = TRUE)


# Extract cluster names from filenames
clust1 <- sub(paste0(".*", code1, "(.*)", "_DGE_results.xlsx", "$"), "\\1", files1)
clust2 <- sub(paste0(".*", code2, "(.*)", "_DGE_results.xlsx", "$"), "\\1", files2)
clust3 <- sub(paste0(".*", code3, "(.*)", "_DGE_results.xlsx", "$"), "\\1", files3)


# Create data frames with filenames and cluster numbers
df1 <- data.frame(file1 = files1, clust = clust1, stringsAsFactors = FALSE)
df2 <- data.frame(file2 = files2, clust = clust2, stringsAsFactors = FALSE)
df3 <- data.frame(file3 = files3, clust = clust3, stringsAsFactors = FALSE)


# Match file pairs by number
# matched_pairs <- merge(df1, df2, by = "clust")
matched_pairs <- purrr::reduce(.x = list(df1, df2, df3), merge, by = c("clust"), all = TRUE)

# Check the matched pairs
matched_pairs %>% View()
nrow(matched_pairs)

# Remove any unmatched rows, and the summary excel
matched_pairs <- matched_pairs %>% filter(!is.na(clust) & clust != "/summary")
nrow(matched_pairs)


############################ Prep data for heatmap

# Empty list to store combined dataframe
allsigdegs <- list()

# Extract and combine dataframe
for (i in seq_len(nrow(matched_pairs))) {

  # Define the filenames of interest
  degfile1 <- matched_pairs[i, "file1"]
  degfile2 <- matched_pairs[i, "file2"]
  degfile3 <- matched_pairs[i, "file3"]

  clustnum <- matched_pairs[i, "clust"]

  # Name of the comparison:
  sheetname <- paste0(compname, "_", clustnum)

  # Read in dataframes from excel sheets, remove any empty dfs
  file_list <- list(degfile1, degfile2, degfile3)
  names(file_list) <- c(code1, code2, code3)
  file_list <- file_list[!is.na(file_list)]

  df_list <- lapply(names(file_list), function(degname) {
    file <- file_list[[degname]]

    read_xlsx(file, sheet = "sig_degs") %>%
      # as.data.frame() %>%
      select(gene, avg_log2FC) %>%
      rename_with(~ degname, .cols = -gene)
  })

  # Combine dataframes for the heatmap
  allsigdegs[[sheetname]] <- purrr::reduce(.x = df_list, merge, by = c("gene"), all = TRUE)

}

allsigdegs %>% head()


#################### For each element of the list(i.e. each cluster), make a heatmap

for (i in seq_along(allsigdegs)) {

  label4heatmap <- names(allsigdegs)[i]

  all_degs <- allsigdegs[[i]] %>% as.data.frame()
  print(label4heatmap)
  print(nrow(all_degs))

  # Remove the genes which don't change in veh-ccpe (i.e. code3)
  rownums2keep <- !is.na(all_degs[, code3])
  all_degs <- all_degs[rownums2keep, ]
  labelmod <- "veh-ccpe-only"

  # # Remove genes which aren't ERK targets
  # source("workflow/ERKgenes_heatmap.R")
  # all_degs <- all_degs %>% filter(gene %in% ourERKgenes_pro)
  # labelmod <- "ERKtargets"

  # Turn NA values to 0
  all_degs[is.na(all_degs)] <- 0

  # Ensure rownames are gene names
  rownames(all_degs) <- all_degs$gene
  head(all_degs)

  # Make color palette
  my_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

  # Prep data for heatmap
  library(pheatmap)
  mat <- all_degs %>%
    subset(., select = -gene) %>%
    MinMax(data = ., min = -2.5, max = 2.5) %>%
    as.matrix()

  breaks <- seq(-2.5, 2.5, length.out = 101)

  print(nrow(all_degs))

  if (nrow(all_degs) > 0) {

    pdf(file.path(paste0("results/testerdump/", label4heatmap, labelmod, ".pdf")),
        height = 0.12 * nrow(all_degs) + 4,
        width = 7,
        bg = "transparent")

    pheatmap(mat,
             color = my_palette,
             cluster_rows = nrow(all_degs) > 2,
             cluster_cols = FALSE,
             fontsize_row = 10,
             fontsize_col = 14,
             labels_row = rownames(mat),
             labels_col = colnames(mat),
             main = label4heatmap,
             angle_col = 45,
             cellwidth = 25,
             treeheight_row = 0,
             breaks = breaks)

    dev.off()

  }

}

######################### Further investigation
intersect(rownames(all_degs), read_xlsx("results/DGE/ALKES/ALKES_DGE_results.xlsx", sheet = "sig_degs") %>% pull(gene))
################################## Violin Plots
VlnPlot(obj3, features = c("NR4A1", "NR4A2", "NR4A3", "STYX"), group.by = "orig.ident")


obj3 <- subset(obj2, named_clusters == c("Epi1", "Activated epi"))
obj3$cluster_condition <- paste0(obj3@meta.data[[clustering]], "_", obj3@meta.data[["orig.ident"]])
obj3$cluster_condition %>% unique()

obj3 <- subset(obj3, cluster_condition == c("Epi1_PM_naive", "Activated epi_14dpsham",
                                            "Activated epi_10dpci_veh", "Activated epi_10dpci_cCpE") %>% rev())

obj3$cluster_condition <- factor(obj3$cluster_condition,
                                 levels = c("Epi1_PM_naive", "Activated epi_14dpsham",
                                            "Activated epi_10dpci_veh", "Activated epi_10dpci_cCpE"))

obj3$cluster_condition %>% unique()

DotPlot(
  obj3,
  scale = FALSE,
  features = rownames(all_degs),
  group.by = "cluster_condition") +
  scale_colour_viridis(option = "viridis") +
  RotatedAxis() +
  coord_flip()

ggsave("results/testerdump/cluster04DEGs_dotplot.png", width = 6)

################## Check which DEGs are in the Top markers for Epi peelout clusters
markerlist <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/custom_res1_markers.xlsx"

# Load marker list
m <- read_xlsx(markerlist)
m %>% group_by(gene) %>% reframe(unique(cluster)) %>% stack()

# Empty DEGs dataframe turned to NA
allsigdegs <- lapply(allsigdegs, function(clusterdegs){
  if (nrow(clusterdegs) > 0) {
    clusterdegs
  } else {
    (NA)
  }
})

# Remove empty DEG dataframe
allsigdegs <- allsigdegs[!is.na(allsigdegs)]
combined <- bind_rows(allsigdegs, .id = "source")
View(combined)

# Save the combined dataframe
write_xlsx(combined, "results/modules/DEgenesEpiPeelout.xlsx")

# How many DEGs are also markers?
intersect(unique(combined$gene), unique(m$gene)) %>% length()

# What are the veh-ccpe DEGs?
vehccpe_degs <- combined %>% filter(!is.na(ALGOL)) %>% pull(gene) %>% unique()

# How many veh-ccpe DEGs are also markers?
intersect(vehccpe_degs, unique(m$gene)) %>% length()

# More information about DEGs and markers
for (goi in vehccpe_degs) {
  if(goi %in% m$gene) {print(goi)}
  # combined %>% filter(gene == goi) %>% print()
  # m %>% filter(gene == goi) %>% print()

}