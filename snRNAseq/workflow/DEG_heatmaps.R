################## README
# Author: Inika Prasad
# Date created: 2025-11-10
# Aim: Read in DEG lists, fold change, make a heatmap to get overview of different DEGs
# Current usage: to visualize DEGs for epi peeout clusters
# Takes input coming from makeheatmaps4DEGs.R

################## Load libraries
library(readxl)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(scales)
library(pals)
library(stringr)

################### Rename data
all_degs <- all_sig_degs_fcvals

#################### Combine DEG lists
# all_degs <- purrr::reduce(all_degs, full_join, by = "gene")

#################### Visualize log2FC values

# Identify minimum & maximum log2 Fold Change
all_degs %>% dplyr::select(where(is.numeric)) %>% 
    pivot_longer(
        cols = everything(),
        names_to = "variable",
        values_to = "value"
    ) %>% 
ggplot(., aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue") +
  facet_wrap(~ variable, scales = "free") +
  theme_gray()

ggsave(log2FC_dist_filename)


################### Prep heatmap data
# Turn NA values to 0
all_degs[is.na(all_degs)] <- 0

all_degs

# Ensure rownames are gene names
rownames(all_degs) <- all_degs$gene
head(all_degs)

# Prep data for heatmap
mat <- all_degs %>%
  subset(., select = -gene) %>%
  as.matrix()

annot_mat <- all_sig_degs_adjpvals %>%
  subset(., select = -gene) %>%
  as.matrix()


################### Prep heatmap colors and scale

# Make a symmetrical scale encompassing the min and max values
scale_end <- max(abs(floor(min(mat, na.rm = TRUE))), abs(ceiling(max(mat, na.rm = TRUE))))
scale_end

my_palette <- colorRampPalette((rev(brewer.rdylbu(40))))(100)
show_col(my_palette, ncol = 10)

# Set the middle of the heatmap (0) to grey.
# Unnecessary if you have the p val annotations
# mid <- ceiling(length(my_palette) / 2)
# my_palette[mid] <- "grey"

# Define the legend
breaks <- seq(-scale_end, scale_end, length.out = 101) # means 100 colors

################## Define annotations
# Make annotation dataframe
col_annot_df <- data.frame(Comparison = Comparison,
                           Cluster = Cluster)
rownames(col_annot_df) <- colnames(mat)
col_annot_df

################## Make heatmap

# filename <- file.path(paste0("results/testerdump/", "STARS_Cluster04_vehccpe", ".png"))
# filename <- file.path(paste0("results/figures/", "AGAINCluster04_vehccpe_RECOMPLOG2FC", ".pdf"))
# filename <- "results/figures_megaminus/CHANGETHIS.pdf"
#

# Convert gene names to sentence case (for Pw naming convention)
rownames(mat) <- str_to_sentence(rownames(mat))

pheatmap(t(mat),
         color = my_palette,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 10,
         fontsize_col = 14,
        #  labels_row = rownames(mat2),
        #  labels_col = colnames(mat2),
        #  main = "Differentially Expressed Genes",
         angle_col = 45,
         gaps_col = mygaps,
        #  cellheight = 15,
         treeheight_row = 0,
         treeheight_col = 0,
         border_color = "white",
        #  na_col = "grey",
        display_numbers = t(annot_mat),
        fontsize = 7,
        #  annotation_col = col_annot_df,
         annotation_row = col_annot_df,
         annotation_colors = ann_colors,
         breaks = breaks,
         cellheight = 30,
         cellwidth = 22,
         legend_breaks = c(-scale_end, -scale_end/2, 0, scale_end/2, scale_end),
         drop_levels = FALSE,
        #  width = 24,
         height = 8,
         filename = filename,
         main = paste0("log2FC cutoff = ", log2FC_cutoff)
        )

dev.off()

################## Make scrunched heatmap

pdf("squished_heatmap.pdf", width = 5, height = 5)
pheatmap(mat,
         color = my_palette,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 10,
         fontsize_col = 14,
         labels_row = rownames(mat),
         show_rownames = FALSE,
         gaps_col = mygaps,
         labels_col = colnames(mat),
         angle_col = 45,
         cellwidth = 15,
         annotation_col = col_annot_df,
         annotation_colors = ann_colors,
         treeheight_row = 0,
         breaks = breaks)
dev.off()
