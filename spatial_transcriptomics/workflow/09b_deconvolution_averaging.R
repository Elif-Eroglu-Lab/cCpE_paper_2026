################################ README
# Author: Inika Prasad, inika.prasad@gmail.com
# Date: 2026-04-26
# Aim: To average prediction scores across different clusters
# for categories like Cardiomyocytes, Erythrocytes, etc.

################################ Load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(writexl)
library(readxl)
library(patchwork)
library(tibble)
library(ggpointdensity)
library(ggnewscale)

set.seed(100)

source("workflow/functions.R")

#################################### Load variables
whichcondition <- "all"
output_rds <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_seurat_object.rds")

# Specify details
# lowerlimit <- 0
lowerlimit <- 0.01
avgpred_allsections <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_averagedpreds_sections.pdf")
avgpreds_umaps <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_averagedpreds_umaps.pdf")
pngfolderhome <- paste0("../results/09_deconvolution/", whichcondition)


#################################### Evaluate unmapped clusters
# NOTE:ELEPHANT. START HERE IF YOU'RE NOT RERUNNING THE ANALYSIS
spatial <- readRDS(output_rds)

# If some of the cell types are not found their rowSums will be 0
spatial %>% rowSums() != 0

spatial@assays$predictions$data %>% rowSums()

######### Group together scores for some clusters

# Specify which clusters to group together
map <- c(
  "01" = "CM",
  "05" = "CM",
  "15" = "CM",
  "20" = "CM",
  "25" = "CM",
  "18" = "CM",
  "07" = "Erythrocytes",
  "14" = "Erythrocytes",
  "10" = "Myeloid_cells",
  "06" = "Myeloid_cells",
  "17" = "Lymphocytes",
  "23" = "Lymphocytes",
  # "16" = "Epi",
  # "26" = "Epi",
  # "09" = "Epi",
  "21" = "SMCs",
  "19" = "Valve",
  "12" = "EC_FB",
  "11" = "EC_FB",
  "22" = "EC_FB",
  "03" = "EC_FB",
  "02" = "EC_FB",
  "04" = "EC_FB",
  "08" = "EC_FB",
  "13" = "EC_FB",
  "24" = "EC_FB",
  "27" = "Glial",
  "activated-epi" = "Epi",
  "Epi-03" = "Epi",
  "EPDC-A" = "Epi",
  "EPDC-02" = "Epi",
  "Epi-07" = "Epi",
  "Epi-musclesig" = "Epi",
  "EPDC-boomerang" = "Epi",
  "EPDC-12" = "Epi",
  "EPDC-11" = "Epi",
  "EPDC-14" = "Epi",
  "EPDC-09" = "Epi",
  "EPDC-16" = "Epi",
  "EPDC-15" = "Epi",
  "max" = "max_vestigial"
)

# Extract prediction scores
avg_preds <- spatial@assays$predictions$data

# Average prediction scores
avg_preds <- avg_preds %>%
  as.data.frame() %>%
  rownames_to_column("rowname") %>%
  mutate(group = map[rowname]) %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("group")

# Set all values less than or equal to the lowerlimit as NA
lowerlimit # check the lower limit
avg_preds[avg_preds <= lowerlimit] <- NA

rownames(avg_preds)

spatial <- AddMetaData(spatial, metadata = avg_preds %>% as.matrix() %>% t())

########################### Make figure: Prediction scores on a UMAP
FeaturePlot(spatial, features = unique(map)) &
    scale_color_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    )

ggsave(avgpreds_umaps, height = 20, width = 25)


############################ View Prediction scores on sections
source("workflow/functions.R")
x <- plotall(spatial,
  crop = TRUE,
  # selected_clustering = myclusters,
  # selected_cols = selected_cols,
  features2plot = unique(map),
  image.alpha = 0.3,
  pt.size.factor = 3.2,
  # cells.highlight = cells
)

x$naive &
    scale_fill_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    )
x$vehicle.5 &
    scale_fill_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    )
x$c.CpE.6 &
    scale_fill_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    )


################  Make pdf with all sections, overview
pdf(avgpred_allsections)
for (i in seq_along((x))) {

  print(x[[i]] &
    scale_fill_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    ))
}

dev.off()

################## Make pngs with the sections of interest
# Specify the section for which you want plots
image <- "vehicle.5"
image <- "vehicle.7"
image <- "c.CpE.6"
image <- "naive.1"

# Print a plot for every predicted class
for (i in unique(map)) {

  # Make a filepath for the image
  filename <- file.path(pngfolderhome,
                        image,
                        paste0(image, "_", i, "_NA", lowerlimit, "AVGPRED.png"))
  filename

  # Make a filepath for the image
  filename_legend <- file.path(pngfolderhome,
                               image,
                               paste0(image, "_", i, "_NA", lowerlimit, "_legend", "AVGPRED.png"))
  filename_legend

  # Create the directory to store images if it doesn't exist already
  if (!dir.exists(dirname(filename))) {
    dir.create(dirname(filename))
  }

  x <- SpatialPlot(
    object = spatial,
    images = image,
    features = i,
    image.alpha = 0.5,
    pt.size.factor = 3.5,
    stroke = 0.5
  ) &
    scale_fill_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    ) &
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      # make legend bigger
      legend.key.height = unit(3, "cm"),
      legend.key.width  = unit(0.6, "cm"),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 30)
    )

  # Extract plot and legend separately
  plot_no_legend <- x +
    theme(legend.position = "none") +
    ggtitle(paste0(i, ", NA threshold: ", lowerlimit))
  legend <- cowplot::get_legend(x)

  # Save the plot and legend
  ggsave(filename, plot_no_legend)
  ggsave(filename_legend, legend)


}
