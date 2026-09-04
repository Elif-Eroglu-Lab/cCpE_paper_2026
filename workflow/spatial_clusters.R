############################## README
# Author: Inika Prasad
# Aim: Plot clusters spatially

############################## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(clustree)
library(gridExtra)
library(grid)

set.seed(100)

source("workflow/functions.R")

############################# Load data
inputfile <- "../results/06_dimred_clustering/allsections/allsections_dimred.rds"
inputfile <- "../results/01_merge_samples/01_allsections_merged.rds"


obj <- readRDS(inputfile)

# To check if coordinate flipping worked
plot_spatial_single <- function(obj,v,label=NULL,prefix=NULL,path=NULL,filepath=NULL,slot="counts",assay="Spatial",cols=NULL,lnrow=NULL,point_size=0.4,height=10,width=20){
  if (length(v) > 1) stop("Argument 'v' must be of length 1.")

  # check if variable exists
  chk <- match(v, c(rownames(obj), colnames(obj[[]])))

  chkg <- FALSE
  # check if variable to color exists
  if (!is.na(match(v,rownames(obj))) & is.na(match(v,colnames(obj[[]])))) chkg <- TRUE

  if (!is.na(chk)) {
    dfr <- FetchData(obj, vars=c("sample","x","y","condition","replicate",v), slot=slot, assay=assay)
    #group_by(id) %>%
    #mutate(mean=mean(!!sym(v))) %>%
    #ungroup() %>%
    #arrange(!!sym(v))

    p <- ggplot(dfr, aes(x=y,y=x,col=!!sym(v),fill=!!sym(v))) +
      geom_point(size = point_size, shape = 20) +
      #ggpp::geom_text_npc(aes(npcx=0.5,npcy=0,label=mean))+
      scale_y_reverse() +
      facet_grid(condition ~ replicate) +
    #   theme_simple(base_size = 8) +
      coord_fixed() +
      theme(plot.background = element_rect(fill = "white", color = NA))

    if(is.numeric(dplyr::pull(dfr, v))) {
      p <- p +
        guides(color = guide_colorbar(barwidth = 20, barheight = 0.5, title = v)) +
        guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5, title = v)) +
        scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
        scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
    }else {
      if(is.null(cols)) cols <- randomcoloR::distinctColorPalette(length(unique(dplyr::pull(dfr,v))))
      p <- p +
        guides(color = guide_legend(title=v,nrow=lnrow,override.aes=list(size=2))) +
        guides(fill = guide_legend(title=v,nrow=lnrow,override.aes=list(size=2))) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols)
    }

    if (is.null(filepath) && (is.null(path))) {
      return(p)
    } else {
      if (is.null(filepath)) {
        if (is.null(label)) label <- tolower(v)
        fname <- paste0(paste(label, "spatial", sep = "-"), ".png")
        if (!is.null(prefix)) fname <- paste0(paste(prefix, "spatial", label, sep = "-"), ".png")
        if (chkg) fname <- paste0(paste("spatial", assay, slot, label, sep = "-"), ".png")
        if (!is.null(prefix) && chkg) fname <- paste0(paste(prefix, "spatial", assay, slot, label, sep = "-"), ".png")
        if (!is.null(path)) filepath <- file.path(path, tolower(fname))

      }
      ggsave(filepath, p, height = height, width = width)
      message(paste(filepath, "exported."))

    }

  }else{
    message(paste0("Variable '", v, "' not found."))
  }
}

for (i in)

rotate_visiumv2(obj, slice = c("vehicle.1", "vehicle.2"))

plot_spatial_single(obj, v = "nCount_Spatial", label = "UMIs")

resolution <- "cca_resolution_1.8"

############################# Load data

library(patchwork)

Idents(obj) <- resolution

# High-contrast color palette
selected_cols <- pals::glasbey(30)
names(selected_cols) <- Idents(obj) %>% levels()

# Desired order of images
new_order <- gsub("-", ".", levels(obj$orig.ident))
obj@images <- obj@images[new_order]

# Plot all images with clustering info
x <- plotall(obj,
            #  selected_clustering = resolution,
             pt.size = 3,
             selected_cols = selected_cols)

x$vehicle.1
x$vehicle.2
x$vehicle.8
x$naive.1
# Plot many images per page

combined <- wrap_plots(x, ncol = 3, guides = "collect")


# Save image
pdf("tester.pdf", width = 8.3, height = 11.7)
combined
dev.off()

spatial_clusters <- list()
spatial_clusters[[1]] <- combined
spatial_clusters[[2]] <- combined

# Save image
pdf("tester.pdf", width = 8.3, height = 11.7)
spatial_clusters
dev.off()

