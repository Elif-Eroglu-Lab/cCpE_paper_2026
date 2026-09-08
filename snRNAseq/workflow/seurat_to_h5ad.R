
################################ Load Libraries
library(Seurat)
library(dplyr)
# library(tidyr)
# library(tibble)
# library(ggplot2)
# library(RColorBrewer)
# library(cowplot)
# library(writexl)
library(readxl)
# library(patchwork)
# library(viridis)
# library(pals)
library(scales)
# library(pheatmap)


set.seed(100)

# source("workflow/functions.R")

################################ Load data
{ # Epi peelout
input_rds <-"results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"
reduction <- "umap_peeled.harmony"
datasetname <- "10dpci_ctrl_Epipeelout"
clustering <- "merged_clusters"
# clustering <- "harmony_peeled_resolution_1.4"
ptsize <- 2
clustering_infofile <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/10dpci_ctrl_EpiPeelout_CUSTOMres.xlsx"

df <- read_xlsx(clustering_infofile, sheet = clustering)

map <- df$Cluster
names(map) <- df$Cluster

cols <- df$cols
names(cols) <- df$Cluster

cols %>% show_col()
}

############################## Load data
obj <- readRDS(input_rds)
obj2 <- JoinLayers(obj)

# make new metadata called named_clusters
obj2@meta.data$named_clusters <- factor(
  map[as.character(obj2@meta.data[[clustering]])],
  levels = unname(map))
clustering <- "named_clusters"

Idents(obj2) <- clustering
levels(obj2)

show_col(cols)

##############################
# Filtering done
# Normalization done
# Log1p done
# Highly variable genes done
# PCA done
# Neighbors done

library(Seurat)
library(Matrix)

subset_and_export <- function(seurat_obj, condition, clusters, outdir) {

  sub <- subset(
    seurat_obj,
    subset = orig.ident == condition & merged_clusters %in% clusters
  )

  sub[["RNA"]] <- JoinLayers(sub[["RNA"]])  # merge split v5 layers

  dir.create(outdir, showWarnings = FALSE)

  # counts and data as sparse matrices, genes x cells (Seurat default)
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  data   <- GetAssayData(sub, assay = "RNA", layer = "data")

  writeMM(counts, file.path(outdir, "counts.mtx"))
  writeMM(data,   file.path(outdir, "data.mtx"))

  write.csv(data.frame(gene = rownames(sub)),
            file.path(outdir, "genes.csv"), row.names = FALSE)
  write.csv(data.frame(barcode = colnames(sub)),
            file.path(outdir, "barcodes.csv"), row.names = FALSE)

  meta <- sub@meta.data[, c("orig.ident", "merged_clusters",
                            "barcode",
                            "nCount_RNA", "nFeature_RNA",
                            "G2M.Score", "S.Score", "Phase",
                            "scrublet_Scores")]
  write.csv(meta, file.path(outdir, "metadata.csv"), row.names = TRUE)
}


subset_and_export(obj2, "10dpci_veh", c("EPDC_A","09","10", "02"), "results/trajectory/cellrank/10dpci_veh_")
subset_and_export(obj2, "10dpci_cCpE", c("EPDC_A","09","10", "02"), "results/trajectory/cellrank/10dpci_cCpE_")
