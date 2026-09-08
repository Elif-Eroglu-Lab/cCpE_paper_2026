# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: To run Harmony, CCA, and RPCA integration

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################## Load libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(harmony)

set.seed(100)

############################## Load Snakemake variables
input_rds <- snakemake@input[["rds"]]
plot0 <- snakemake@input[["plot0"]]
output_rds <- snakemake@output[["rds"]]

############################## Load merged (unintegrated) object
obj <- readRDS(input_rds)

################################ Harmony integration
# Integrate samples with Harmony
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony"
)
# ################################ CCA integration
# # Integrate samples with CCA
# obj <- IntegrateLayers(
#   object = obj, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca"
# )
# ################################ RPCA integration
# # Integrate samples with RPCA
# obj <- IntegrateLayers(
#   object = obj, method = RPCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.rpca"
# )

#################################### Save object post-integrations

saveRDS(obj, output_rds)

#################################### Close log file
sink(type = "message")
sink()
close(logfile)