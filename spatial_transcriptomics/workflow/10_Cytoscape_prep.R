########################### README ###########################
# Author: Inika Prasad, inika.prasad@gmail.com
# Date: 2026-04-20

# Get human identifiers for the CLuster 10 (border muscle DEGs),
# for use in a Cytoscape network

########################## Load libraries ########################
library(readxl)
library(tidyverse)
library(writexl)

########################### Load data ##############

# Read in DEG data
df <- read_xlsx("../results/08c_refinedDEGS/ZANIA_ZAURAK_refined.xlsx",
                sheet = "ZANIA_ZAURAK_10")

######################### Load identifiers ######################
identifiers <- read_xlsx("../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.xlsx")

#################### Get identifiers for DEGs ###################

df2 <- left_join(df, identifiers, join_by(gene == genename_Hs_PwNameFormat))

# Compare human gene names and Pw-compatible gene names
df2 %>% dplyr::select(gene, genename_Hs) %>% View()

# Human gene names
df2 %>% dplyr::pull(genename_Hs) %>% cat()

################### Filter #####################################
df2 %>% filter(abs(avg_log2FC) > 0.5) %>% write_xlsx("../results/10_Cytoscape/BorderZone/log2FC0.5_intersections_ZANIA_ZAURAK.xlsx")
df2 %>% write_xlsx("../results/10_Cytoscape/BorderZone/noFCcutoff_intersections_ZANIA_ZAURAK.xlsx")


df2 %>% filter(abs(avg_log2FC) > 0.5) %>% dplyr::pull(genename_Hs) %>% cat()
df2 %>% dplyr::pull(genename_Hs) %>% cat()
