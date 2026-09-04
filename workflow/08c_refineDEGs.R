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

############################### Load Snakemake variables

vehccpe <- snakemake@input[["vehccpe"]]
veh1veh2 <- snakemake@input[["veh1veh2"]]

output_refdegs <- snakemake@output[["refdegs"]]
############################### Load data

vehccpe_df <- read_xlsx(vehccpe, sheet = "sig_degs")
veh1veh2_df <- read_xlsx(veh1veh2, sheet = "sig_degs")

vehccpe_info <- read_xlsx(vehccpe, sheet = "run_details")
veh1veh2_info <- read_xlsx(veh1veh2, sheet = "run_details")

############################### Process differences

# Identify DEGs that are only in vehicle vs ccpe
truediff <- setdiff(vehccpe_df$gene, veh1veh2_df$gene)

# extract all info about the refined gene list in vehicle-ccpe
refined_degs <- vehccpe_df %>% filter(gene %in% truediff)

# Add info about up or downregulated
refined_degs <- refined_degs %>%
  mutate(direction = case_when(
    avg_log2FC < 0 ~ "downregulated",
    avg_log2FC > 0 ~ "upregulated",
    TRUE        ~ "neither"
  ))

############################### Save results

write_xlsx(list(refined_degs = refined_degs,
                info = rbind(vehccpe_info, veh1veh2_info)),
           path = output_refdegs)

############################### Close log file
sink(type = "message")
sink()
close(logfile)