library(readxl)
library(dplyr)
library(openxlsx)

##########################
unique_code <- c("ALGORAB", "ALGOL")

# Get the relevant files
folder1 <- paste0("results/DGE/", unique_code)
files1 <- list.files(folder1, pattern = "\\.xlsx$", full.names = TRUE)
files1 <- files1[!files1 %in% paste0("results/DGE/", unique_code, "/summary_DGE_results.xlsx")]
files1

# Specify order in which conditions should be considered/appear in heatmap
df <- read_xlsx("results/figures/DEG_sheetkey.xlsx")
order_conditions <- df$unique_code

################## Collect log2FC information in the correct order
# Empty list to collect information about DEGs
# Create a new workbook
wb <- createWorkbook()

# For every condition
for (i in seq_along(order_conditions)) {

    # Determine which file has info about that condition
    degname <- order_conditions[i]
    sheetname <- df[df$unique_code == degname,]$`Sheet name`
    degfile <- files1[grep(degname, files1)]

    # if there actually is a file
    if (length(degfile) != 0) {
        # read in the significant degs for that file/condition
        sig_degs <- read_xlsx(path = degfile,
                              sheet = "sig_degs") %>%
        as.data.frame()

    # Add a worksheet (sheetname can include spaces)
    addWorksheet(wb, sheetname)

    # Write data to that sheet
    writeData(wb, sheet = sheetname, sig_degs)
    }


}

saveWorkbook(wb, "results/figures/EpiPeelout_all_degs.xlsx", overwrite = TRUE)
