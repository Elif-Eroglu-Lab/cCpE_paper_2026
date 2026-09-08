################## README
# Author: Inika Prasad
# Aim: To make the output of CellBender compatible with Scrublet
# i.e. from .h5 format to .h5ad

################### Setup logfile
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ################### Load libraries
    from cellbender.remove_background.downstream import load_anndata_from_input_and_output

    #################### Define input & output 
    # cb_input = "../../../data/raw/DDSdownload/10dpci_tester/10dpci_DDSdelivery/files/P36258/01-SC-Analysis/E.Eroglu_25_03/P36258_101/outs/raw_feature_bc_matrix.h5"
    # cb_output = "../results/01_remove_ambient_RNA/10dpci_veh_filtered.h5"
    # conversion_output = "../results/CB2Scrublet/counts_matrix.h5ad"

    cb_input = snakemake.input["cb_input"]
    cb_output = snakemake.input["cb_output"]

    conversion_output = snakemake.output["conversion_output"]

    print(cb_input)
    print(cb_output)
    print(conversion_output)

    #################### Load data
    adata = load_anndata_from_input_and_output(input_file = cb_input,
                                            output_file = cb_output)

    ################### Write data
    # Write to an h5ad format compatible with Scrublet
    adata.write(conversion_output)