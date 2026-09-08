# Combine QC metrics from metrics.csv

# View ambient RNA profile

# Pull relevant information from the log file

from cellbender.remove_background.report import generate_summary_plots
import os
input_file = "../../../data/raw/dataNCBI_lighter/E.Eroglu_25_01/P34662_102/outs/raw_feature_bc_matrix.h5"
output_file = "./../results/01_remove_ambient_RNA/PM_naive.h5"

os.path.exists(output_file)
os.path.exists(input_file)

generate_summary_plots(input_file=input_file, 
                       output_file=output_file,
                       truth_file=truth_file)