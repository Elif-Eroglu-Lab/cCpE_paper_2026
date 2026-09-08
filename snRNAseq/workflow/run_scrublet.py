################# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: Run Scrublet on each sample

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ################# Load libraries
    import scrublet as scr
    import scipy.io
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import h5py
    import scipy.sparse as sp
    import scanpy
    import pandas as pd

    # Set seed for reproducibility
    np.random.seed(42)

    ############### Load Snakemake variables
    # input_h5ad = "../results/CB2Scrublet/counts_matrix.h5ad"
    # output_h5ad = "../results/CB2Scrublet/scrubbed_counts_matrix.h5ad"
    # output_tsv = '../results/CB2Scrublet/scrublet_results.csv'

    input_h5ad = snakemake.input["input_h5ad"]

    output_h5ad = snakemake.output["output_h5ad"]
    output_csv = snakemake.output["output_csv"]

    plot1 = snakemake.output["plot1"]
    plot2 = snakemake.output["plot2"]

    expected_doublet_rate = 0.08
    min_gene_variability_pctl = 85

    ################# Load counts matrix
    # Load the raw counts matrix as a scanpy sparse matrix with 
    # cells as rows and genes as columns

    adata = scanpy.read_h5ad(input_h5ad)

    # Compute QC metrics on the full matrix
    cell_sums = np.array(adata.X.sum(axis=1)).ravel()        # nCount_RNA
    cell_nnz  = np.array((adata.X > 0).sum(axis=1)).ravel()  # nFeature_RNA
    gene_nnz  = np.array((adata.X > 0).sum(axis=0)).ravel()  # cells per gene

    # Build masks
    cell_mask = (cell_sums < 50000) & (cell_nnz > 200)
    gene_mask = gene_nnz >= 3

    #  Apply both masks
    adata = adata[cell_mask, :][:, gene_mask]

    # Get quick summary of the adata object (especially the dimensions)
    adata

    counts_matrix = adata.layers["cellbender"]

    ################## Run Scrublet

    # Initialize Scrublet object
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)

    # Simulate doublets
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                            min_cells=3, 
                                                            min_gene_variability_pctl=min_gene_variability_pctl, 
                                                            n_prin_comps=30)

    # Plot histogram
    scrub.plot_histogram();
    plt.savefig(plot1)

    # Run UMAP
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')

    # Plot scores on histogram
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(plot2)

    # Store results in AnnData
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    # Write the new adata object
    adata.write(output_h5ad)

    # Extract relevant information
    cells = adata.obs_names # Cells
    barcodes = pd.DataFrame(adata.obs_names.to_list(), columns=["barcode"])
    results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
    scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")

    # Save information 
    dataframe = pd.concat([barcodes, results, scores], axis=1)
    dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "Doublet")
    dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "Singlet")

    dataframe.to_csv(output_csv, sep = ";", index = False)