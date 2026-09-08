# README
Author: Inika Prasad
Aim: Attempt to construct Snakemake pipeline to ease the parallel & iterative analyses of single-nuclei RNA sequencing data

## Environments:
The global conda environment and each rule's specific conda environment are specified in the Snakefile. They are created and activated when the pipeline is run.
But to try dry runs, you need the environment with Snakemake up and running:
mamba env create --name snakemake_snRNAseq --file envs/global.yaml
conda activate snakemake_snRNAseq

## Dry run:
snakemake -n -p --configfile config/config.yaml

## Execute the job
Modify core number based on your computer. Good practice to use no more than 75-90% of total cores.
I have 24 cores, so I'll try with 17 cores
snakemake --cores 17 --configfile config/config.yaml --use-conda --sdm conda
snakemake results/01_remove_ambient_RNA/PM_naive.h5 --cores 17 --use-conda --sdm conda 
snakemake results/01_remove_ambient_RNA/Early_active.h5 --cores 17 --use-conda --sdm conda 
snakemake results/01_remove_ambient_RNA/14dpSham.h5 --cores 17 --use-conda --sdm conda 


## Visualize the workflow
snakemake filenameX.whatever --dag | dot -Tsvg > DAGforhowfilenameXwasproduced.svg
snakemake --dag | dot -Tsvg > DAG.svg

## Versions of the workflow
CellBender takes a lot of time to run, so for the sake of testing & troubleshooting, I made Snakefile_nocellbender
and config_nocellbender.yaml. They use the filtered feature matrices. 

### for a dry run
snakemake -s Snakefile_nocellbender -n -p --configfile config/config_nocellbender.yaml

### for all samples
snakemake -s Snakefile_nocellbender --cores 17 --configfile config/config_nocellbender.yaml --use-conda --sdm conda

### for only one sample, to test
snakemake -s Snakefile_nocellbender --cores 17 --configfile config/config_nocellbender_1sampletest.yaml --use-conda --sdm conda

## Structure of Snakemake rules:
rule myrulename:
    input:
    output:
    conda:
    log:
    script: