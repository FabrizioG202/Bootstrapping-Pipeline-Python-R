# Bootstrapping Pipeline for cluster enrichment assessment.
The following pipeline allows to assess the enrichment of a cluster with respect to a list of provided features. It can work on multilple cell lines and multiple feature types at the same time.
It is currently split between a python and an R version, to run either one, run _snakemake compute_pvalue_ in the corresponding directory.

## Python: 
The python version is based upon PyRanges and NCLS.
## R:
The R version is based upon the R package GenomicRanges.

## General workflow:
Each Version of the pipeline contains a config file to specify the number of bootstraps. _The 2 different pipelines produce different temporary files right now, so they are not interoperable._

# Data Structure:
The main data folder has the following structure:

    data:
        features:
            - cell_line_1:
                - feature_1:
                    - file_1.csv
                - feature_2:
                    - file_1.csv
            - cell_line_2:
                - feature_1:
                    - file_1.csv
                - feature_2:
                    - file_1.csv
            - ...

        clusters:
            - cell_line_1:
                - chr1....json
                - chr2....json
                - ...
            - cell_line_2:
                - chr1....json
                - chr2....json
                - ...
            - ...
        annotations:
            - fn_BED:
                - chr1: 
                    - utr3.bed
                    - ...
                - chr2:
                    - utr3.bed
                    - ...
                - ...
            - hg19.genome

The pipeline will generate a table in the following format "CELL_LINE.FEATURE.tsv". The table will contain the following columns by default: chromosome, cluster_name, p_value.

# Available Analysis workflows:
The analysis folder contains script which can be used to perform downstream analysis of the results produced with the pipeline. Each script expects by default the scripts to be present in the ./results/ folder (the same as the pipeline output one) but a custom path can be passed via command line arguments. (Any of the script is available both as a command-line python script and as a jupyter notebook).
## FDR and filtering (_fdr_and_filtering_):
This script allows to parse the tables and perform a multiple testing correction. It also allows to filter the results based on a given p-value, q-value (after applying FDR) and chromosome threshold.
## Gene Set Enrichment Analysis (_gene_set_enrichment_):
This script allows to perform a gene set enrichment analysis on the results. It expects the gene set to be present in the ./analysis/data/ folder.
## Plotting (_TODO_):
This script allows to plot the results of the pipeline. It expects the results to be present in the ./results/ folder.
## Enhance/TSS enrichment analysis (_TODO_):
This script allows to perform a peak category enrichment analysis on the results. It expects the gene set and the peak/TSS-Enhancer table to be present in the ./analysis/data/ folder.