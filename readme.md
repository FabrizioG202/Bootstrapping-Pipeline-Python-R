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