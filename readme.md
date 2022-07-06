# Cluster Analysis Pipelines based on Bootstrapping
This repository contains the following pipelines:
## Boostrapping pipeline 1 (CAGE Enrichment):
Identify clusters which are enriched in a particular feature set using bootstrapping. To do so, a large amount of random genomes are generated.
Outputs to: _results/enrichment/_

### Filtering:
It is possible to choose the data on which the pipeline is run. To do so edit the config.yaml file in the bootstrapping folder. 
Check that all the feature you need are being printed out when running the pipeline.

## Boostrapping pipeline 2 (CTCF Enrichment):
Produce CTCF (or any other feature) enrichment along the 3D genome.
Outputs to: _results/ctcf_enrichment/_
## P-Value Correction (DAGGER):
Correct the p-values of the bootstrapping pipeline 1 using DAGGER.
Outputs to: _results/qvalues/_
## Downstram Analysis (Downstream Analysis):
Useful scripts for performing downstream analysis on the data produced by the 2 previous pipelines.

# Results file structure:
For the structure of the results files, see the [results_file_structure.png](results_file_structure.png) file.