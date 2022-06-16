# Cluster Analysis Pipelines based on Bootstrapping
This repository contains the following pipelines:
## Boostrapping pipeline 1 (CAGE Enrichment):
Check for clusters which are enriched in a particular feature set using bootstrapping. To do so, a large amount of random genomes are generated.
Outputs to: _results/pvalues/_
## Boostrapping pipeline 2 (CTCF Enrichment):
Produce CTCF (or any other feature) enrichment along the 3D genome.
Outputs to: _results/ctcf_enrichment/_
## P-Value Correction (DAGGER):
Correct the p-values of the bootstrapping pipeline 1 using DAGGER.
Outputs to: _results/qvalues/_
## Downstram Analysis (Downstream Analysis):
Useful scripts for performing downstream analysis on the data produced by the 2 previous pipelines.