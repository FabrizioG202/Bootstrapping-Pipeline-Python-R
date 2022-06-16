# Unibind enrichment
This pipeline provides a way to generate the files necessary to perform Transcription Factor Binding Site enrichment analysis using the [Unibind Enrichment Tool ](https://unibind.uio.no/enrichment/).

It takes as an input...

And outputs 3 files, one containing mapping to enhancer (in the FANTOM5 database) which are contained inside putative transcriptional Hubs and the other containing all the enhancer active in the cell line, that is all the enhancers which map to a peak in the FANTOM5 database. 