# Unibind Enrichment Pipeline:
This pipeline creates the necessary files to run the enrichment analysis using Unibind Enrichment tool (https://unibind.uio.no/enrichment/). It creates 3 files for each cell line:
- Foreground (The reads inside the given trascriptional hubs (enriched clusters) which belong to enhancer regions)
- Active (All the enhancer reads in the given cell line)
- Full (All the enhancer reads in the FANTOM5 data)

(The TSS pipeliene is not enabled right now)