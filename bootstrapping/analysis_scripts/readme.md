# Available Analysis workflows:
## Clusters Data I/O and conversions 
This file provides basic utilities to convert between different data formats to store cluster information.
1. python clus_io.py tsv2clus <source_path> [out_path = in.CLUS] [chromo_column = 0] [name_column = 2]
2. python clus_io.py clus2bed <cluster_folder> <source_path> [out_path = in.BED]

## Cluster Utilities:
1. python clus_utils.py intersectClus <out_path> file1 file2 ...
2. python clus_utils.py sortByChromo <clus_file>

## Lift Over:
This file provides utilities to perform lift over between different genome assemblies.
1. python lift_over.py <sourc_path> <from_assembly> <to_assembly> [out_path = in_assembly.BED]

## Parse Results:
This file contains utilities for parsing the results of the pipelien as well as filtering them by any metric.
1. python parse_results.py applyFDR <source_path> <FDR_MODE> [save_path]
2. python parse_results.py filter <source_path> <parameter> <action> <cut_off> [save_path]
3. python parse_results.py 2Clus <source_path> [save_path]

## Filter enhancers:
This file will output 2 files one containing the enhancers that are in the clusters and one containing the enhancers that are not in the clusters.
1. python filter_enhancers.py fromBed <source_cluster_bed> <features_file> <peak_type_file> [out_folder]
