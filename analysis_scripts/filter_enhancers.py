#%%
import pyranges as pr
import sys
import pandas as pd
from src.crayon import *

def read_peak_type(path : str, ) -> tuple[list[str], list[str]]:
    """
    Read the peak type file.
    """
 
    ### The file has the following structure: its a TSV file with 2 column: id and peak type, separated by a tab.
    ### We want to return two lists, one with the ids which correspond to peak type "enhancer" and one with the ids which correspond to peak type "tss".

    data : pd.DataFrame = pd.read_csv(path, sep="\t")
    data = data.set_index("ID")

    enhancer_ids = data.loc[data["type"] == "enhancer"].index.tolist()
    tss_ids = data.loc[data["type"] == "tss"].index.tolist()

    return enhancer_ids, tss_ids

def save_fantom_ids_as_bed_file(ids : list[str], save_path : str) -> None:
    ### Each peak comes in the form chr10:3510310-3510716, save it to a bed file.
    with open(save_path, "w") as f:
        for id in ids:
            chromo = id.split(":")[0]
            start, end = id.split(":")[1].split("-")

            f.write("{}\t{}\t{}\t{}\n".format(chromo, start, end, id))
#%%


if __name__ == "__main__":  
    argvs = sys.argv
    
    match argvs:
        case [_, "fromBed", bed_file, feature_file, peak_type_file,]:
            ### Reading the peak type file.
            peak_type_ids, tss_ids = read_peak_type(peak_type_file)
            
            ### Reading the clusters
            clusters : pr.PyRanges = pr.read_bed(bed_file)
                
            ### Reading the features (to overlap)
            _features : pr.PyRanges = pr.read_bed(feature_file)
            enhancers : pr.PyRanges = _features[_features.Name.isin(peak_type_ids)]
            all_enhancers_ids = set(enhancers.Name.values.flatten().tolist())

            ### Overlapping them = finding the features which are in the clusters
            enhancers_in_clusters : pr.PyRanges = enhancers.overlap(clusters) 
            enhancers_in_clusters_ids = enhancers_in_clusters.Name.values.flatten().tolist() # getting only the values
            enhancers_in_clusters_ids = set(enhancers_in_clusters_ids) # removing the duplicates

            ### Log the number of enhancers in the clusters vs the number of enhancers in the features
            print(f"Number of enhancers in the clusters: {len(enhancers_in_clusters)}, number of enhancers in the features: {len(all_enhancers_ids)}") 

            ### Get the enhancers which are not in the clusters:
            enhancers_not_in_clusters = all_enhancers_ids - enhancers_in_clusters_ids

            ### Convert both to bed files:
            in_clusters_save_path = bed_file.replace(".bed", "_enhancers_in_clusters.bed")
            not_in_clusters_save_path = bed_file.replace(".bed", "_enhancers_not_in_clusters.bed")

            save_fantom_ids_as_bed_file(enhancers_in_clusters_ids, in_clusters_save_path)
            print(f"Saved {blue(in_clusters_save_path)}")
            save_fantom_ids_as_bed_file(enhancers_not_in_clusters, not_in_clusters_save_path)
            print(f"Saved {blue(not_in_clusters_save_path)}")

        case [_, "info" | "i" | "--info" | "--i" | "--help" | "-h"]:
            print("Filter enhancers from a bed file.\n")
            print("Usage:")
            print("\tpython filter_enhancers.py fromBed <bed_file> <feature_file> <peak_type_file>")
            print("\n")
            print("Output:")
            print("\t<bed_file>_enhancers_in_clusters.bed")
            print("\t<bed_file>_enhancers_not_in_clusters.bed")
            print("\n")