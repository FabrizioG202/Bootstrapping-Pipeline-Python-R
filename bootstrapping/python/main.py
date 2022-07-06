import subprocess
import os
import sys
import argparse
from utils import *  # Import the colored printing functions

def __chr_sort_part(x : str):
    p = x.replace("chr", "")
    if p.isdigit():
        return int(p)
    else:
        return p

def sort_by_chromosome(x : dict | list):
    if type(x) == dict:
        # Return a sorted dict.
        return {k: v for k, v in sorted(x.items(), key=lambda item: __chr_sort_part(item[0]))}
    elif type(x) == list:
        return sorted(x, key=lambda x: __chr_sort_part(x))
    else:
        raise TypeError(f"The type of the argument is not a dict or a list. It is: {type(x)}")

if __name__ == "__main__":

    # Read the command line arguments
    parser = argparse.ArgumentParser(description="Compute the p-value of a feature in a genome.")
    parser.add_argument("-f", "--features", help="The path to the file containing the data.", required=True)
    parser.add_argument("-a", "--annotations", help="The path to the file containing the annotation.", required=False, default="../data/counts/")
    parser.add_argument("-g", "--genome", help="The path to the genome file.", required=False, default="../../data/hg19.genome")
    parser.add_argument("-c", "--clusters", help="The path to the cluster folder.", required=False, default="../../data/clusters/")
    parser.add_argument("-n", "--n", help="The number of bootstrap iterations.", required=False, default=1000)
    parser.add_argument("-t", "--t", help="The number of threads.", required=False, default=1)
    parser.add_argument("-u", "--function", help="The path to the folder containing the functions.", required=False, default="../../data/annotations/")
    parser.add_argument("-o", "--output", help="The path to the folder containing the results.", required=False, default="../results/")
    args = parser.parse_args()

    # check that the passed files exist
    if not os.path.isfile(args.genome):
        sys.exit("The genome file does not exist.")
    if not os.path.isdir(args.clusters):
        sys.exit("The cluster folder does not exist.")
    if not os.path.isfile(args.features):
        sys.exit("The features file does not exist.")

    MULTITHREADING = args.t > 1

    # Importing this here since it's heavy and we want the user to known of any previous errors before importing them.
    from src import compute_p_value
    import pandas as pd

    # The path to the feature folder
    feature_path = args.features
    cluster_folder  = args.clusters

    # Annotation file is there.
    annotation_counts_path = args.annotations

    # The resulting dataframe.
    results_df = None

    # the chromosome list is jusst all the stuff from the filen name chr1_spec_res.json from the cluster folder
    chromosomes = [f.split("_")[0] for f in os.listdir(args.clusters) if f.endswith(".json")]
    chromosomes = sort_by_chromosome(chromosomes)

    # Compute the p-value for each chromosome.
    for chromosome in chromosomes:

        # the chromosome file name is like this: chr10_spec_res.json
        # extract the chromosome name (chr10)
        chromo = chromosome.split("_")[0]

        # cluster file is just the first file in the folder starting with the chromosome name
        cluster_file = os.path.join(cluster_folder, chromo + "_spec_res.json")
        
        # the p_values
        p_value: pd.DataFrame = compute_p_value(int(args.n), feature_path, chromo, cluster_file, genome_path=args.genome,
                                                annotation_path=annotation_counts_path, function_path=args.function, batch_size=10, enable_tqdm=True)

        # Add the chromosome to the p_values.
        p_value["chromosome"] = chromo

        # Add the cell type to the p_values.
        results_df = p_value if results_df is None else pd.concat([results_df, p_value])

    # Save the results.
    # create the results folder if it doesn't exist.
    if not os.path.isdir(os.path.dirname(args.output)):
        os.mkdir(os.path.dirname(args.output))


    # keep only the cluster name (name), the p-value and the chromosome and move the chromosome to the first column.
    results_df = results_df[["chromosome", "p_value", "name"]]
    # move the chromosome to the first column.
    results_df = results_df[["chromosome", "name", "p_value"]]
    results_df.to_csv(args.output, sep="\t", index=False)
 