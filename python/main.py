import subprocess, os, sys
import argparse
from utils import * ### Import the colored printing functions

if __name__ == "__main__":

    ### Read the command line arguments
    parser = argparse.ArgumentParser(description="Compute the p-value of a feature in a genome.")
    parser.add_argument("-f", "--feature", help="The path to the file containing the data.", required=False, default="../data/features/")
    parser.add_argument("-a", "--annotation", help="The path to the file containing the annotation.", required=False, default="../data/counts/")
    parser.add_argument("-g", "--genome", help="The path to the genome file.", required=False, default="../data/hg19.genome")
    parser.add_argument("-c", "--cluster", help="The path to the cluster folder.", required=False, default="../data/clusters/")
    parser.add_argument("-n", "--n", help="The number of bootstrap iterations.", required=False, default=1000)
    parser.add_argument("-t", "--t", help="The number of threads.", required=False, default=1)
    parser.add_argument("-u", "--function", help="The path to the folder containing the functions.", required=False, default="../data/annotations/")
    parser.add_argument("-r", "--result", help="The path to the folder containing the results.", required=False, default="../results/")
    args = parser.parse_args()

    # check that the genome file exists
    if not os.path.isfile(args.genome):
        print(red("The genome file does not exist."))
        sys.exit(1)

    # check that the cluster folder exists
    if not os.path.isdir(args.cluster):
        print(red("The cluster folder does not exist."))
        sys.exit(1)

    # check that the feature file exists
    if not os.path.isdir(args.feature):
        print(red("The feature folder does not exist."))
        sys.exit(1)

    # The folder structure is the following:
    # the cluster folder contains a subfolder for each cell type.
    cluster_cell_types = os.listdir(args.cluster)

    # The feature folder contains one subfolder for each cell type and each cell type subfolder contains one
    # subfolder for each feature type.
    feature_cell_types = os.listdir(args.feature)
    feature_types = {cell_type: os.listdir(os.path.join(args.feature, cell_type)) for cell_type in feature_cell_types}

    # Check that each feature cell type is also in the cluster cell types
    for cell_type in feature_cell_types:
        if cell_type not in cluster_cell_types:
            print("No cluster folder for {}".format(cell_type))
            sys.exit(1)

    MULTITHREADING = args.t > 1

    # Importing this here since it's heavy and we want the 
    from src import compute_p_value
    import pandas as pd

    # Main loop.
    for cell_type in feature_cell_types:

        # The cluster path.
        cluster_path = os.path.join(args.cluster, cell_type)

        # The available chromosomes.
        available_chromosomes = os.listdir(cluster_path)

        # Log the amount of chromosomes found.
        log_info(f"Found {len(available_chromosomes)} chromosomes for {cell_type}, the p-values will be computed for each chromosome.")

        # Looping over the path.
        for feature_type in feature_types[cell_type]:

            # Log what we are doing
            print("Computing p-value for feature type {} in cell type {}".format(blue(feature_type), yellow(cell_type)))

            # The path to the feature folder
            feature_path = os.path.join(args.feature, cell_type, feature_type)

            # Look for the annotation counts file
            # It is the fist file ending in tsv.
            annotations_folder = os.path.join(args.annotation, cell_type, feature_type)
            annotation_file = [file for file in os.listdir(annotations_folder)][0]

            # Annotation file is there.
            annotation_counts_path = os.path.join(annotations_folder, annotation_file)

            # the feature file is the first in the feature_path folder which is not the annotation counts file.
            feature_filename = [file for file in os.listdir(feature_path)][0]
            feature_file = os.path.join(feature_path, feature_filename)

            results_df = None

            if not MULTITHREADING:
                # Compute the p-value for each chromosome.
                for chromosome in available_chromosomes:
                    

                    # the chromosome file name is like this: chr10_spec_res.json
                    # extract the chromosome name (chr10)
                    chromo = chromosome.split("_")[0]

                    # the path to the full chromosome.
                    full_cluster_path = os.path.join(cluster_path, chromosome)

                    # the p_values
                    p_value: pd.DataFrame = compute_p_value(int(args.n), feature_file, chromo, full_cluster_path, genome_path=args.genome, annotation_path=annotation_counts_path, function_path=args.function, batch_size=10, enable_tqdm=True)

                    # Add the chromosome to the p_values.
                    p_value["chromosome"] = chromo

                    # Add the cell type to the p_values.
                    results_df = p_value if results_df is None else pd.concat([results_df, p_value])
            else:
                # exit, multiprocessing is not supported on windows
                print(red("Multiprocessing is not supported on windows."))
                sys.exit(1)

            # Save the results.
            # create the results folder if it doesn't exist.
            if not os.path.isdir(args.result):
                os.mkdir(args.result)
            
            # THE SAVE filename should be the following:
            # cell_type.feature_type.tsv
            save_filename = "{}.{}.tsv".format(cell_type, feature_type)
            save_path = os.path.join(args.result, save_filename)

            # keep only the cluster name (name), the p-value and the chromosome and move the chromosome to the first column.
            results_df = results_df[["chromosome", "p_value", "name"]]
            # move the chromosome to the first column.
            results_df = results_df[["chromosome", "name", "p_value"]]
            results_df.to_csv(save_path, sep="\t", index=False)
            log_info("Saved results to {}".format(save_path))