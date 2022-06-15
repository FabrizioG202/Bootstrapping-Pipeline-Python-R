# %%
import json
import pandas as pd
from utils import extract_bins, RESOLUTIONS


def save_cluster_to_bed(clusters_path : str, cluster_name : str, chromosome : str):
    # Read the bins from the cluster file
    with open(clusters_path, "r") as f:
        clusters = json.load(f)

        # Resolutions and bins
        resolution = RESOLUTIONS[cluster_name.split("_")[0]]
        bins = clusters["cl_member"][cluster_name]
        bins = extract_bins(bins)

        chromosomes = []
        starts = []
        ends = []

        for bin_ in bins:
            start = bin_
            end = int(bin_ + resolution)

            chromosomes.append(chromosome)
            starts.append(start)
            ends.append(end)

    # Create the bed file
    bed_df = pd.DataFrame(
        {
            "chrom": chromosomes,
            "start": starts,
            "end": ends,
        }
    )

    # save the bed file
    bed_df.to_csv(f"{cluster_name}.bed", sep="\t", index=False, header=False)

#%%
# If the file is run as a script, run the main function.
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Convert the clusters to bed format."
    )
    parser.add_argument(
        "clusters_path",
        type=str,
        help="The path to the clusters file.",
    )
    parser.add_argument(
        "cluster_name",
        type=str,
        help="The name of the cluster.",
    )
    parser.add_argument(
        "chromosome",
        type=str,
        help="The chromosome of the cluster.",
    )
    args = parser.parse_args()

    save_cluster_to_bed(args.clusters_path, args.cluster_name, args.chromosome)