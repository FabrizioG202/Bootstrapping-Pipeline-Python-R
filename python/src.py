# %%
import pandas as pd
import numpy as np
import pyranges as pr
import numba
import json
from typing import Literal
import random
import ncls
import time
import os
from tqdm.autonotebook import tqdm

# %%
# The resolutions
RESOLUTIONS = {
    "5kb": 5e3,
    "10kb": 1e4,
    "50kb": 5e4,
    "100kb": 1e5,
    "500kb": 5e5,
    "1Mb": 1e6,
}

# %%
### Joins the ranges which are close to each other
def __join_contiguous(starts: list[int], resolution: int):
    starts_out = []
    ends_out = []
    index = 0
    while index < len(starts):
        starts_out.append(starts[index])
        index += 1
        while index < len(starts) and starts[index] - starts[index - 1] <= resolution:
            index += 1
        ends_out.append(starts[index - 1] + resolution)

    return starts_out, ends_out


### Imports the clusters.
def import_clusters(path: str, *, merged: bool = False, keep_clusters: list[str] = None, keep_indices: list[int] = None, join_contiguous: bool = False) -> tuple[pd.DataFrame, np.ndarray]:
    """Reads the clusters json file and returns the clusters.

    Args:
        path (str): The path to the file.
        merged (bool, optional): Wether to return the clusters as merged (one line per cluster) or non-merged "exploded" (one line per bin). Defaults to False.

    Returns:
        pr.PyRanges: The clusters.
    """

    # Convert the last column to a list of string to a list of int
    # "1" -> [1]
    # "[1,2,3]" -> [1,2,3]
    def __extract_bins(x: str | list) -> list[int]:
        if type(x) == str:
            return [int(x)]
        elif type(x) == list:
            return [int(i) for i in x]
        else:
            raise TypeError(f"The type of the argument is not a string or a list. It is: {type(x)}")

    # Loads the json and retrieve only the information about the clusters.
    clusters: dict[str, object] = json.load(open(path))["cl_member"]

    # Convert the members to a list of int
    clusters: dict[str, list[int]] = {k: __extract_bins(v) for k, v in clusters.items()}

    # convert the clusters to a dataframe
    clusters_df = pd.DataFrame(columns=["name", "members"])
    clusters_df["name"] = clusters.keys()
    clusters_df["members"] = clusters.values()  # the start values
    # Transform the clusters in ranges
    clusters_df["resolution"] = clusters_df["name"].apply(lambda x: int(RESOLUTIONS[x.split("_")[0]]))

    # Keep only the clusters that are requested
    if keep_clusters is not None:
        clusters_df = clusters_df[clusters_df["name"].isin(keep_clusters)]
        # reset the index
        clusters_df.reset_index(drop=True, inplace=True)
    elif keep_indices is not None:
        clusters_df = clusters_df.iloc[keep_indices]
        # reset the index
        clusters_df.reset_index(drop=True, inplace=True)

    # If the clusters are requested to be merged, return the merged clusters
    if merged:
        return clusters_df, np.array([])

    if not join_contiguous:
        # Otherwise, explode the clusters
        clusters_df = clusters_df.explode("members")
        indices = clusters_df.index.values  # the indices needed to re-merge the clusters

        #  # Create a column for the resolution (extracted from the name)
        clusters_df["end"] = clusters_df["members"] + clusters_df["resolution"]
        clusters_df.drop(columns=["resolution"], inplace=True)  # Drop the resolution column
        clusters_df.reset_index(drop=True, inplace=True)  # Reset the index

        # Rename the columns
        clusters_df.rename(columns={"members": "start"}, inplace=True)

    else:
        # Join the contiguous clusters
        starts, ends = zip(*clusters_df["members"].apply(lambda x: __join_contiguous(x, clusters_df["resolution"].iloc[0])))
        clusters_df["start"] = starts
        clusters_df["end"] = ends
        clusters_df.drop(columns=["members"], inplace=True)

        # explode the clusters
        clusters_df = clusters_df.explode(["start", "end"])
        indices = clusters_df.index.values  # the indices needed to re-merge the clusters

        # Drop the resolution column
        clusters_df.drop(columns=["resolution"], inplace=True)

        # reset the index
        clusters_df.reset_index(drop=True, inplace=True)

    return clusters_df, indices


# %%
def unique_void_view(a):
    return (
        np.unique(a.view(np.dtype((np.void, a.dtype.itemsize * a.shape[1]))))
        .view(a.dtype)
        .reshape(-1, a.shape[1])
    )

def compress_ints(a, b):
    """Compresses 2 ints together using bitshift"""
    return (a << 32) | b

def decompress_ints(a):
    """Decompresses 2 ints from compressed int using bitshift"""
    return (a >> 32), a & 0xFFFFFFFF
    

### the overlap counter.
class OverlapCounter:

    cluster_ncls: ncls.NCLS64
    step: int
    unique_cluster_count: int
    cluster_indexer: np.ndarray

    def __init__(self, clusters: pd.DataFrame, cluster_indexer: np.ndarray, *, step: int = 25):

        # Check that the clusters dataframe contains the columns start and end
        assert "start" in clusters.columns and "end" in clusters.columns, "Clusters must contain the columns start and end"

        starts = clusters.start.values.flatten().astype(np.int64)
        ends = clusters.end.values.flatten().astype(np.int64)
        indices = np.arange(len(starts), dtype=np.int64)

        # assert that the clusters array has 3 columns
        assert clusters.shape[1] == 3, "Please provide a clusters array with 3 columns [start, end, index]"

        self.cluster_ncls = ncls.NCLS(starts, ends, indices)
        self.step = step
        self.unique_cluster_count = np.max(cluster_indexer) + 1
        self.cluster_indexer = cluster_indexer

    def count_overlaps(self, features: np.ndarray) -> np.ndarray:
        assert features.shape[1] == 4, "Please provide a features array with 4 columns [start, end, index, batch]"

        # Getting the number of batches by finding the maximum batch index
        n_batches = features[-1][3] + 1

        # Initializing the overlap counts
        counts_ = np.zeros((self.unique_cluster_count, n_batches), dtype=np.int64)  # [cluster, batch]

        # split the features based on the batch index (the last column)
        batch_splits = np.split(features, np.where(np.diff(features[:, 3]) != 0)[0] + 1, axis=0)
        for i, batch in enumerate(batch_splits):
            l_idx_, r_idx_ = self.cluster_ncls.all_overlaps_both(batch[:, 0].T, batch[:, 1].T, batch[:, 2].T)
    
            r_idx_ = self.cluster_indexer[r_idx_]

            # this is a better way of finding the unique indices.
            cc_idxs = compress_ints(l_idx_, r_idx_) # compress the indices
            uu_idsx = pd.unique(cc_idxs) # get the unique indices
            u_l_idxs, u_r_idxs = decompress_ints(uu_idsx)

            OverlapCounter.__update_counts(counts_, u_l_idxs, u_r_idxs, batch[0, 3])

        return counts_

    @staticmethod
    @numba.jit(nopython=True)
    def __update_counts(counts_: np.ndarray, l_idx: np.ndarray, r_idx: np.ndarray, batch: int) -> None:
        for _, r_id in zip(l_idx, r_idx):

            # transforming the cluster index
            counts_[r_id, batch] += 1


# %%
### The Shuffler.


@numba.experimental.jitclass()
class RegionDelimiter:

    # The data requested for the shuffler
    intervals: numba.int64[:, :]
    weights: numba.float64[:]
    cumsum_weights: numba.float64[:]
    count: numba.int64

    # The display only data
    # The mode: -1 is exclusion, 1 is inclusion.
    mode: numba.int64

    def __init__(self, intervals, *, mode, count):
        # The intervals are in the shape (start, end)
        self.intervals = intervals
        lenghts: np.ndarray = intervals[:, 1] - intervals[:, 0]
        self.weights: np.ndarray = lenghts / lenghts.sum()
        self.cumsum_weights = self.weights.cumsum()

        # The mode.
        self.mode = mode
        self.count = count


class BedShuffler:

    ### the data to shuffle
    features: np.ndarray

    ### The information about the chromosome.
    chromosome: str
    chromosome_length: int

    ### The delimiters to use for the shuffling
    delimiters: list[RegionDelimiter]

    @staticmethod
    def pyranges_to_numpy(bed: pr.PyRanges, *, last_column: Literal["index", "lenght"] | None | int):
        data = np.asarray(np.stack((bed.Start.values.flatten(), bed.End.values.flatten()), axis=1), order="f", dtype=np.int64)

        if last_column == "index":
            data = np.column_stack((data, np.arange(data.shape[0])))
        elif last_column == "lenght":
            data = np.column_stack((data, bed.End.values.flatten() - bed.Start.values.flatten()))
        elif isinstance(last_column, int):
            data = np.column_stack((data, np.full(len(bed), last_column)))

        # Make the data f-contiguous
        data = np.asarray(data, order="f")
        return data

    def __init__(self, features: np.ndarray, chromosome: str, chromosome_length: int):
        self.features = features
        self.chromosome = chromosome
        self.chromosome_length = chromosome_length

        self.delimiters = numba.typed.List()
        #self.delimiters = []
        
    ### Add an include file. (-incl)
    def add_include(self, path: str, count: int = 1):
        self.delimiters.append(
            RegionDelimiter(
                BedShuffler.pyranges_to_numpy(pr.read_bed(path).merge(), last_column=None),
                mode=1,
                count=count,
            )
        )

    ### Add an exclude file (-excl)
    def add_exclude(self, path: str, count: int = 1):
        incl = pr.read_bed(path).merge()
        incl = pr.PyRanges(chromosomes=[self.chromosome], starts=[0], ends=[self.chromosome_length]).subtract(incl)
        self.delimiters.append(
            RegionDelimiter(
                BedShuffler.pyranges_to_numpy(incl, last_column=None),
                mode=-1,
                count=count,
            )
        )

    ### Shuffle the file.
    def shuffle(self, *, batch_size: int, max_tries: int = 100):
        return BedShuffler.shuffle_impl(self.features, self.delimiters, batch_size, max_tries)

    @staticmethod
    @numba.njit
    def shuffle_impl(features: np.ndarray, delimiters: list[RegionDelimiter], batch_size: int, max_tries: int):

        # the shuffled entries, in the shape
        # (start, end, index, batch_id)
        shuffled_entries = np.zeros((len(features) * batch_size, 4)[::-1], dtype=np.int64).T
        
        index = 0
        # Prepare the intervals:
        for j in range(batch_size):
            for delim in delimiters:

                # sample count indices from the features.
                _count = delim.count
                replace = _count > len(features)

                # Generating the indices
                indices = np.random.choice(np.arange(len(features)), _count, replace=replace)

                # Sampling the features
                sampled_features = features[indices]

                # looping over the sampled features
                for entry in sampled_features:

                    # The entry to shuffle.
                    length = entry[1] - entry[0]

                    # Whether the entry is valid.
                    inbounds = False

                    # The current try.
                    tries_ = 0

                    # Looping
                    while not inbounds and tries_ < max_tries:
                        tries_ += 1

                        # Sample the interval
                        interval_idx = np.searchsorted(delim.cumsum_weights, np.random.uniform(0, 1), side="left")
                        interval = delim.intervals[interval_idx]

                        # Sample the new start.
                        start = random.randint(interval[0], interval[1])

                        # Checking that the interval fits.
                        enclosed = interval[0] <= start and interval[1] >= start + length
                        if not enclosed:
                            continue

                        inbounds = True

                    if inbounds:
                        shuffled_entries[index] = np.array([start, start + length, index, j])
                    else:
                        shuffled_entries[index] = np.array([entry[0], entry[1], index, j])

                    index += 1

        return shuffled_entries


# %%
### The Odds Ratio Calculator and the Compute P-Value Function.
@numba.njit()
def calculate_odds_ratio(counts: np.ndarray, observed: np.ndarray) -> np.ndarray:
    # The resulting array of odds ratios
    out = np.zeros(counts.shape[0])

    # Looping over the matrix:
    # The shape is currently this.
    # batch 0 | batch 1 | batch 2 | batch 3 | batch 4 | batch 5 | batch 6 | batch 7 | batch 8 | batch 9
    # 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0
    # 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0
    for r in range(counts.shape[0]):

        # We store in the ith position of the array the number of times the observed overlaps are more than the
        # randomly generated ones.
        out[r] = np.sum(counts[r, :] >= observed[r])

    # Returning the array of odds ratios
    return out


def read_annotations(path: str, _chromo: str):
    # Reading the annotations from the annotation.tsv file - (boilerplate code)
    # + Renaming the columns to remove the placeholder j
    _annotations = pd.read_csv(path, sep="\t")  # reading the file
    _annotations.rename(
        columns={k: k.replace("j", "").replace(".BED", "") for k in _annotations.columns},
        inplace=True,
    )
    _annotations.set_index("chr", inplace=True)  # setting to be indexed by the chromosome
    _annotations = _annotations.loc[_chromo]  # Selecting the annotations relative to the chromosome
    # Making it a dictionary for easier access.
    _annotations = _annotations.to_dict()
    return _annotations

def red(text: str):
    return "\033[31m" + text + "\033[0m"

def compute_p_value(n: int, features_path: str, chromo: str, cluster_path: str, *, batch_size=1, genome_path: str, annotation_path: str, function_path: str, enable_tqdm: bool = True):

    # Print all the parameters:
    # print("Features:", features_path)
    # print("Chromosome:", chromo)
    # print("Cluster:", cluster_path)
    # print("Genome:", genome_path)
    # print("Annotation:", annotation_path)
    # print("Function:", function_path)
    # print("Batch size:", batch_size)
    # print("Enable tqdm:", enable_tqdm)

    ### General Variables.
    verbose = True

    ### Reading the annotations from the annotation.tsv file - (boilerplate code)
    _annotations = read_annotations(path=annotation_path, _chromo=chromo)

    ### Reading the genome
    genome = pd.read_csv(genome_path, sep="\t", header=None, names=["chr", "size"]).set_index("chr").to_dict()["size"]
 
    ### Importing the features
    features_ = pr.read_bed(features_path)
    features_ = features_[features_.Chromosome == chromo]  ### Filtering out by chromosome.
    features_numpy = BedShuffler.pyranges_to_numpy(features_, last_column="index")
    features_numpy[:, 1] += 1

    ### importing the clusters, making them into arrays.
    clusters, cluster_indexer = import_clusters(cluster_path)

    ### Creating the overlap counter.
    clusters_overlaps = OverlapCounter(clusters, cluster_indexer)

    ### Importing the clusters as merged.
    merged_clusters, _ = import_clusters(cluster_path, merged=True)

    ### Add a column of all zeros
    features_numpy = np.column_stack(
        [features_numpy, np.zeros(features_numpy.shape[0], dtype=np.int64)],
    )
    observed_overlaps = clusters_overlaps.count_overlaps(features_numpy)

    ### Add a column to the clusters dataframe containing the number of overlaps
    merged_clusters["observed"] = observed_overlaps

    ### Keep the clusters with at least one overlap
    #merged_clusters = merged_clusters[merged_clusters.observed > 0]
    observed_overlaps = merged_clusters.observed.values.flatten()

    ### Get the names of the clusters (these are the clusters to keep).
    cluster_names = merged_clusters.name.values.flatten()

    ### Reimport the Requested clusters.
    clusters, cluster_indexer = import_clusters(cluster_path, merged=False, keep_clusters=cluster_names)
    merged_clusters, _ = import_clusters(cluster_path, merged=True, keep_clusters=cluster_names)

    # Adding a column for the odds
    merged_clusters["odds_ratio"] = 0

    # The number of big loops we need to do.
    # So if batch_count == 4, we will do 4 loops.
    batch_count = n // batch_size

    

    # Creating the overlap counter
    overlap_counter = OverlapCounter(clusters, cluster_indexer)

    # Creating the shuffler
    shuffler = BedShuffler(BedShuffler.pyranges_to_numpy(features_, last_column=None), chromo, genome[chromo])

    # Adding all the inclusions and exclusions
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_3utr.BED"), _annotations["3utr"])
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_5utr.BED"), _annotations["5utr"])
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_down.BED"), _annotations["down"])
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_exon.BED"), _annotations["exon"])
    shuffler.add_exclude(os.path.join(function_path, chromo, f"{chromo}_inter_no.BED"), _annotations["inter_no"])

    # include intron, prom, prom12, prom23:
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_intron.BED"), _annotations["intron"])
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_prom.BED"), _annotations["prom"])
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_prom12.BED"), _annotations["prom12"])
    shuffler.add_include(os.path.join(function_path, chromo, f"{chromo}_prom23.BED"), _annotations["prom23"])

    enable_tqdm = True
    # The progress bar, it's manual as it is a nested for loop.
    if enable_tqdm:
        progress_bar = tqdm(total=n, desc=f"{red(chromo)} ", position=0)

    # Looping
    for _ in range(batch_count):
  

        ### Shuffling the features
        shuffled_entries = shuffler.shuffle(batch_size=batch_size)
        shuffled_entries[:, 1] += 1  ### Adding the features.


        ### Counting the Overlaps
        overlaps = overlap_counter.count_overlaps(shuffled_entries)
    
        ### (Updating the progress bar)
        if enable_tqdm:
            progress_bar.update(batch_size)


        ### updates the odds ratio.
        odds_ratio = calculate_odds_ratio(overlaps, observed_overlaps)


        ### Adding the odds ratio to the dataframe.
        merged_clusters["odds_ratio"] += odds_ratio

    ### Closing the progress bar.
    if enable_tqdm:
        progress_bar.close()

    ### Adding the p_value column: (unbiased p.value)
    merged_clusters["p_value"] = (merged_clusters["odds_ratio"] + 1) / (n + 1)

    return merged_clusters


if __name__ == "__main__":
    FEATURE_PATH = "../data/features/H1/CTCF/ENCFF402JJK.bed"
    GENOME_PATH = "../data/hg19.genome"
    ANNOTATION_PATH = "../data/counts/H1/CTCF/ENCFF402JJK.bed.counts.tsv"
    FUNCTION_PATH = "../data/annotations/"


    CLUSTER_PATH = "../data/clusters/H1/chr10_spec_res.json"
    p_values = compute_p_value(10000, chromo = "chr10", features_path= FEATURE_PATH, genome_path= GENOME_PATH, annotation_path= ANNOTATION_PATH, function_path= FUNCTION_PATH, cluster_path= CLUSTER_PATH)
