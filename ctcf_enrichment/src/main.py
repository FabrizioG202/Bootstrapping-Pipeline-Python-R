from utils import *
import numpy as np
import ncls
from crayon import *
import pandas as pd
import numba
from PIL import Image
from tqdm.autonotebook import trange, tqdm
from shuffler import BedShuffler
import cuda_implementations
import sys


def flatten_bins(bins: list[list[int]]):
    """
    Flattens the bins array into a 1D numpy array and returns the indices at which each bin starts.
    """
    idxs = []
    out = []
    offset = 0
    for bin_ in bins:
        idxs.append(offset)
        out.extend(bin_)
        offset += len(bin_)

    # add the last index
    idxs.append(len(out))
    return np.array(out), np.array(idxs)


@numba.njit(
    boundscheck=False,
    fastmath=True,
)
def fast_bins_to_matrix_cpu(__flatten, __access, size: int) -> np.ndarray:

    out = np.zeros((size, size), dtype=np.uint16)

    for row in range(size):
        for col in range(size):
            s_1_idx = __access[row]
            e_1_idx = __access[row+1]

            s_2_idx = __access[col]
            e_2_idx = __access[col+1]

            # Count the unique values in the sorted arrays s_1 and s_2
            i = s_1_idx
            j = s_2_idx
            shared = 0
            while i < e_1_idx and j < e_2_idx:
                a = __flatten[i]
                b = __flatten[j]
                if a == b:
                    shared += 1
                    i += 1
                    j += 1
                elif a < b:
                    i += 1
                else:
                    j += 1

            out[row, col] = (e_1_idx - s_1_idx) + (e_2_idx - s_2_idx) - shared

    return out

# Handles the overlap


class OverlapCounter:
    """
    Handles the counting of the overlaps between the bin pair and the features.
    """
    bins_ncls: ncls.NCLS64
    bins_starts: np.ndarray
    resolution: int
    bin_count: int

    gpu: bool

    def __init__(self, chromosome_length: int, resolution: int, *, log_warning=True):

        # Initiating the range values.
        starts = np.arange(0, chromosome_length, resolution, dtype=np.int64)
        ends = starts + resolution - 1
        indices = np.arange(len(starts), dtype=np.int64)

        # Print the warning:
        if log_warning:
            print(
                f"Creating bins from a chromosome of length: {yellow(chromosome_length)}, with a resolution of {yellow(resolution)}. Total bins: {green(len(starts))}")

        # creating the ncls object
        self.bins_ncls = ncls.NCLS64(starts, ends, indices)
        self.bin_count = len(starts)

        # Creating attributes
        self.resolution = resolution
        self.bins_starts = starts

        # Check if there is a GPU
        self.gpu = numba.cuda.is_available()
        if self.gpu:
            print(f"GPU is available: {green('Running on GPU!')}")
        else:
            print(
                f"No GPU available: {red('Running on CPU!')}")

    @staticmethod
    def create_bins(l_idx_: np.ndarray, r_idx_: np.ndarray, bin_count: int) -> np.ndarray:
        # note to self: the lenght of the indices is always the same, independing of the resolution.
        # thus this function should not be the bottleneck and can be left in pure python.

        # Creating the bins
        out: list[list[int]] = [[] for _ in range(bin_count)]
        for l, r in zip(l_idx_, r_idx_):
            out[r].append(l)

        return out

    def count_overlaps(self, features: np.ndarray) -> np.ndarray:

        # L indices and R indices
        l_idx_, r_idx_ = self.bins_ncls.all_overlaps_both(
            features[:, 0].T, features[:, 1].T, features[:, 2].T)

        # Creating the bin_to_overlap
        bins = OverlapCounter.create_bins(l_idx_, r_idx_, self.bin_count)

        return bins

    def __get_overlap_matrix_CPU(self, features: np.ndarray) -> np.ndarray:

        # create the bins
        bins = self.count_overlaps(features)

        # flatten them and keep track of the indices
        flat_bins, idxs = flatten_bins(bins)

        # Create the matrix.
        matrix = fast_bins_to_matrix_cpu(flat_bins, idxs, self.bin_count)
        return matrix

    def __get_overlap_matrix_GPU(self, features: np.ndarray) -> np.ndarray:
        # create the bins
        bins = self.count_overlaps(features)

        # flatten them and keep track of the indices
        flat_bins, idxs = flatten_bins(bins)

        # Create the matrix.
        matrix = cuda_implementations.calculate_overlap_matrix(
            flat_bins, idxs, self.bin_count)
        return matrix

    def get_overlap_matrix(self, features: np.ndarray) -> np.ndarray:
        if self.gpu:
            return self.__get_overlap_matrix_GPU(features)
        else:
            return self.__get_overlap_matrix_CPU(features)


@numba.njit(
    boundscheck=False,
    fastmath=True,
)
def update_counts(_expected: np.ndarray, _shuffled: np.ndarray, _counts: np.ndarray, bin_count: int):
    for i in range(bin_count):
        for j in range(bin_count):
            if _shuffled[i, j] >= _expected[i, j]:
                _counts[i, j] += 1


def get_enrichment_matrix(
    chromo: str,
    resolution: str | int,
    features: np.ndarray,
    counts_: dict[str, int],
    genome_file: str,
    region_folder: str,
    N: int
) -> dict[int, int]:

    # Reading the genome
    genome = read_genome_file(genome_file)

    # Converting the resolution
    resolution = resolution if isinstance(
        resolution, int) else resolution_string_to_int(resolution)

    # Creating the shuffler
    shuffler: BedShuffler = BedShuffler(features, chromo, genome[chromo])

    # adding the counts and all.
    shuffler.load_regions(region_folder,
                          counts_)  # Load the regions

    # Creating the overlap counter
    overlap_counter = OverlapCounter(
        genome[chromo], resolution, log_warning=True)

    # Make the features arrays f contiguous
    features = np.asfortranarray(features)

    # Count the observed overlaps
    observed_matrix = overlap_counter.get_overlap_matrix(features)

    # create the count matrix
    count_matrix = np.zeros(
        (overlap_counter.bin_count, overlap_counter.bin_count), dtype=np.uint16)

    ### Shuffling and counting
    for _ in trange(N):
        # Shuffle the features
        features = shuffler.shuffle(batch_size=1)

        # Count the random overlaps
        random_matrix = overlap_counter.get_overlap_matrix(features)

        # Update the counts
        update_counts(observed_matrix, random_matrix,
                      count_matrix, overlap_counter.bin_count)

    # Fast-compute the pvalue
    pvalues_matrix = count_matrix / N

    return pvalues_matrix


if __name__ == "__main__":
    args = sys.argv[1:]

    match args:
        case [chromo, resolution, features_path, counts_path, genome_path, regions_path, N, out_path]:
            # dealing with N
            # if N is only numerical, convert it to int
            if resolution.isdigit():
                resolution = int(resolution)
            N = int(N)

            features = read_bed_as_array(features_path, chromo)
            counts_: dict[str, int] = read_counts_file(
                counts_path, chromo)

            matrix = get_enrichment_matrix(
                chromo, resolution, features, counts_, genome_path, regions_path, N)

            np.save(out_path, matrix)
            print(green(f"Saved the matrix to {out_path}"))
            

        case ["--help" | "--h" | "-h"] | _:
            print(
                "Usage: python main.py <chromo> <resolution> <features_path> <counts_path> <genome_path> <regions_path> <N> <out_path>")
