import numba
import numpy as np
import pyranges as pr
from typing import Literal
import random
import os

from utils import read_genome_file


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

    @staticmethod
    def empty(chromosome: str, chromosome_length: int) -> "BedShuffler":
        return BedShuffler(np.empty((0, 3), dtype=np.int64), "", 0)

    def load_regions(self, regions_path : str, _counts : dict[str, int]):
        """
        Loads the regions given the folder containing them and the provided counts.
        """
    
        # Adding all the inclusions and exclusions
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_3utr.BED"), _counts["3utr"])
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_5utr.BED"), _counts["5utr"])
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_down.BED"), _counts["down"])
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_exon.BED"), _counts["exon"])

        # Excluding the intergenic non-coding
        self.add_exclude(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_inter_no.BED"), _counts["inter_no"])

        # include intron, prom, prom12, prom23:
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_intron.BED"), _counts["intron"])
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_prom.BED"), _counts["prom"])
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_prom12.BED"), _counts["prom12"])
        self.add_include(os.path.join(regions_path, self.chromosome, f"{self.chromosome}_prom23.BED"), _counts["prom23"])

        
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
        total_counts = 0
        for delimiter in delimiters:
            total_counts += delimiter.count
        array_size = total_counts  # Used to be features.shape[0]
        
        shuffled_entries = np.zeros((array_size * batch_size, 4)[::-1], dtype=np.int64).T
        
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