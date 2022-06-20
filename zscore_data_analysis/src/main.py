import tqdm
from unicodedata import numeric
import numpy as np
import numba
import os
from crayon import *
import sys


def read_genome_file(file_name: str) -> dict[str, int]:
    """
    Reads a genome file and returns a dictionary with the genome name and the length of the genome.
    """
    genome_dict = {}
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip()
            chromo, length = line.split('\t')
            genome_dict[chromo.strip()] = int(length)
    return genome_dict


class DataSplitter:

    # Hubs contains the bed file as a list of intervals like this:
    # [start, end]
    # [start, end]
    # ...
    # [start, end]
    hubs: np.ndarray

    # The parent chromosome
    chromo: str

    def __parse_bed_as_array(self, bed_path: str) -> np.ndarray:
        """
        Parses the bed file and returns an array of intervals.
        """
        intervals = []
        with open(bed_path, "r") as bed_file:
            for line in bed_file:
                line = line.strip()
                if line == "":
                    continue

                # Split the line into a list of values
                chr_, start, end, *rest = line.split("\t")

                # Check if the chromosome is the same as the parent chromosome
                if chr_ != self.chromo:
                    continue

                # Convert the start and end to int
                start = int(start)
                end = int(end)

                # Add the interval to the list
                intervals.append((start, end))

        return np.array(intervals)

    def __init__(self, chromo: str, hubs_path: str, resolution: int):
        self.chromo = chromo
        self.resolution = resolution
        self.hubs = self.__parse_bed_as_array(hubs_path)

    @staticmethod
    @numba.njit
    def check_single_pair_nb(hubs: np.ndarray, pair: tuple[int, int], resolution: int):
        """
        Checks wether the given pair is a valid pair (belongs inside the hubs).
        """
        for hub in hubs:
            contained = True
            for p in pair:

                start = p
                end = p + resolution

                if start < hub[0] or end > hub[1]:
                    contained = False
                    break

            # If the pair is inside the hub, return True
            if contained:
                return True

        return False

    def check_single_pair(self, pair: tuple[int, int]):
        return DataSplitter.check_single_pair_nb(self.hubs, pair, self.resolution)

    def split_data(self, tsv_path: str, matrix_path: str, out_folder: str):
        """
        Splits the matrix data into two files:
        - one containing the pairs that are inside the hubs
        - one containing the pairs that are outside the hubs
        """

        # the file paths.
        inside_file_path = out_folder + f"{self.chromo}_inside.tsv"
        outside_file_path = out_folder + f"{self.chromo}_outside.tsv"
        # create the output folder if it doesn't exist
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)

        # The inside counts
        inside_count = 0
        outside_count = 0

        # The matrices.
        matrix = np.load(matrix_path)

        # Opening the files.
        with open(inside_file_path, "w") as inside_file, open(outside_file_path, "w") as outside_file, open(tsv_path, "r") as tsv_file:
            # Add amn header to the inside and outside files
            header = "bin1\tbin2\tzscore\tpvalue\n"
            inside_file.write(header)
            outside_file.write(header)

            # Looping over the lines.
            for line in tsv_file:
                line = line.strip()
                if line == "":
                    continue

                # Split the line into a list of values
                start1, start2, *rest = line.split("\t")

                # Convert the start and end to int
                start1 = int(start1)
                start2 = int(start2)

                # get the index coefficients
                i = start1 // self.resolution
                j = start2 // self.resolution

                # retrieve the p-value from the matrix
                pvalue = matrix[i, j]

                line = line + "\t" + str(pvalue)

                # Check if the pair is inside the hubs
                if self.check_single_pair((start1, start2)):
                    inside_file.write(line + "\n")
                    inside_count += 1
                else:
                    outside_file.write(line + "\n")
                    outside_count += 1
        #print(f"Split the data between {green(inside_count)} inside and {red(outside_count)} outside pairs.")


def join_all_files(out_folder: str):
    """
    Joins all the files in the given folder.
    """
    files = os.listdir(out_folder)

    # sort the files by chromosome:
    files.sort(key=lambda x: int(x.split("_")[0].replace("chr", "")))
    inside_file = os.path.join(out_folder, "inside.tsv")
    outside_file = os.path.join(out_folder, "outside.tsv")

    #Create the progress bar
    pbar = tqdm.tqdm(total=len(files))

    # Open the files and write stuff to them.
    with open(inside_file, "w") as inside_file, open(outside_file, "w") as outside_file:
        for file in files:
            
            # Set the description of the progress bar
            pbar.set_description(f"Joining {yellow(file)}")

            # Open the file
            if file.endswith("inside.tsv"):
                with open(out_folder + file, "r") as f:
                    # skip the header
                    next(f)
                    for line in f:
                        inside_file.write(line)
            elif file.endswith("outside.tsv"):
                with open(out_folder + file, "r") as f:
                    # skip the header
                    next(f)
                    for line in f:
                        outside_file.write(line)

            # Update the Progress Bar
            pbar.update(1)

    # Remove the files
    for file in files:
        os.remove(out_folder + file)



# Smart parse a boolean string
def smart_bool(str: str):
    return str.lower() in ("yes", "true", "t", "1")


AUTOSOMAL_CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
ALL_CHROMOSOMES = AUTOSOMAL_CHROMOSOMES + ["chrX", "chrY"]


def smart_resolution(value: str):
    # check if the vlaue is a number
    if value.isdigit():
        return int(value)

    else:
        prefixes = ["kb", "Mb", "Gb", "Tb"]
        for prefix in prefixes:
            if value.lower().endswith(prefix.lower()):
                return int(value[:-len(prefix)]) * 1000**prefixes.index(prefix)


if __name__ == "__main__":
    # Read the command line arguments
    argv = sys.argv[1:]

    match argv:

        case [enrichment_path, zscores_path, hubs_path, output_folder, resolution, *only_autosomal]:

            # Parse the parameters
            only_autosomal = smart_bool(only_autosomal[0]) if only_autosomal else False
            resolution = smart_resolution(resolution)

            # Define the looping variables
            loop_over = AUTOSOMAL_CHROMOSOMES if only_autosomal else ALL_CHROMOSOMES
            pbar = tqdm.tqdm(total=len(loop_over), desc="Splitting data: ")

            # Loop over the chromosomes.
            for chromo in loop_over:

                # Set the name of the chromosome as the progress bar
                pbar.set_description(f"Splitting data: {yellow(chromo)}")

                # Creating the splitter
                splitter = DataSplitter(chromo, hubs_path, resolution=100_000)

                # Splitting the data
                splitter.split_data(
                    matrix_path=f"{enrichment_path}/{chromo}.npy",
                    tsv_path=f"{zscores_path}/{chromo}.txt",
                    out_folder=output_folder
                )
                pbar.update(1)

            # Close the progress bar
            pbar.close()

            join_all_files("../results/hubs_data/100kb/")
            print(green("Successfully joined all the files and removed the temporary files."))

        case _:
            print(red("Usage: python3 split_data.py enrichment_path zscores_path hubs_path output_folder resolution [only_autosomal]"))
            print(red("Only_autosomal: True or False"))

