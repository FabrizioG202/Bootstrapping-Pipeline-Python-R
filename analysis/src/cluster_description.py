# %%
# The table for converting resolutions to string.
RESOLUTIONS = {
    "5kb": 5e3,
    "10kb": 1e4,
    "50kb": 5e4,
    "100kb": 1e5,
    "500kb": 5e5,
    "1Mb": 1e6,
}

# Extracts bins from the strange json storage method.
def extract_bins( x: str | list) -> list[int]:
    if type(x) == str:
        return [int(x)]
    elif type(x) == list:
        return [int(i) for i in x]
    else:
        raise TypeError(f"The type of the argument is not a string or a list. It is: {type(x)}")

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
    
import json
import pyranges as pr

class Cluster:
    name : str
    resolution : str
    bins : str

    def __init__(self, name : str, bins : list[int]):
        self.name = name
        self.bins = bins
        self.resolution = RESOLUTIONS[self.name.split("_")[0]]

    def __str__(self):
        return f"{self.name}, resolution: {self.resolution}, bins: {self.bins}"

    def __repr__(self):
        return f"{self.name}, resolution: {self.resolution}, bins: {self.bins}"

    # provides an iterator over the bins starts and ends
    def __iter__(self):
        return iter(zip(self.bins, [x + int(self.resolution) for x in self.bins]))

    # get the total length
    @property
    def total_length(self):
        return self.resolution * len(self.bins)

    def __len__(self):
        return len(self.bins)

    def find_overlaps(self, chromosome : str, features : pr.PyRanges):

        # build the pyranges for the cluster
        starts, ends = zip(*self)
        cluster_range = pr.PyRanges(chromosomes = [chromosome] * len(starts), starts = starts, ends = ends)

        # Check the overlaps:
        overlaps = features.overlap(cluster_range)
        
        return overlaps

# %%
class ClustersDescription:
    # The path to the clusters description file
    path : str 

    # The bins:
    members : dict[str, Cluster]

    # The chromosome to which this description belongs.
    chromosome : str

    # The tree
    def __init__(self, path : str, chromosome : str):
        self.path = path
        self.chromosome = chromosome
        self.members = {}
        with open(self.path) as f:
            clusters = json.load(f)
            cl_members = clusters["cl_member"]
            for cluster, bins in cl_members.items():

                self.members[cluster] = Cluster(cluster, extract_bins(bins))

    # gets the item by name
    def __getitem__(self, key):
        return self.members[key]

    @property
    def size (self):
        return len(self.members)
    
    def __str__(self):
        return f"Clusters description: {self.path}"

    def bins_of(self, cluster : str):
        return self.members[cluster].bins

