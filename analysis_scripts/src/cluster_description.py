# %%
# The table for converting resolutions to string.
try: from utils import RESOLUTIONS, extract_bins    
except: from .utils import RESOLUTIONS, extract_bins
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

