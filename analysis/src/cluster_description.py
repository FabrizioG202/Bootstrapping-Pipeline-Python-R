# %%
import src.utils
import json

class Cluster:
    name : str
    resolution : str
    bins : str

    def __init__(self, name : str, bins : list[int]):
        self.name = name
        self.bins = bins
        self.resolution =src.utils.RESOLUTIONS[self.name.split("_")[0]]

    def __str__(self):
        return f"{self.name}, resolution: {self.resolution}, bins: {self.bins}"

    def __repr__(self):
        return f"{self.name}, resolution: {self.resolution}, bins: {self.bins}"


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

                self.members[cluster] = Cluster(cluster, src.utils.extract_bins(bins))

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

