#%%
from matplotlib.cbook import flatten


def intersect_clusters(clus_files  : list[dict[str, list[str]]]) -> dict[str, list[str]]:
    """
    Intersects two clusters.
    """
    # the output.
    out = {}
    
    # the chromosomes are just the set of flattened lists of keys from the dict.
    chromosomes = set(flatten([list(d.keys()) for d in clus_files]))
    
    # for each chromosome.
    for chromo in chromosomes:
        # the clusters for this chromosome.
        clusters = [d[chromo] for d in clus_files if chromo in d] 
        
        # the intersection of the clusters.
        out[chromo] = list(set.intersection(*map(set, clusters)))

    # remove empty clusters.
    out = {k:v for k,v in out.items() if len(v) > 0}
        
    return out

def save_clus_file(path : str, clusters : dict[str, list[str]]) -> None:
    """
    Saves a clus file.
    """
    with open(path, "w") as f:
        for chromo, cluster_ids in clusters.items():
            f.write("{}\t{}\n".format(chromo, "\t".join(cluster_ids)))

def parse_clus_file(path : str) -> dict[str, list[str]]:
    results = {}

    with open(path, "r") as f:
        for line in f:
            # Skip comments (#)
            if line.startswith("#"): continue

            # Read the rest.
            words = line.split("\t")
            chromosome = words[0].strip()
            clusters = [w.strip() for w in words[1:]]

            results[chromosome] = clusters

    return results

#%%
import sys
if __name__ == "__main__":

    match sys.argv:
        case [_, "intersectClus", out_path, *clus_files]:
            clus_files = [d for d in clus_files if d != ""]
            clusters = intersect_clusters([parse_clus_file(f) for f in clus_files])
            save_clus_file(out_path, clusters)

        case [_, "help"] | _:
            print("""
            Usage:
                clus_utils.py [intersectClus] [out_path] [clus_files...]
                clus_utils.py help
            """)

        
        
