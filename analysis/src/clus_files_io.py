### Parses a CLUS file in the form:
### chr1    cluster_1   cluster_2   ...
### chr2    cluster_1   cluster_2   ...
### ...
### And return a dictionary of the form:
### {
###     "chr1" : ["cluster_1", "cluster_2", ...],
###     "chr2" : ["cluster_1", "cluster_2", ...],
###     ...
### }
def parse_clus_file(path : str) -> dict[str, list[str]]:
    results = {}

    with open(path, "r") as f:
        for line in f:
            words = line.split("\t")
            chromosome = words[0].strip()
            clusters = [w.strip() for w in words[1:]]

            results[chromosome] = clusters

    return results