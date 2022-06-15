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
            # Skip comments (#)
            if line.startswith("#"): continue

            # Read the rest.
            words = line.split("\t")
            chromosome = words[0].strip()
            clusters = [w.strip() for w in words[1:]]

            results[chromosome] = clusters

    return results

def clus_file_to_tsv(path : str, out_path : str) -> None:
    clus_dict = parse_clus_file(path)

    with open(out_path, "w") as f:
        for chromosome in clus_dict:
            for cluster in clus_dict[chromosome]:
                f.write("{}\t{}\n".format(cluster, chromosome))
            

def tsv_to_clus_file(path : str, out_path : str, chromo_column_idx : int = 0, name_column_idx : int = 2, has_header  : bool = True) -> None:
    clus_dict = {}

    with open(path, "r") as f:
        if has_header: next(f) # Skip header.
        for line in f.readlines():
            words = line.split("\t")
            chromosome = words[chromo_column_idx].strip()
            cluster = words[name_column_idx].strip()

            if chromosome not in clus_dict:
                clus_dict[chromosome] = []

            clus_dict[chromosome].append(cluster)

    with open(out_path, "w") as f:
        for chromosome in clus_dict:
            f.write("{}\t{}\n".format(chromosome, "\t".join(clus_dict[chromosome])))

    