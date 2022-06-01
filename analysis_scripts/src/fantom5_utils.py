import pyranges as pr
import numpy as np
from src.clus_files_io import parse_clus_file

### Parses the GMT file returining a dict containing for each 
### Returns a dictionary mapping each gene entrez id to the ontology group it belongs to.
### Like this:
### {
###     "ENTREZ_ID_1": ["GO:0005737", "GO:0005739"],
###     "ENTREZ_ID_2": ["GO:0005737", "GO:0005739", "GO:0005740"],
###     ...
### }
def parse_gmt(path: str) -> dict[str, list[str]]:
    # Create a dictionary to store the results
    results = {}

    # Open the GMT file
    with open(path, "r") as f:
        # For each line in the file
        for line in f:
            # Split the line into a list of words
            words = line.split("\t")

            # The motif name
            gene_set = words[0]
            link = words[1]
            genes = words[2:]
            # strip each gene
            genes = [gene.strip() for gene in genes]

            # Add the gene set to the dictionary
            results[gene_set] = [link, genes]

    # invert the dictionary
    gmt_inverse : dict[str, list[str]] = {}

    # For each gene set
    for gene_set, (link, genes) in results.items():
        for gene in genes:
            if gene not in gmt_inverse:
                gmt_inverse[gene] = []
            gmt_inverse[gene].append(gene_set)
        
    return gmt_inverse, list(results.keys())

### read the table from the path
### Reads the table of the FANTOM5 peaks and returns a dictionary mapping each peak to its entrez id.
### Like this:
### {
###     "CAGE_peak_1": "ENTREZ_ID_1",
###     "CAGE_peak_2": "ENTREZ_ID_2",
###     ...
### }
def parse_fantom_data(path: str) -> dict[str, str]:
    # Create a dictionary to store the results
    results: dict[str, str] = {}

    header_found = False

    # Open the file
    with open(path, "r") as f:
        # For each line in the file
        for line in f:
            # Go to the next line if it starts with a #
            if line.startswith("#"):
                continue

            # deal with the header
            if not header_found:
                header_found = True
                continue

            # Split the line into a list of words
            words = line.split("\t")
            cage_id = words[0].strip()
            entrez_id = words[1].strip()

            # Add the gene to the dictionary
            results[cage_id] = entrez_id
        
    # Return the results
    return results

### Maps the given entry id to the gene entrez id
### Returns a list containing the entrez ids
def overlaps_to_entrez(overlaps: list[str], peak_to_entrez : dict[str, str]) -> list[str]:
    entrez_ids = []
    for overlap in overlaps:
        entrez_ids.append(peak_to_entrez.get(overlap, None))
    
    # Remove duplicates and NAs
    entrez_ids = [entrez_id for entrez_id in entrez_ids if entrez_id not in ["NA", None]]
    return entrez_ids

### Return the ontology group to which each entrez id belongs to.
### Maps each entrez id to the ontology group it belongs to.
### returns a list of ontology groups.
def map_entrez_to_ontology(entrez_ids: list[str], gene_ontology : dict[str, str]) -> list[str]:
    ontology_hist = {}
    unknown_entrez_ids = []
    for entrez_id in entrez_ids:
        if entrez_id not in gene_ontology:
            unknown_entrez_ids.append(entrez_id)
            continue
        ontology_groups = gene_ontology[entrez_id]

        for ontology_group in ontology_groups:
            if ontology_group not in ontology_hist:
                ontology_hist[ontology_group] = 0
            ontology_hist[ontology_group] += 1
    
    return ontology_hist, unknown_entrez_ids


### Maps the provided features (passed in as the list of names to the the given ontology groups) and returns an histogram of the counts of each group.
### Like this:
### {   
###     "GO:0005737": 1,
###     "GO:0005739": 1,
###     "GO:0005740": 3,
###     ...
### }
### The features are expected to be already void of duplicates as it would lead to a double counting.
def get_groups_for_features(features: list[str], peak_to_entrez: dict[str, str], gene_ontology : dict[str, str]) -> dict[str, int]:
    # Get the entrez ids
    entrez_ids = overlaps_to_entrez(features, peak_to_entrez)

    # Get the ontology
    groups, _ = map_entrez_to_ontology(entrez_ids, gene_ontology)
    
    # Return the groups
    return groups, len(entrez_ids)