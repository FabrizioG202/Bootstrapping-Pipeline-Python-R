#%%
import sys
from matplotlib.cbook import flatten
import pyranges as pr
import numpy as np
import src.fantom5_utils as f5u
import pandas as pd

from scipy.stats import hypergeom

def features_to_entrez(features, peak_to_entrez):
    entrez_ids = []
    for feature in features:
        entrez_ids.append(peak_to_entrez.get(feature, None))

    # Remove NAs and None
    entrez_ids = [x for x in entrez_ids if x not in [None, "NA"]]

    # Remove duplicates
    entrez_ids = (set(entrez_ids))

    return entrez_ids

def read_gene_ontolgy(path : str):
    out = {}
    for line in open(path):
        line = line.strip().split("\t")
        out[line[0]] = set(line[2:])
    return out

# Todo import the bed file
if __name__ == "__main__":

    # args = sys.argv
    args = "./gene_set_enrichment.py ../analysis_results/CAGE_enriched.bed ../data/features/HMEC/CAGE/features.bed ../data/fantom5/peak_to_entrez.tsv ../data/fantom5/ontology.gmt".split(" ")

    match args:
        case [_, clusters_bed_file, features_bed_file, peak_to_entrez, ontology_path]:
            # Reading the enriched clusters file as ranges.
            clusters_ranges = pr.read_bed(clusters_bed_file)
            peak_to_entrez = f5u.parse_fantom_data(peak_to_entrez)
            features_ranges = pr.read_bed(features_bed_file)

            # Reading the ontology data.
            ontology_data = read_gene_ontolgy(ontology_path)
            
            # Getting the ids of the genes in the enriched clusters.
            features_in_clusters = features_ranges.overlap(clusters_ranges)
            features_in_clusters = set(features_in_clusters.Name.values.flatten().tolist())
            foreground_ids = features_to_entrez(features_in_clusters, peak_to_entrez)

            # Calculate the theoretical groups.
            all_features = features_ranges.Name.values.tolist()
            background_ids = features_to_entrez(all_features, peak_to_entrez)

            # the pvalues
            p_values = {}

            for set_name, set_genes in ontology_data.items():
                hits_in_sample = len(foreground_ids & set_genes)
                sample_size = len(foreground_ids)
                hits_in_population = len(background_ids & set_genes)
                fail_in_population = len(background_ids) - hits_in_population

                # Print the parameter to this point:
                print("Hits in sample:", hits_in_sample, "Sample size:", sample_size, "Hits in population:", hits_in_population, "Population size:", len(background_ids))
                p_values[set_name] = (hypergeom.sf(hits_in_sample - 1, len(background_ids), fail_in_population, sample_size), 0)


            # Make the pvalues into a dataframe
            p_values_df = pd.DataFrame(p_values).T
            p_values_df.columns = ["p_value", "odds_ratio"]

            # Sort the pvalues
            p_values_df.sort_values(by = "p_value", inplace = True, ascending = True)

        case [_, "--help" | "help"] | _:
            print("""Usage: python gene_set_enrichment.py <clusters_bed_file> <features_bed_file> <id_to_entrez_file> <ontology_path>""")
            

#%%
p_values_df.p_value.mean()