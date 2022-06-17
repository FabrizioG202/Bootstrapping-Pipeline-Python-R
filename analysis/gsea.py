from scipy.stats import hypergeom
from tqdm import tqdm
import sys
from src.crayon import *


def read_gene_list(path):
    """
    Reads a list of genes from a file.
    """
    with open(path, "r") as f:
        genes = f.read().splitlines()
    genes = [gene.strip().upper() for gene in genes]
    return set(genes)


def read_gene_ontology(path):
    """
    Reads a list of genes from a file.
    """

    out = dict()
    with open(path, "r") as f:
        for line in f:
            if line.startswith("!"):
                continue

            _set, _link, *genes = line.split("\t")
            genes = set([gene.strip().upper() for gene in genes])
            out[_set.strip()] = genes

    return out


if __name__ == "__main__":
    argv = sys.argv[1:]

    match argv:
        case [foreground_path, background_path, ontology_path, output_path]:

            # Read the foreground genes
            foreground_genes = read_gene_list(foreground_path)

            # Read the background genes
            background_genes = read_gene_list(background_path)
            # Reading the ontology
            ontology = read_gene_ontology(ontology_path)

            enrichment = {}
            for set_name, genes in tqdm(ontology.items(), desc="Enrichment"):
                hits_in_sample = len(foreground_genes.intersection(genes))
                sample_size = len(foreground_genes)
                hits_in_background = len(background_genes.intersection(genes))
                background_size = len(background_genes)
                fails_in_background = background_size - hits_in_background
                p_value = hypergeom.sf(hits_in_sample - 1, hits_in_background + fails_in_background, hits_in_background, sample_size)

                # If this gene is not in the background, it will not be in the sample either.
                if hits_in_background == 0:
                    continue

                # odds_ratio <-(hitInSample/sampleSize)/(hitInPop/length(cage_active_genes_vec))
                odds_ratio = (hits_in_sample / sample_size) / (hits_in_background / background_size)

                # Adding the value
                enrichment[set_name] = (p_value, odds_ratio, hits_in_sample)

            with open(output_path, "w") as f:
                # Add an header
                f.write("Set\tp-value\todds-ratio\tsample\n")
                for name, (p_value, odds_ratio, hits_in_sample) in sorted(enrichment.items(), key=lambda x: x[1][0]):
                    f.write("{}\t{}\t{}\t{}\n".format(name, p_value, odds_ratio, hits_in_sample))

            print(green("Saved to {}".format(output_path)))


        case _:
            print("Usage: python3 gsea.py <foreground_path> <background_path> <ontology_path> <output_path>")
