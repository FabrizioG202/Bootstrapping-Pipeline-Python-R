#%%
import os, json
import pandas as pd
import enum 
import statsmodels.stats.multitest as sm
from src.utils import sort_by_chromosome, RESOLUTIONS
from src .cluster_description import Cluster, ClustersDescription

### The values representing the columns in the results files
class ValueKeys(enum.Enum):
    CHROMOSOME = "chromosome"
    NAME = "name"
    P_VALUE = "p_value"
    Q_VALUE = "q_value"

### Whether to apply FDR for each chromosome separately or to the whole genome at the same time.
class FDRMode(enum.Enum):
    PER_CHROMOSOME = "per_chromosome"
    WHOLE_GENOME = "whole_genome"
    DAGGER = "dagger"

### Which plotting library to use (not used right now)
class PlotLibrary(enum.Enum):
    SEABORN = "seaborn"
    MATPLOTLIB = "matplotlib"
    PLOTLY = "plotly"

### How to serialize the BED file representing the clusters.
### If "multiple_clusters" is selected, a single track will contain multiple clusters.
### If "single_cluster" is selected, a single track will contain a single cluster and thus the file will contain multiple tracks.
class BedSerializationMode(enum.Enum):
    MULTIPLE_CLUSTERS = "multiple_clusters"
    SINGLE_CLUSTER = "single_cluster"

### A single results table.
class ResultsTable:
    path : str
    df : pd.DataFrame

    def __init__(self, path, df):
        self.path = path
        self.df = df

    @staticmethod 
    def from_file(path) -> "ResultsTable":
        df = pd.read_csv(path, sep="\t")
        # add a column for the resolution
        df["resolution"] = df["name"].apply(lambda x: RESOLUTIONS[x.split("_")[0]])
        return ResultsTable(path, df)

    # Utilities for querying the results and filtering them.
    def where(self, key : str, *, less_than=None, greater_than=None, is_in=None,equals=None, not_equals=None, ) -> "ResultsTable":
        _df = self.df.copy()
        if less_than is not None:
            _df = _df[_df[key] < less_than]
        if greater_than is not None:
            _df = _df[_df[key] > greater_than]
        if is_in is not None:
            _df = _df[_df[key].isin(is_in)]
        if equals is not None:
            _df = _df[_df[key] == equals]
        if not_equals is not None:
            _df = _df[_df[key] != not_equals]
        return ResultsTable(self.path, _df)

    @property
    def names(self):
        return self.df["name"].tolist()

    @property
    def size(self):
        return len(self.df)

    ### Applies FDR
    def apply_fdr(self, mode : FDRMode) -> "ResultsTable":

        # FDR applied using DAGGER.
        if mode == FDRMode.DAGGER:
            raise NotImplementedError("DAGGER-based FDR not implemented yet!")

        # FDR applied on a single chromosome
        if mode == FDRMode.PER_CHROMOSOME:

            # Get how many chromosomes there are
            chromosomes = self.df["chromosome"].unique()
            _df = self.df.copy()

            # loop over them correcting for FDR
            for chromosome in chromosomes:
                _df.loc[_df["chromosome"] == chromosome, "q_value"] = sm.fdrcorrection(_df.loc[_df["chromosome"] == chromosome, "p_value"], alpha=0.05)[1]
            
            # return the new table
            return ResultsTable(self.path, _df)
        
        # Do the same for the whole genome
        if mode == FDRMode.WHOLE_GENOME:
            _df = self.df.copy()
            _df["q_value"] = sm.fdrcorrection(_df["p_value"], alpha=0.05)[1]
            return ResultsTable(self.path, _df)
    
    ### plots the dataframe:
    def plot_chromosome_counts(self, key : str, plot_library : PlotLibrary, *, title=None, xlabel=None, ylabel=None, legend=True, **kwargs):
        if plot_library == PlotLibrary.PLOTLY:
            import plotly.graph_objs as go
            import plotly.express as px
            
            # plots an histogram of the entries for each chromosome
            fig = go.Figure()
            counts = self.df.groupby("chromosome")[key].count().to_dict()
            counts : dict = sort_by_chromosome(counts)
            fig.add_trace(go.Bar(x=list(counts.keys()), y=list(counts.values()), name="Counts"))
            fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel, legend_title_text="Chromosome")
            return fig

        else:
            raise NotImplementedError("Plotting library not implemented")

    ### Writes the bed to a file.
    def write_bed(self, path  :str, clusters : ClustersDescription, mode : BedSerializationMode = BedSerializationMode.MULTIPLE_CLUSTERS):
        
        if mode == BedSerializationMode.MULTIPLE_CLUSTERS:
            chromosomes = []
            starts = []
            ends = []
            names = []

            for row in self.df.itertuples():
                chromo = row.chromosome
                assert chromo == clusters.chromosome, "Chromosome mismatch, consider filtering before"

                cluster = clusters[row.name]
                bins = cluster.bins
                _ends = [s + int(cluster.resolution) for s in bins]

                chromosomes.extend([chromo] * len(bins))
                starts.extend(bins)
                ends.extend(_ends)
                names.extend([row.name] * len(bins))

            __df = pd.DataFrame({"chromosome": chromosomes, "start": starts, "end": ends, "name": names})
            __df.to_csv(path, sep="\t", index=False, header=False)

        else:
            with open(path, "w") as f:
                for row in self.df.itertuples():
                    chromo = row.chromosome
                    assert chromo == clusters.chromosome, "Chromosome mismatch, consider filtering before"
                    cluster = clusters[row.name]
                    bins = cluster.bins
                    _ends = [s + int(cluster.resolution) for s in bins]

                    
                    # Write the track name and the description, which contains the name of the cluster and the resolution
                    f.write(f"track name={row.name} description=\"{row.name}_{cluster.resolution}\"\n")

                    for s, e, n in zip(bins, _ends, [row.name] * len(bins)):
                        f.write(f"{chromo}\t{s}\t{e}\n")
                        
    ### Intersect 2 results tables.
    def intersect(self, other : "ResultsTable") -> "ResultsTable":
        return ResultsTable(self.path, self.df.merge(other.df, on=["name", "chromosome"], ))

    ### Writes the table to a CLUS file like this:
    ### chr1   cluster1_name  cluster2_name cluster3_name ...
    ### chr2   cluster1_name  cluster2_name cluster3_name ... 
    def write_clus_file(self, path : str):
        chromosomes = self.df["chromosome"].unique()
        TAB = "\t" #as a tab cannot be used in f strings.
        with open(path, "w") as f:
            for chromosome in chromosomes:
                _df = self.df.loc[self.df["chromosome"] == chromosome]
                f.write(f"{chromosome}\t{TAB.join(_df['name'].tolist())}\n")


    # for jupyter notebooks
    def __repr__(self):
        return self.df.to_string()
    
    def _repr_html_(self):
        return self.df.to_html()
            

    def __len__(self):
        return len(self.df)

    # Save the table to a file
    def save(self, path : str):
        # do not save the resolution column
        self.df[self.df.columns.difference(["resolution"])].to_csv(path, sep="\t", index=False)

#%%
import sys

def replace_extension(path, extension : str):
    parts = path.split(".")
    parts[-1] = extension
    return ".".join(parts)

if __name__ == "__main__":

    match sys.argv:
        case [_, "applyFDR", file_path, fdr_mode, *args]:
            table = ResultsTable.from_file(file_path)
            _fdr_mode = None
            match fdr_mode:
                case "dagger":
                    _fdr_mode = FDRMode.DAGGER
                case "chr" | "per_chromosome" | "perChromosome":
                    _fdr_mode = FDRMode.PER_CHROMOSOME
                case "genome" | "whole_genome" | "wholeGenome":
                    _fdr_mode = FDRMode.WHOLE_GENOME
            
            # apply FDR
            table = table.apply_fdr(_fdr_mode)

            # save_path 
            save_path = args[0] if len(args) != 0 else file_path.replace(".tsv", f"_FDR_{fdr_mode}.tsv")
            table.save(save_path)

        case [_, "filter", file_path, parameter, action, cutoff, *args]:
            table = ResultsTable.from_file(file_path)
            _parameter = None
            cut_off = float(cutoff)
            
            match parameter:
                case "p_value" | "p" | "pValue":
                    _parameter = ValueKeys.P_VALUE.value
                case "q_value" | "q" | "qValue":
                    _parameter = ValueKeys.Q_VALUE.value
            
            # filter
            match action:
                case "less_than" | "lt":
                    table = table.where(_parameter, less_than=cut_off)
                case "greater_than" | "gt":
                    table = table.where(_parameter, greater_than=cut_off)
                case "equal_to" | "eq":
                    table = table.where(_parameter, equal_to=cut_off)

            # save_path
            save_path = args[0] if len(args) != 0 else file_path.replace(".tsv", f"_filtered_{parameter}_{action}_{cut_off}.tsv")
            table.save(save_path)

        case [_, "toClus" | "2Clus", file_path, *args]:
            table = ResultsTable.from_file(file_path)
            save_path = args[0] if len(args) != 0 else replace_extension(file_path, "clus")
            table.write_clus_file(save_path)

        case [_, "--help" | "--h" | "-h" | "-help"]:
            print("""
            Usage: python ResultsTable.py [OPTIONS]
            Options:
                applyFDR <file_path> <fdr_mode> [save_path]
                filter <file_path> <parameter> <action> <cutoff> [save_path]
                --help | -h | -help
            """)
            exit(0)