import sys
# supported operartions: 
# - tsv to clus file
# - clus file to tsv
# - clus file to bed.
import sys

def replace_extension(path : str, new_extension : str) -> str:
    parts = path.split(".")
    parts[-1] = new_extension
    return ".".join(parts)

import src.clus_files_io, src.cluster_description

# usage script conversion_type
if __name__ == "__main__":

    match sys.argv:
        case [_, "tsv2clus", source_path, *args]:

            # the output path, if not provided, just replace the the extension.
            out_path = (args[0]) if len(args) > 0 else replace_extension(source_path, "clus")
            chromo_column = int(args[1]) if len(args) > 0 else 0
            clus_column = int(args[2]) if len(args) > 1 else 2

            src.clus_files_io.tsv_to_clus_file(source_path, out_path, chromo_column_idx=chromo_column, name_column_idx = clus_column)

        case [_, "clus2bed", clusters_folder, source_path, *args]:
            print("Converting clus file to bed file.")

            # read the clusters from the file.
            clusters : dict[str, list[str]] = src.clus_files_io.parse_clus_file(source_path)
            import os

            out = ""
            for chromo, cluster_ids in clusters.items():
                description = src.cluster_description.ClustersDescription(os.path.join(clusters_folder, chromo + "_spec_res.json"), chromo)

                for _idx in cluster_ids:
                    cluster = description[_idx]
                    for start, end in cluster:
                        out += "{}\t{}\t{}\t{}\n".format(chromo, start, end, _idx)

            out_path = args[0] if len(args) > 0 else replace_extension(source_path, "bed")

            # save the content of out to a file
            with open(out_path, "w") as f:
                f.write(out)
            
            # Print that the file was created.
            print("Done.")

        case [_, "help"] | _:
            #Print the help:
            print("Usage:")
            print("\tpython3 clus_io.py tsv2clus <source_path> <out_path> <chromo_column> <clus_column>")
            print("\tpython3 clus_io.py clus2bed <clusters_folder> <source_path> <out_path>")