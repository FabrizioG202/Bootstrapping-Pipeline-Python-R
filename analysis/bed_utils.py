import sys
from src.crayon import *

AUTOSOMAL_CHROMOS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

if __name__ == "__main__":
    argv = sys.argv

    match argv:
        
        case [_, "intersect" | "overlap", source_path, other_path, out_path, *args]:
            import pyranges as pr
            data = pr.read_bed(source_path)
            data_length_before = len(data)
            
            # The ranges defining the region in which we want the features to be contained.
            other = pr.read_bed(other_path)
            
            # Overlapping the data.
            data = data.overlap(other)

            # save the data
            data.to_bed(out_path)

            data_length_after = len(data)
            
            # Print diagnostics
            retained_percentage = data_length_after / data_length_before * 100
            print(f"Number of retained features: {green(data_length_after)}/{yellow(data_length_before)} ({italic(int(retained_percentage))}%)")

            # Print the number of features that were removed.
            print(green("Successfully intersected, saved to: " + out_path))
        
        case [_, "collapse" | "merge", source_path, *args]:
            import pyranges as pr
            data = pr.read_bed(source_path)
            data_length_before = len(data)

            data = data.merge()

            data_length_after = len(data)

            # Print diagnostics
            retained_percentage = data_length_after / data_length_before * 100
            print(f"Reduced number of features from {yellow(data_length_before)} to {green(data_length_after)} ({italic(int(retained_percentage))}%)")
            # Overwrite the source file.
            data.to_bed(source_path)

            print(green("Successfully collapsed, saved to: " + source_path))

        
        case [_, "keepAutosomal" | "autosomal", source_path, *args]:

            # get the outpath, if none is provided, replace
            out_path = args[0] if len(args) > 0 else source_path

            length_before = 0
            length_after = 0

            out = ""
            with open(source_path, "r") as source_file:
                for line in source_file.readlines():
                    
                    # Skip comments
                    if line.startswith("#"):
                        out += line
                        continue
                    
                    # Read the chromosome
                    chromo = line.split("\t")[0]
                    
                    # Check that the chromosome is autosomal.
                    if not chromo.startswith("chr"):
                        raise ValueError(f"Chromosome {chromo} is not in the expected format.")

                    if chromo in AUTOSOMAL_CHROMOS:
                        out += line
                        length_after += 1

                    length_before += 1

            # Rewrite the file.
            with open(out_path, "w") as source_file:
                source_file.write(out)
            
            # Print diagnostics
            retained_percentage = length_after / length_before * 100
            print(f"Number of retained features: {green(length_after)}/{yellow(length_before)} ({italic(int(retained_percentage))}%)")

            # Print the number of features that were removed.
            print(green("Successfully kept autosomal features, saved to: " + source_path))

        case [_, "help"] | _:
            print("Usage:")
            print("\tpython3 bed_utils.py collapse/merge <source_path>")
            print("\tpython3 bed_utils.py intersect/overlap <other_path> <source_path> <out_path>")
            print("\tpython3 bed_utils.py keepAutosomal <source_path> <out_path>")

