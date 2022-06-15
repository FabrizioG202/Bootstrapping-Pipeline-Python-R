import sys
from src.crayon import *

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
            data = data.merge()
            
            # Overwrite the source file.
            data.to_bed(source_path)

            print(green("Successfully collapsed, saved to: " + source_path))
        case [_, "help"] | _:
            print("Usage:")
            print("\tpython3 bed_utils.py collapse/merge <source_path>")
            print("\tpython3 bed_utils.py intersect/overlap <other_path> <source_path> <out_path>")