import sys
from src.crayon import *

if __name__ == "__main__":
    argv = sys.argv

    match argv:
        
        case [_, "intersect" | "overlap", source_path, other_path, out_path, *args]:
            import pyranges as pr
            data = pr.read_bed(source_path)
            other = pr.read_bed(other_path)

            data = data.overlap(other)

            # save the data
            data.to_bed(out_path)

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
            print("\tpython3 bed_utils.py collapse <source_path>")
            print("\tpython3 bed_utils.py keepIn <other_path> <source_path> <out_path>")