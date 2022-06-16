from src.crayon import *
import sys
from src.fantom5_io import parse_fantom_peak_type_data

def replace_extension(path : str, new_extension : str) -> str:
    parts = path.split(".")
    parts[-1] = new_extension
    return ".".join(parts)

def extract_ranges_from_name(name :str) -> tuple[str, int, int]:

    name = name.split(",")[0]
    chromo, range = name.split(":")

    # Split either by - or ..
    if ".." in range:
        start, end = range.split("..")
    else:
        start, end = range.split("-")
    
    return chromo, int(start), int(end)


if __name__ == "__main__":
    argv = sys.argv

    match argv:
        
        case [_, "filter", source_bed_path, fantom_5_path, *args]:
            # read the peak type as a dict.
            peak_types = parse_fantom_peak_type_data(fantom_5_path)
            
            out_file_path = args[0] if len(args) > 0 else replace_extension(source_bed_path, "fantom5.types")
            
            with open(source_bed_path, "r") as source_bed_file, open(out_file_path, "w") as out_file:
                
                # Looping over the lines
                for i, line in enumerate(source_bed_file.readlines()):
                    
                    # skip headers
                    if line.startswith("#"):
                        continue
                
                    # split the line and get the peak name.
                    try:
                        peak_name = line.split("\t")[3].strip()
                    except IndexError:
                        print("Error: line {} is not in the expected format.".format(i))
                        print(source_bed_path)
                    peak_type = peak_types.get(peak_name, None)
                    if peak_type is None:
                        print(peak_name)
                        continue

                    out_file.write(f"{peak_name}\t{peak_types[peak_name]}\n")

            print(green("Successfully filtered the peaks and wrote to:" +  out_file_path))


        case [_, "getRanges", "enhancer" | "tss" as peak_type, source_types_path, out_path, *args]:

            if peak_type == "tss":
                if len(args) == 0:
                    print(red("Please provide a radius for the TSS extraction!"))
                    exit(1)
                else:
                    tss_radius = int(args[0])
            
            count = 0
            
            with open(source_types_path, "r") as source_types_file, open(out_path, "w") as out_file:
                for i, line in enumerate(source_types_file.readlines()):
                    if line.startswith("#"):
                        continue
                    
                    peak_name, pt = line.split("\t")
                    peak_name = peak_name.strip()
                    pt = pt.strip()

                    if peak_type != pt:
                        continue

                    chromo, start, end = extract_ranges_from_name(peak_name)
                    
                    if peak_type == "tss":
                        start = start - tss_radius
                        end = end + tss_radius

                    out_file.write(f"{chromo}\t{start}\t{end}\n")
                    count += 1

            
            print(green(f"Successfully extracted {count} ranges to: " +  out_path))

        case [_, "countTypes", source_types_path, *args]:
            file = open(source_types_path, "r")
            enhancer_count = 0
            tss_count = 0
            for line in file.readlines():
                if line.startswith("#"):
                    continue
                
                peak_name, pt = line.split("\t")
                peak_name = peak_name.strip()
                pt = pt.strip()

                if pt == "enhancer":
                    enhancer_count += 1
                elif pt == "tss":
                    tss_count += 1

                
            print(f"Enhancers: {blue(enhancer_count)}, TSS: {green(tss_count)}")

        case [_, "background", "enhancer" | "tss" | "all" as peak_type, source_types_path, out_path]:
            
            file = open(source_types_path, "r")
            out_file = open(out_path, "w")

            with open(source_types_path, "r") as source_types_file, open(out_path, "w") as out_file:
                for i, line in enumerate(source_types_file.readlines()):
                    if line.startswith("#"):
                        continue
                    
                    peak_name,  pt = line.split("\t")
                    peak_name = peak_name.strip()
                    pt = pt.strip()

                    if peak_type != "all" and peak_type != pt:
                        continue

                    chromo, start, end = extract_ranges_from_name(peak_name)
                    out_file.write(f"{chromo}\t{start}\t{end}\n")


        case [_, "--help"] | _:
            print("Usage:")
            print("\tpython3 fantom5_utils.py filter <source_bed_path> <fantom_5_path>")
            print("\tpython3 fantom5_utils.py getRanges <peak_type> <source_types_path> <out_path> [TSS radius]")
            print("\tpython3 fantom5_utils.py countTypes <source_types_path>")
            print("\tpython3 fantom5_utils.py background <peak_type> <source_types_path> <out_path>")