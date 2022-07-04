from src.crayon import *
import sys
from src.fantom5_io import parse_fantom_peak_type_data, parse_fantom_data

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

def count_peak_types(path : str) -> tuple[int, int]:
    """
    Count the number of peaks in a peak file.
    """
    tss_ = 0
    enhancer_ = 0
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            else:
                type_ = line.split("\t")[1].strip()
                if type_ == "tss":
                    tss_ += 1
                elif type_ == "enhancer":
                    enhancer_ += 1
                else:
                    raise ValueError("Unknown peak type: {}".format(type_))
    return tss_, enhancer_

if __name__ == "__main__":
    argv = sys.argv

    match argv:
        
        # Gets the ranges in the provided [source_bed_path] which are of the asked type (tss or enhancer)
        # And writes the filtered one to output path
        case [_, "getRanges", "enhancer" | "tss"  as peak_type, source_bed_path, fantom_5_peak_types_path, out_path, *args]:

            #Parse the fantom peak type data
            peak_types = parse_fantom_peak_type_data(fantom_5_peak_types_path)

            # ask for a radius if the user required the TSS.
            tss_radius = None
            if peak_type == "tss":
                if len(args) == 0:
                    print(red("Please provide a radius for the TSS extraction!"))
                    exit(1)
                else:
                    tss_radius = int(args[0])
            
            # Diagnostics:
            not_in_fantom_5 = 0
            count = 0
            with open(source_bed_path, "r") as source_bed, open(out_path, "w") as out_file:
                for i, line in enumerate(source_bed.readlines()):
                    
                    # Skip header
                    if line.startswith("#"):
                        continue
                    
                    # read the bed line
                    chromo, start, end, id_, *rest = line.split("\t")
                    chromo = chromo.strip()
                    start = int(start)
                    end = int(end)

                    # Extract the id of the peak
                    id_ = id_.strip()

                    # Get the peak type
                    this_peak_type = peak_types.get(id_, None)

                    # If the peak type is not in the fantom 5 peak types, skip it
                    if this_peak_type is None:
                        not_in_fantom_5 += 1
                        continue

                    # now check that the peak is of the requested type.
                    if this_peak_type != peak_type:
                        continue

                    if peak_type == "tss":
                        start = start - tss_radius
                        end = end + tss_radius
                    
                    # Write the line to the output file
                    out_file.write("{}\t{}\t{}\t{}\n".format(chromo, start, end, id_))
                    count += 1

            
            print(f"Successfully extracted {green(count)} ranges to: " +  out_path  + " with {} peaks not mapped in fantom 5".format(red(not_in_fantom_5)))

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

        # Saves a list of the entrez ids of the peaks in the provided [source_bed_path]
        case [_, "getEntrez", source_bed_path, peak_to_entrez_path, out_path]:
        
            

            # the resulting ids
            ids = []

            # Reading the peaks to entrez ids
            peaks_to_entrez = parse_fantom_data(peak_to_entrez_path)

            # diagnostics
            mapped = 0
            unmapped = 0
            nas = 0
            
            # loop over the lines of the file (they are in a type-like format (CAGE.ID, peak_type)))
            with open(source_bed_path, "r") as source_bed:
                    
                for line in source_bed.readlines():

                    # skip headers
                    if line.startswith("#"):
                        continue
                        
                    # split the line
                    chromo, start, end, peak_id, *rest = line.split("\t")
                    peak_id = peak_id.strip()
                    
                    # get the entrez id for the peak name
                    entrez_id = peaks_to_entrez.get(peak_id, None)
                    
                    if entrez_id is None:
                        unmapped += 1
                        continue

                    elif entrez_id == "NA":
                        nas += 1
                        continue

                    # add the entrez id to the list
                    ids.append(entrez_id)
                    mapped += 1

            # Remove the duplicates
            ids = list(set(ids))

            # write the ids to a file
            with open(out_path, "w") as out_file:
                for id in ids:
                    out_file.write(f"{id}\n")


            # Print some diagnostics
            print(f"Mapping Finished: Mapped: {green(mapped)}, Unmapped: {red(unmapped)}, NA: {yellow(nas)}, file saved to {out_path}")
       
        case [_, "fullBackground", "enhancer" | "tss" as peak_type, fantom_5_peak_types_path, out_path, ]:
            
            with open(fantom_5_peak_types_path, "r") as fantom_5_peaks, open(out_path, "w") as out_file:
                for i, line in enumerate(fantom_5_peaks.readlines()):
                    if line.startswith("#"):
                        continue

                    peak_name,  pt = line.split("\t")
                    peak_name = peak_name.strip()
                    pt = pt.strip()

                    # Get only the required peak types
                    if peak_type != pt:
                        continue
                    
                    # Get the chromo, start and end
                    chromo, start, end = extract_ranges_from_name(peak_name)

                    # Write the line to the file
                    out_file.write(f"{chromo}\t{start}\t{end}\t{peak_name}\n")


        case [_, "typesTable", foreground_types_path, background_types_path, out_path]:
            foreground_types = count_peak_types(foreground_types_path)
            background_types = count_peak_types(background_types_path)

            with open(out_path, "w") as out_file:
                # Add the header
                out_file.write("\t".join(["type", "foreground", "background"]) + "\n")
                out_file.write("\t".join(["enhancer", str((foreground_types[1])), str((background_types[1]))]) + "\n")
                out_file.write("\t".join(["tss", str((foreground_types[0])), str((background_types[0]))]) + "\n")

        case [_, "--help"] | _:
            print("Usage:")
            print("\tpython3 fantom5_utils.py filter <source_bed_path> <fantom_5_path>")
            print("\tpython3 fantom5_utils.py getRanges <peak_type> <source_types_path> <out_path> [TSS radius]")
            print("\tpython3 fantom5_utils.py countTypes <source_types_path>")
            print("\tpython3 fantom5_utils.py background <peak_type> <source_types_path> <out_path> [TSS radius]")