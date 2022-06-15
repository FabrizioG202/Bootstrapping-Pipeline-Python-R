import liftover
import sys
import os
from src.crayon import *

def read_bed(path : str) -> list[tuple]:
    """
    Reads a bed file and returns a list of tuples.
    """
    bed_list = []
    with open(path, 'r') as bed:
        for line in bed:
            # Skip header, and empty lines and lines containing UCSC data
            if line.startswith('#') or line.startswith('browser') or line.startswith('track') or line.startswith('ucsc'):
                continue
            line = line.strip().split('\t')
            bed_list.append((line[0], int(line[1]), int(line[2]), *line[3:]))
    return bed_list

def lift_bed_list(bed_list : list[tuple], lifter) -> list[tuple]:
    """
    Takes a list of tuples and lifts them using the lifter object.
    """
    lifted_list = []
    failed_entries = 0
    one_to_many = 0
    for bed in bed_list:
        chrom, start, end, *rest = bed
        new_start = lifter[chrom][start]
        new_end = lifter[chrom][end]
        if len(new_start) == 0 or len(new_end) == 0:
            failed_entries += 1
            continue
        elif len(new_start) > 1 or len(new_end) > 1:
            if len(new_start) != len(new_end):
                continue
            else:
                one_to_many += 1
        
      
        for new_start_pos, new_end_pos in zip(new_start, new_end):
            if new_start_pos > new_end_pos:
                failed_entries += 1
                continue

            lifted_list.append((chrom, new_start_pos[1], new_end_pos[1], *rest))
    
    return lifted_list, failed_entries, one_to_many

def save_bed(bed_list : list[tuple], path : str):
    """
    Saves the bed list to a file.
    """
    with open(path, 'w') as bed_file:
        for bed in bed_list:
            bed = '\t'.join([str(s) for s in bed])
            bed += '\n'
            bed_file.write(bed)

if __name__  == "__main__":
    # Get the arguments
    argv = sys.argv

    match argv:

        case [_, file_path, from_assembly, to_assembly, *args]:
            replace = "-r" in args
            bed_data = read_bed(file_path)
            
            if replace:
                out_path = file_path
            else:
                out_path = args[0] if len(args) > 0 else os.path.splitext(file_path)[0] + '_' + to_assembly + '.bed'

            # Lift over the data
            try:
                lifter = liftover.get_lifter(from_assembly, to_assembly)
            except:
                print(red(f"cannot lift over {from_assembly} to {to_assembly} (are you sure those exists?)"))
                sys.exit(1)

            lifted_data, failed, one_to_many = lift_bed_list(bed_data, lifter)

            # Save the lifted data
            save_bed(lifted_data, out_path)
            
            # Log the results: "Successfully lifted j (green) entries with x (red) failures and y (yellow) one to many"
            print(f"Successfully lifted {green(len(lifted_data))} entries with {red(failed)} failures and {yellow(one_to_many)} one to many")


        case [_, "help"| "h"| "--help"] | _:
            print(
                "Usage: python lift_over.py <bed_file> <from_assembly> <to_assembly>\n"
            )