from struct import pack, unpack

def count_file_rows(file_name : str) -> int:
    """
    Counts the number of rows in a file without reading the whole file, useful for large files.
    """
    with open(file_name, 'r') as f:
        return sum(1 for line in f)

def read_genome_file(file_name : str) -> dict[str, int]:
    """
    Reads a genome file and returns a dictionary with the chromosome names as keys and the chromosome lengths as values.
    """
    with open(file_name, 'r') as f:
        return {line.split("\t")[0]: int(line.split("\t")[1]) for line in f}

def read_counts_file(path: str, chromosome : str) -> dict[str, int]:
    """
    Reads a file with counts and returns a dictionary with the chromosome names as keys and the chromosome lengths as values.
    """
    with open(path, 'r') as f:
        columns = f.readline().strip().split("\t")[1:]
        columns = [c.replace("j", "").replace(".BED", "") for c in columns]

        for line in f.readlines():
            _chr = line.split("\t")[0]
            parts = line.split("\t")[1:]
            if _chr == chromosome:
                return {columns[i]: int(parts[i]) for i in range(len(parts))}

        raise ValueError("Chromosome not found in file.")

def resolution_string_to_int(resolution : str) -> int:
    """
    Converts a resolution string to an integer value:
    eg: 10kb -> 10000
    """
    resolution = resolution.lower()
    if resolution.endswith("kb"):
        return int(resolution[:-2]) * 1000
    elif resolution.endswith("mb"):
        return int(resolution[:-2]) * 1000000
    elif resolution.endswith("gb"):
        return int(resolution[:-2]) * 1000000000
    else:
        raise ValueError("Unknown resolution: " + resolution)

import numpy as np

def read_bed_as_array(path : str, chromo : str):
    # create the array
    data = np.zeros((0, 3), dtype=np.int64, order='F')

    # open the file
    with open(path) as f:
        # skip the header
        f.readline()
        
        # read the data
        for i, line in enumerate(f):
            _chromo, start, end, *_ = line.split("\t")

            # Chromosome data.
            if _chromo != chromo:
                continue
            
            # Add the data.
            data = np.append(data, [[int(start), int(end), i]], axis=0)
            
    # return the data
    return data

def format_bytes(num : int) -> str:
    """
    Formats a number of bytes to a human readable string.
    """
    for x in ['B', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0
    return "%3.1f %s" % (num, 'PB')


def save_matrix(array : np.ndarray, path : str):
    with open(path, 'wb') as file:
        shape = array.shape
        #write the shape as the first bytes
        file.write(pack('i', shape[0]))
        file.write(pack('i', shape[1]))
        file.write(pack('d' * len(array.flatten()) , *(array.flatten().tolist())))
        
def read_matrix(path : str) -> np.ndarray:
    with open(path, 'rb') as file:
        shape_x = unpack('i', file.read(4))[0]
        shape_y = unpack('i', file.read(4))[0]
        
        packed = file.read()[8:]
        array = unpack('d' * (len(packed) // 8), packed) # 8 bytes per double
        return np.array(array).reshape(shape_x, shape_y)