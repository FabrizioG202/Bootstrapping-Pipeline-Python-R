# The table for converting resolutions to string.
RESOLUTIONS = {
    "5kb": 5e3,
    "10kb": 1e4,
    "50kb": 5e4,
    "100kb": 1e5,
    "500kb": 5e5,
    "1Mb": 1e6,
}

# Extracts bins from the strange json storage method.
def extract_bins( x: str | list) -> list[int]:
    if type(x) == str:
        return [int(x)]
    elif type(x) == list:
        return [int(i) for i in x]
    else:
        raise TypeError(f"The type of the argument is not a string or a list. It is: {type(x)}")

def __chr_sort_part(x : str):
    p = x.replace("chr", "")
    if p.isdigit():
        return int(p)
    else:
        return p

def sort_by_chromosome(x : dict | list):
    if type(x) == dict:
        # Return a sorted dict.
        return {k: v for k, v in sorted(x.items(), key=lambda item: __chr_sort_part(item[0]))}
    elif type(x) == list:
        return sorted(x, key=lambda x: __chr_sort_part(x))
    else:
        raise TypeError(f"The type of the argument is not a dict or a list. It is: {type(x)}")