### Checks the likely N at which the file was run at.
import pandas as pd
import numpy as np
import sys
from src.crayon import *

def eco_format(n : int) -> str:
    """
    Formats the number like it's used in economics.
    """
    suffixes = ['', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    for i in range(len(suffixes)):
        if n < 1000:
            return f'{n}{suffixes[i]}'
        n /= 1000
    return f'{n:.2f}Y'

if __name__ == "__main__":
    argv = sys.argv

    match argv:
        case [_, path]:
            df = pd.read_csv(path, sep="\t")
            p_values = df["p_value"]
            min_p_value = np.min(p_values)
            print("Minimum pvalue: {}".format(min_p_value))
            print("Deduced N: {}".format((blue(eco_format(int(1/min_p_value))))))
        case _:
            print("Usage: check_n.py <path>")