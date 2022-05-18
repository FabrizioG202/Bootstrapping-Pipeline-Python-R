# %%
from pyjaspar import jaspardb
import requests
import xml.etree.ElementTree as ET

# %%
def extract_sequence_from_position(chromosome : str, start : int, end : int) -> str:
    
    """
    Returns the sequence from the given position in the given chromosome using UCSC.
    """
    
    url = "https://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment={}:{},{}".format(chromosome, start, end)
    response = requests.get(url).text
    
    # Decode the response as XML
    
    root = ET.fromstring(response)
    
    # Print all the elements in the XML
    for element in root.iter():
        if element.tag == "DNA":
            return element.text.strip()


# %%
import os
def _cache_file_from_position(chromosome : str, start : int, end : int) -> str:    
    return "cache/{}_{}_{}.fa".format(chromosome, start, end)

# %%
class SequenceUtils:

    @staticmethod
    def get_coordinates(chromosome : str, start : int, end : int) -> str:
        """ Checks if the given coordinates are cached and returns the sequence if they are, otherwise it downloads the sequence from UCSC and caches it. """

        # Check if the coordinates are cached
        cached_file = _cache_file_from_position(chromosome, start, end)

        #Check that the cache directory exists
        if not os.path.exists("cache"):
            os.makedirs("cache")

        # Try to retrieve it.
        if os.path.isfile(cached_file):
            # If they are, return the sequence
            with open(cached_file, "r") as f:
                return f.read()
        else:
            # If they are not, download the sequence and cache it
            sequence = extract_sequence_from_position(chromosome, start, end)
            with open(cached_file, "w") as f:
                f.write(sequence)
            return sequence

    @staticmethod
    def query_human_motifs():
        jaspar_database = jaspardb(release="JASPAR2022")
        motifs = jaspar_database.fetch_motifs(species=["9606"])
        return motifs