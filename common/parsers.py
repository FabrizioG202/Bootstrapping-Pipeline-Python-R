#%%
#Parses the exclusion or filtering files
class ExclusionPattern:
    excluded_lines : str
    included_lines : str

    excluded_types : str
    included_types : str

    excluded_target : str
    included_target : str

    def __init__(self, excluded_lines, included_lines, excluded_types, included_types, excluded_target, included_target):
        self.excluded_lines = excluded_lines
        self.included_lines = included_lines

        self.excluded_types = excluded_types
        self.included_types = included_types

        self.excluded_target = excluded_target
        self.included_target = included_target

    def __str__(self):
        return f"Excluded lines: {self.excluded_lines}\nIncluded lines: {self.included_lines}\nExcluded types: {self.excluded_types}\nIncluded types: {self.included_types}\nExcluded target: {self.excluded_target}\nIncluded target: {self.included_target}"

    def filter(self, requests: list[tuple[str, str, str]]):
        # The passed requests are a list of tuples of the form (line, type, filename)
        # filter out the excluded lines and add the included lines

        results = []
        
        for request in requests:
            # filter out based on the active target (highest priority):
            target_like = f"{request[0]}/{request[1]}"
            if target_like in self.excluded_target:
                continue
            elif target_like in self.included_target:
                results.append(request)
                continue

            # filter out based on the active types (second priority):
            type_like = f"{request[1]}"
            if type_like in self.excluded_types:
                continue
            elif type_like in self.included_types:
                results.append(request)
                continue
            
            # filter out based on the active lines (lowest priority):
            line_like = f"{request[0]}"
            if line_like in self.excluded_lines:
                continue
            else:
                results.append(request)
        return results

# Parses the exclusion file, this is an example of it.
# #exclude:
# line: HMEC
# HMEC/CAGE
# type: CAGE
def parse_filter_file(path : str):
    with open(path, "r") as f:
        lines = f.readlines()
        excluded_lines = []
        included_lines = []
        excluded_types = []
        included_types = []
        excluded_target = []
        included_target = []

        mode = "exclude"

        for line in lines:
            # skip empty lines
            if line.strip() == "":
                continue

            # read the mode
            if line.startswith("#"):
                mode = line.strip()[1:].lower().replace(":", "")
                continue
            
            match line.split(" "):
                case ["line:", a]:
                    if mode == "exclude":
                        excluded_lines.append(a.strip())
                    else:
                        included_lines.append(a.strip())
                case ["type:", a]:
                    if mode == "exclude":
                        excluded_types.append(a.strip())
                    else:
                        included_types.append(a.strip())
                case other:
                    if mode == "exclude":
                        excluded_target.append(other[0].strip())
                    else:
                        included_target.append(other[0].strip())
    return ExclusionPattern(excluded_lines, included_lines, excluded_types, included_types, excluded_target, included_target)

PATH = "./filter.txt"

print(parse_filter_file(PATH))

requests = [
    ("HMEC", "CAGE", "hmec_cage.txt"),
    ("H1", "CAGE", "h1_cage.txt"),
    ("HMEC", "CTCF", "hmec_ctcf.txt"),
    ("H1", "CTCF", "h1_ctcf.txt"),
    ("HMEC", "DNase", "hmec_dnase.txt"),
    ("H1", "DNase", "h1_dnase.txt"),
]

print("Filtering: ")
print("\n".join(map(str, parse_filter_file(PATH).filter(requests))))