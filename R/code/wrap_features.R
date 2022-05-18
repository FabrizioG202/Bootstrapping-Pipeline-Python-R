# Parsing the provided arguments.
#* This is done at the start of the file as the import of the required libraries below is really slow and
# we want to give the user the feedback immediately
library("optparse")

# For load_data
source("code/common.R")

# The arguments
option_list <- list(
    make_option(c("-d", "--data"),
        type = "character", default = NULL,
        help = "The input dataset [Required]", metavar = "character"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "The output file name. The path is relative to the directory of the input. [Required]", metavar = "character"
    ),

    # the cell line
    make_option(c("-c", "--line"),
        type = "character", default = NULL,
        help = "The cell line to use. [Required]", metavar = "character"
    ),

    # the feature type
    make_option(c("-f", "--feature"),
        type = "character", default = NULL,
        help = "The feature type to use. [Required]", metavar = "character"
    ),

    # The file extension (Passed through snakemake just because it is easier to use)
    make_option(c("-e", "--extension"),
        type = "character", default = NULL,
        help = "The file extension to use. [Required]", metavar = "character"
    )
)

# The argument parser
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



# Chacking that the input file has been passed
if (is.null(opt$data) | is.null(opt$output)) {
    stop("Please specify an input and an output file.")
}

# and that the cell line and feature type have been passed
if (is.null(opt$line) | is.null(opt$feature)) {
    stop("Please specify a cell line and a feature type!")
}

# The input file name
suppressMessages(library(crayon))
message("\tImporting libraries: ",  yellow("..."))


# If no errors are found in the input parameters, proceed with the script.
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
options(scipen = 999999999)
suppressMessages(library(valr))
suppressMessages(library(optparse))
message("\tImporting libraries: ",  green("Done"))

# the resolution datatable
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

#-----------------------------------------
## Utility function
wrap_features_fn <- function(input_file, extension) {

    # If the file is a bed file, then we can directly load it.
    if (extension == "bed") {
        bed_table <- read_bed(input_file)
        bed_table <- bed_table %>% filter(!(is.na(start)))
        features_grange <- with(bed_table, GRanges(chrom, IRanges(start, end)))
    } else if (extension == "rda") {
        # Load the data from the CAGE_file as an object called cage_coord_tbl
        features_table <- load_data(input_file)

        # Filtering the dataframe removing the rows with missing values in the field "start"
        features_table <- features_table %>% filter(!(is.na(start)))

        # Converting the data to genomic ranges
        features_grange <- GRanges(
            seqnames = features_table$chr,
            ranges = IRanges(
                start = as.numeric(features_table$start),
                end = as.numeric(features_table$end)
            )
        )
    } else {
        stop("The extension is not supported")
    }

    return(features_grange)
}

# The input file, read from the arguments
input_filepath <- "data/features/HMEC/CAGE/features.Rda"
input_filepath <- opt$data

# Parsing the output path from the provided argument
output_filepath <- opt$output

# the file extension
ext <- opt$extension
ext <- "rda"

# Running the function
wrapped_features <- wrap_features_fn(input_filepath, ext)

# Saving the file to disk.
save(wrapped_features, file = output_filepath)

# Print a success message and the output filepath
message("\tData saved to: ", blue(output_filepath))