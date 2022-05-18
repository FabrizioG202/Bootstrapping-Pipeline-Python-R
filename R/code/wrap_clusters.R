# Parsing the provided arguments.
#* This is done at the start of the file as the import of the required libraries below is really slow and
# we want to give the user the feedback immediately
library("optparse")

# The arguments
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    help = "The input directory"
  ),
  make_option(c("-o", "--output"),
    type = "character",
    help = "The output directory"
  ),
  # The cell line
  make_option(c("-l", "--line"),
    type = "character",
    help = "The cell line to use"
  ),
  make_option(c("-p", "--progress"), default =  FALSE,
    help = "Show progress (Showing it will block parallelism)"
  )
)

# The argument parser
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that none of the passed parameters are empty
if (is.null(opt$input) || is.null(opt$output)) {
  stop("Please provide the input and output directories")
}

# Make the directory if it does not exist
if (!file.exists(opt$input)) {
  dir.create(opt$input, showWarnings = F)
}


# Getting the parameters
input_directory <- "data/clusters/HMEC/"
output_directory <- "data/cluster/wrapped"
input_directory <- opt$input
output_directory <- opt$output
progress <- opt$progress

# The input filename
suppressMessages(library(crayon))
message("\tImporting libraries: ",  yellow("..."))

# Produce Coordinate tables for BHiCect clusters
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(furrr))
suppressMessages(library(crayon))
message("\tImporting libraries: ",  green("Done"))

has_progress <- TRUE
# Try to import the library progress
tryCatch(library(progress), error = function(e) {
  has_progress <- FALSE
})

# Show the progress only if it is requested
has_progress = F
has_progress = has_progress && progress

library(stringr)
options(scipen = 999999999)

# The resolution datatable
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set


# Planning the execution on 3 worker sessions
if (!has_progress) plan(multisession, workers = 3)

file <- list.files(input_directory)[[1]]

# Looping over the chromosomes
for (file in list.files(input_directory)) {

  # The chromosome name
  chr <- strsplit(file, ".", fixed = TRUE)[[1]][1]

  # grab the name of the chromosome from the file name.
  # the filename is like this: chr10_spec_res.txt
  # the chromosome name is chr10
  # the RegExp grabs everything from the start of the file to the first _
  chromo <- strsplit(chr, "_", fixed = TRUE)[[1]][1]

  filepath <- paste0(input_directory, "/", file)
  filepath = "./data/clusters/HMEC/chr21_spec_res.Rda"

  # The output file
  output_file <- paste0(output_directory, "/", chr, "_wrapped.Rda")

  # loading the chromosome file (loads an object called chr_spec_res)
  load(filepath)

  # Converting the chromosome file to a tibble
  chr_cl_tbl <- tibble(chr = chromo, cl = names(chr_spec_res$cl_member), bins = chr_spec_res$cl_member)
 
  # Add a table representing the resolution at which the cluster was detected to the table
  chr_cl_tbl <- chr_cl_tbl %>% mutate(res = map_chr(cl, function(x) {
    return(strsplit(x, split = "_")[[1]][1])
  }))

  # Adding a GRange object to the table
  cat(paste0("\tChromosome: ", blue(chromo), "\n"))

  # If the user has the library progress, print the progress bar, otherwise use the future package and do not print the progress.
  if (has_progress) {
    pb <- progress_bar$new(format = "  wrapping [:bar] :percent eta: :eta", clear = FALSE, total = nrow(chr_cl_tbl) + 1)
    pb$tick(0)

    chr_cl_tbl <- chr_cl_tbl %>% mutate(GRange = pmap(list(chr, bins, res), function(chr, bins, res) {
     # pb$tick()

      ranges = (GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = as.numeric(bins),
          end = as.numeric(bins) + res_num[res] - 1
        )
      ))
    }))
  } else {
    #cat("\t(Running in Parallel)", "\n")
    chr_cl_tbl <- chr_cl_tbl %>% mutate(GRange = future_pmap(list(chr, bins, res), function(chr, bins, res) {
      return(GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = as.numeric(bins),
          end = as.numeric(bins) + res_num[res] - 1
        )
      ))
    }))
  }

  # Saving the file to disk
  save(chr_cl_tbl, file = "./data/tmp/wrapped_clusters/chr21_wrapped.Rda")
}

if(!has_progress) plan(sequential)