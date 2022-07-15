# create a file containing the words abc, def.
library(optparse)

# the option list contains the following options:
# - input file
# - output file
# - annotations folder.
option_list <- c(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "", metavar = "character"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "", metavar = "character"
    ),
    make_option(c("-a", "--annotations"),
        type = "character", default = NULL,
        help = "", metavar = "character"
    )
)

# The optparser object
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# Check that none of the given arguments are empty
if (is.null(opt$input) | is.null(opt$output)) {
    message(opt$input)
    message(opt$output)
    stop("The input and output paths must be provided.")
}

suppressMessages(library(ChIPseeker))
library(crayon)
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))

# Calculates the number of features of a particular type.
# Requires the trascription database, the chromosome and the function file
annotate_features <- function(trascription_database, features, function_files) {

    # Annotates the features using the database
    #* suppress messages is there as the function throws a warning: 'select()' returned 1:many mapping between keys and columns # nolint
    peak_annotations <- suppressMessages(annotatePeak(features, tssRegion = c(-3000, 3000), TxDb = trascription_database, annoDb = "org.Hs.eg.db", verbose = FALSE))

    # Random annotations
    set.seed(seed)
    random_annotation <- sample(peak_annotations@annoStat$Feature, size = length(features), prob = peak_annotations@annoStat$Frequency / 100, replace = T)

    # check the number of peaks for each category
    n_5 <- length(grep("5'", as.character(random_annotation))) # 5' utr
    n_3 <- length(grep("3'", as.character(random_annotation))) # 3' utr
    n_exon <- length(grep("Exon", as.character(random_annotation))) # exon
    n_intron <- length(grep("Intron", as.character(random_annotation))) # intron
    n_1kb <- length(grep("1kb", as.character(random_annotation))) # 1kb
    n_2kb <- length(grep("2kb", as.character(random_annotation))) # 2kb
    n_3kb <- length(grep("3kb", as.character(random_annotation))) # 3kb
    n_down <- length(grep("Down", as.character(random_annotation))) # downstream
    n_inter <- length(grep("Inter", as.character(random_annotation))) # intergenic
    n_vec <- c(n_3, n_5, n_down, n_exon, n_inter, n_intron, n_1kb, n_2kb, n_3kb) # vector of counts
   
    # Adding the file names to the table
    names(n_vec) <- function_files

    # Clearing the local environment to prevent shadowing, the chipseeker environment may contain duplicates of the current # nolint
    # objects in the environment.
    rm(ChIPseekerEnv, envir = globalenv())
    rm(n_5, n_3, n_exon, n_intron, n_1kb, n_2kb, n_3kb, n_down, n_inter)

    # Keep only the features which are greater than 0
    #n_vec <- n_vec[n_vec > 0]
    return(n_vec)
}

# # (For the opt-parse version)
# feature_filepath <- "./data/features/HMEC/CAGE/features.bed"
# features <- valr::read_bed(feature_filepath)
# wrapped_features <- GRanges(seqnames = features$chrom, ranges = IRanges(features$start, features$end))
# annotations_folder <- "./data/annotations/"
# out_file <- opt$output

feature_filepath <- opt$input
features <- valr::read_bed(feature_filepath)
wrapped_features <- GRanges(seqnames = features$chrom, ranges = IRanges(features$start, features$end))
annotations_folder <- opt$annotations
out_file <- opt$output


seed <- 60

# All the counts of the vectors
all_counts <- c()

# The database containing the transcription data.
trascription_database <- TxDb.Hsapiens.UCSC.hg19.knownGene

# the available list of chromosomes
chromosome_list <- list.files(annotations_folder)

# Iterate over the chromosomes
chromo = "chr13"
for (chromo in chromosome_list) {

    # Filters the features retaining only the ones belonging to the given chromosome.
    features <- wrapped_features[seqnames(wrapped_features) == chromo]
    
    # The folder for the functional annotation
    annotation_folder <- paste0(annotations_folder, chromo, "/")

    # Get all the files in the directory containing the word BED in the title
    # The dollar as it needs to match at the end of the string
    function_files <- grep("BED$", list.files(annotation_folder), value = T)

    # Batch-parse the BED files.
    # It returns a list containing the parsed BED files
    function_data <- lapply(function_files, function(f) {
        valr::read_bed(paste0(annotation_folder, f), n_fields = 3)
    })

    # Updates the names.
    names(function_data) <- function_files

    # The transcription database
    chromosome_txdb <- trascription_database
    seqlevels(chromosome_txdb) <- chromo

    # Logging the start of functional annotation
    cat("Annotating chromosome: ", green(chromo), "\n", sep = "")

    # The number of annotations for each feature type
    count_vector <- annotate_features(chromosome_txdb, features, function_files)

    # Appending the counts to the list
    all_counts <- rbind(all_counts, count_vector)

    # resetting the transcript database
    seqlevels(trascription_database) <- seqlevels0(trascription_database)
}

# Creating the counts table.
# It stores counts in the given form
# chr  | 5'   |  3'  | Exon | Intron | 1kb  | 2kb  | 3kb  | Down | Inter | Total |
# chr1 | 1000 | 1000 | 1000 | 1000   | 1000 | 1000 | 1000 | 1000 | 1000  | 1000  |
counts_table <- data.frame(chr = chromosome_list, all_counts)

# Rename the column names to remove chr1 from them
# The j is added as column names cannot start with _ or numbers and j
# is a letter which does not appear in any other name in the columns so it should be safe to use.
# It's then replaced in the script compute p_value with the chromosome name and an _.
colnames(counts_table) <- c("chr", "j3utr.BED", "j5utr.BED", "jdown.BED", "jexon.BED", "jinter_no.BED", "jintron.BED", "jprom.BED", "jprom12.BED", "jprom23.BED")

# Save the data frame as a tsv file.
# The column names are preserved meaning that the headers are stored.
write.table(counts_table, out_file, sep = "\t", row.names = F, col.names = T, quote = F)

