# Empirical p-value process re-factored as functions
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(valr))
suppressMessages(library(crayon))
suppressMessages(library(future))
suppressMessages(library(future.apply))
library(optparse)

# For the load_data function
source("code/common.R")

# R-Specific stuff
options(scipen = 999999999)
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

# The options for the command line
option_list <- c(
    make_option(c("-d", "--clusters"),
        type = "character", default = NULL,
        help = "The folder containing the clusters", metavar = "character"
    ),
    make_option(
        c("-a", "--annotations"),
        type = "character", default = NULL,
        help = "The path to the annotations counts file", metavar = "character"
    ),
    make_option(
        c("-f", "--feature"),
        type = "character", default = NULL,
        help = "The path to the feature file", metavar = "character"
    ),
    make_option(
        c("-g", "--genome"),
        type = "character", default = NULL,
        help = "The path to the genome file", metavar = "character"
    ),
    make_option(
        c("-u", "--function_path"),
        type = "character", default = NULL,
        help = "The path to the function files", metavar = "character"
    ),
    make_option(
        c("-o", "--output"),
        type = "character", default = NULL,
        help = "The path to the output file", metavar = "character"
    ),
    make_option(
        c("-w", "--workers"),
        type = "integer", default = 3,
        help = "The number of workers to dispatch the work on", metavar = "integer"
    ),

    # The bootstrap count.
    make_option(
        c("-b", "--bootstraps"),
        type = "integer", default = 10,
        help = "The number of bootstrap iterations", metavar = "integer"
    ),

    # The randoms seed.
    make_option(
        c("-s", "--seed"),
        type = "integer", default = 1,
        help = "The random seed for the bootstrapping", metavar = "integer"
    ),

    # The chromosomes.
    make_option(
        c("-c", "--chromosomes"),
        type = "character", default = NULL,
        help = "The chromosomes to process", metavar = "character"
    )
)

# Creating the options to parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Parsing the arguments from the console
clusters_folder <- opt$clusters
annotations_counts_path <- opt$annotations
feature_file <- opt$feature
out_file <- opt$output
genome_filepath <- opt$genome
function_path <- opt$function_path
bootstraps <- opt$bootstraps
seed <- opt$seed
chromosomes <- opt$chromosomes

# the number of workers to dispatch the work on.
WORKERS_NUM <- opt$workers

# How many actual random genomes to generate:
random_genomes_num <- 100

# Generates the random functions.
rn_wrapped_features_build_fn <- function(chromo_counts, fn_bed_l, hg19_coord, tmp_cage_tbl, fn_file) {

    # The environment from the function
    fn_env <- environment()

    # The empty vector of functions
    rn_fn_coord_l <- vector("list", length(chromo_counts))
    rn_fn_coord_l
    names(rn_fn_coord_l) <- names(chromo_counts)

    for (f in names(chromo_counts)) {
        # How many feature do we need to generate?
        tmp_n <- chromo_counts[[f]]


        if (f == fn_file[5]) {
            future::plan(future::multisession, workers = WORKERS_NUM)

            # applies from 1 to 100 to
            rn_fn_coord_l[[f]] <- future_lapply(1:random_genomes_num, future.seed = seed, future.packages = c("dplyr", "valr"), future.envir = fn_env, function(x) {
                # Try pattern to not abort at the shuffling step
                ## Sampleing prior to shuffling ensures a greater likelihood of success
                # 1) Sample n features from the features table.
                # 2) tell the sampler how big is the chromosome (so that we don't sample outside it)
                # 3) tell the sampler the region in which to sample those features
                # 4) within as we want to shuffle within the chromosome
                rn_pol <- valr::bed_shuffle(tmp_cage_tbl %>% dplyr::sample_n(tmp_n), genome = hg19_coord, excl = fn_bed_l[[f]], within = T, max_tries = 1e6)
                return(rn_pol)
            })

            future::plan(future::sequential)
        }

        # If the feature is not chr22_inter_no
        if (f != fn_file[5]) {
            future::plan(future::multisession, workers = WORKERS_NUM)

            # Create a list of 100 random situations
            rn_fn_coord_l[[f]] <- future_lapply(1:random_genomes_num, future.seed = seed, future.packages = c("dplyr", "valr"), future.envir = fn_env, function(x) {
                # Try pattern to not abort at the shuffling step
                ## Sampleing prior to shuffling ensures a greater likelihood of success
                # 1) Sample n features from the features table.
                # 2) tell the sampler how big is the chromosome (so that we don't sample outside it)
                # 3) tell the sampler the region in which to sample those features
                # 4) within as we want to shuffle within the chromosome
                # out a tibble of size n
                # (this basically moves the intervals from tmp_cage_tbl in places which are contained in the regions defined by fn_bed_l[[f]])
                rn_pol <- try(valr::bed_shuffle(tmp_cage_tbl %>% sample_n(tmp_n), genome = hg19_coord, incl = fn_bed_l[[f]], within = T, max_tries = 1e6), silent = TRUE)
                return(rn_pol)
            })

            future::plan(future::sequential)

            # Collect the successful shuffling by eliminating the shuffles producing try-error objects
            rn_fn_coord_l[[f]] <- rn_fn_coord_l[[f]][!(unlist(lapply(rn_fn_coord_l[[f]], function(x) any(class(x) %in% "try-error"))))]
        }

        cat("\t\t-Shuffled ", crayon::blue(f), " features", "\n")
    }

    warnings()

    # Assemble these blocks into n random results
    bootstrap_count <- bootstraps

    # Create 10 random results by:
    # looping over 10 and joining the results from the random sampling of one of the features.
    # For example, loop over chr22_inter_no, chr22_inter_no_1, chr22_inter_no_2, ..., chr22_inter_no_10
    # And join the results obtained from sampling one of the 100 available chr22_inter_no features
    future::plan(future::multisession, workers = WORKERS_NUM)
    rn_peak_coord_tbl_l <- future_lapply(1:bootstrap_count, future.seed = seed, future.envir = fn_env, future.packages = c("dplyr"), function(x) {
        return(do.call(dplyr::bind_rows, lapply(rn_fn_coord_l, function(f) f[[sample(1:length(f), 1)]])))
    })
    future::plan(future::sequential)


    # Remove the list of random functional features from this environment
    rm(rn_fn_coord_l, envir = fn_env)

    # convert the list of random functional features into GRanges
    rn_Grange_l <- lapply(rn_peak_coord_tbl_l, function(x) {
        rnp_Grange <- GRanges(
            seqnames = x$chrom,
            ranges = IRanges(
                start = x$start,
                end = x$end
            )
        )
        return(rnp_Grange)
    })

    # Remove the original list
    rm(rn_peak_coord_tbl_l, envir = fn_env)

    return(rn_Grange_l)
}

## The main function
# chromo is the chromosome name, so for example chr1
# clusters_folder is the folder where the clusters are contained
# clusters_file_suffix is the suffix of the cluster file
# wrapped_features is the GRange object of the features
# function_path is "../data/annotation/", the folder where the annotation files are
empirical_pval_compute_fn <- function(chromo, clusters_folder, wrapped_features, function_path, hg19_coord, annotations_counts) {

    # The main environment of the function
    main_fn_env <- environment()

    # the features of this chromosome
    chr_wrapped_features <- wrapped_features[seqnames(wrapped_features) == chromo]

    # Logging the start of the bootstrappinf
    cat("\tBootstrapping for chromosome ", crayon::green(chromo), "\n")

    # The folder containing the annotations about this chromosome: "../data/annotation/fn_BED/chr1/"
    fn_folder <- paste0(function_path, chromo, "/")

    # All the bed files for this chromosome
    fn_file <- grep("BED$", list.files(fn_folder), value = T)

    # Parsing all the bed files
    # The bed files are parsed and put in a dictionary where the key is the file name and the value is the bed file
    fn_bed_l <- lapply(fn_file, function(f) {

        # Just read the chromosome name, the start and the end
        valr::read_bed(paste0(fn_folder, f), n_fields = 3)
    })
    names(fn_bed_l) <- fn_file

    # Select the row where chromosome is chromo
    chromo_counts <- annotation_counts[annotation_counts$chr == chromo, ]

    # Keep only the element starting from the second one (remove chr)
    chromo_counts <- chromo_counts[, 2:ncol(chromo_counts)]

    # Replace X with chromo in the names of chromo_counts (for some reason while importing the column names)
    # switch from utr3.BED to chr1_utr3.BED
    names(chromo_counts) <- gsub("j", paste0(chromo, "_"), names(chromo_counts))

    # Loading the clusters for this chromosome
    # find the first file in the cluster directory which contains the chromosome name
    clusters_file <- NA
    for (f in list.files(clusters_folder)) {
        if (grepl(paste0(chromo, "_"), f)) {
            clusters_file <- f
            break
        }
    }
    if (is.na(clusters_file)) {
        cat("\tNo clusters file found for chromosome ", red(chromo), "in ", yellow(clusters_folder), ", Skipping it!\n")
        return(NA)
    }

    #cl_chr_tbl <- load_data(paste0(clusters_folder, "/", clusters_file))
    cl_chr_tbl = load_data("./data/tmp/wrapped_clusters/chr22_wrapped.Rda")

    # Logging that the counting of the observed overlaps has started.
    cat("\tCunting observed overlaps for chromosome ", crayon::green(chromo), "\n")

    # table collecting the observed CAGE-peak coordinates
    # counting the number of overlaps between the features and the clusters
    # looping over each list of GRanges for each cluster
    future::plan(future::multisession, workers = WORKERS_NUM)
    
    cl_inter_vec <- unlist(future_lapply(cl_chr_tbl$GRange, future.envir = main_fn_env, future.packages = c("GenomicRanges"), function(x) {

        # Counting the overlaps between a given feature and a cluster
        sum(IRanges::countOverlaps(GRangesList(x), chr_wrapped_features))
    }))
    future::plan(future::sequential)

    cl_chr_rbl_list = GRangesList(cl_chr_tbl$GRange)
    overlaps = countOverlaps(cl_chr_rbl_list, chr_wrapped_features)
    
    cl_inter_vec = unlist(lapply(cl_chr_tbl$GRange, function(x){
        sum(IRanges::countOverlaps(GRangesList(x), chr_wrapped_features))
    }))

    sum(cl_inter_vec)

    # Adds a column containing the number of overlaps to the cluster table and filters for the clusters having more than zero overlap
    cl_chr_tbl <- cl_chr_tbl %>%
        mutate(feature_n = cl_inter_vec) %>%
        filter(feature_n > 0)
    

    # Converting the Granges to lists.
    cl_list <- GRangesList(cl_chr_tbl$GRange)

    # Create this tmp_cage_tbl by:
    # -1) Converting the chr_wrapped_features to a tibbble
    # -2) Selecting only the seqname, start and end columns
    # -3) Renaming the column seqname to chrom
    tmp_cage_tbl <- chr_wrapped_features %>%
        as_tibble() %>%
        dplyr::select(seqnames, start, end) %>%
        dplyr::rename(chrom = seqnames)

    # Logging that the generation of the random features
    cat("\tGenerating random features for chromosome ", crayon::green(chromo), "\n")

    # Building the effective random coordinates on which the p-value will be computed
    rn_Grange_l <- rn_wrapped_features_build_fn(chromo_counts, fn_bed_l, hg19_coord, tmp_cage_tbl, fn_file)

    # loop over all the files in random_genomes
    beds = c()
    overlaps = matrix(rep(NA, nrow(cl_chr_tbl)* 10000), nrow = nrow(cl_chr_tbl), ncol = 10000)
    
    
    for (f in list.files("./random_genomes/")){
        bed_data = valr::read_bed(paste0("./random_genomes/", f))
        for (i in 1:100){
            i = 1
            bed_data_i = bed_data[((i-1)*1926):(i*1926),]
            ranges_ = GRanges(seqnames = bed_data_i$chrom, ranges = IRanges(start = bed_data_i$start, end = bed_data_i$end))
            overlaps = countOverlaps(cl_list, ranges_)
        }        
    }

    # Logging the started of the p-value computation
    cat("\tComputing p-values for chromosome ", crayon::green(chromo), "\n")

    # counting the overlaps between the clusters and the newly created random features
    # The resulting vector will be a list of 10 (number of bootstrap) lists of overlaps count for each cluster.
    future::plan(future::multisession, workers = WORKERS_NUM)
    rn_pval_l <- future_lapply(rn_Grange_l, future.packages = c("GenomicRanges"), function(x) {
        IRanges::countOverlaps(cl_list, x)
    })
    future::plan(future::sequential)

    rn_pval_l <- lapply(beds, function(x) {
        IRanges::countOverlaps(cl_list, x)
    })

    # Removing the useless granges.
    rm(rn_Grange_l)

    # Stack the pvalues one after the other one column is one tow of the rnpval_l
    rn_count_mat <- do.call(cbind, rn_pval_l)

    # Remove the useless list
    rm(rn_pval_l)

    # unlist converts a list to a vector.
    # Apply over the list of clusters, the function which calculates the p-value
    #
    cl_emp_pval <- unlist(lapply(1:nrow(cl_chr_tbl), function(x) {

        # x is the index
        # sums the number of times the overlaps are greater in the randomly generated one than in the observed ones and dividing
        # everything by the number of bootstrap samples.
        # So if the random have a value which is greater than the observed one in 5/10 cases, this will return 0.5
        (sum(rn_count_mat[x, ] >= cl_chr_tbl$feature_n[x]) + 1) / (ncol(rn_count_mat) + 1)
    }))

    # Add the pvalue as a column named emp.pval in the cluster table.
    return(cl_chr_tbl %>% mutate(emp.pval = cl_emp_pval))
}

# The features in the shape of a GRange object
wrapped_features <- load_data(feature_file)

# the coordinates of the h19 genome (chromosome names and sizes)
hg19_coord <- read_delim(genome_filepath,
    "\t",
    escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE,
    #show_col_types = FALSE # to quiet the mesage (crashes on the cluster with R 4.0.2)
)
names(hg19_coord) <- c("chrom", "size")

# If we provided a list of chromosomes to liimit the choice to, we do it
if (!is.null(chromosomes)) {

    #The chromosome list is expected in the form chr1, chr2, chr3, ...
    chr_set <- strsplit(chromosomes, ",")[[1]]
} else {

    # Otherwise just list the folders in the function directory
    chr_set <- list.files(function_path)
}

# Log the chromosomes for which we will compute the p-values
cat("\tChromosomes to compute p-values on: ", chr_set, "\n")

# Getting the annotation counts form the annotation file
annotation_counts <- read.csv(annotations_counts_path, header = T, stringsAsFactors = F, sep = "\t")

# Computing the p-values for all the chromosomes.
cl_chr_emp_pval_l <- lapply(chr_set, function(chromo) {
    chromo = chr22
    cl_chr_tbl <- empirical_pval_compute_fn(chromo, clusters_folder, wrapped_features, function_path, hg19_coord, annotation_counts)
    return(cl_chr_tbl)
})

# Remove the elements which are NAs (NAs can be produced at just one step right now which is the user is trying to comput the p-value on a chromosome
# for which there is no data available).
cl_chr_emp_pval_l <- cl_chr_emp_pval_l[!is.na(cl_chr_emp_pval_l)]

# Print the number of chromosomes which have been removed.
cat("\tNumber of chromosomes removed: ", nrow(cl_chr_emp_pval_l) - nrow(chr_set), "\n")

# "Stack" the dataframe one on top of each other.
cl_chr_emp_pval_tbl <- do.call(bind_rows, cl_chr_emp_pval_l)

# Save the table as a tsv
# Drop the GRange and bins columns
cl_chr_emp_pval_tbl <- cl_chr_emp_pval_tbl %>%
    dplyr::select(chr, cl, feature_n, emp.pval)

# Save the Dataframe to a table.
write.table(cl_chr_emp_pval_tbl, file = out_file, sep = "\t", row.names = F, col.names = T)