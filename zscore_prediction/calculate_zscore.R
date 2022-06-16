library(vroom)
library(mgcv)
library(tidyverse)
options(scipen = 999999999)
library(crayon)

# Read and filter sparse HI-C Data.
load_sparse_data <- function(dat_file) {

    # Just read the data
    # This produces a 3 Columned file:
    # X1 = Bin 1 start
    # X2 = Bin 2 start
    # X3 = Count (Score)
    chr_dat <- vroom(dat_file, delim = "\t", col_names = F, trim_ws = T, escape_double = F)

    # ensure HiC score variable is a numeric, not NA and filter self interactions.
    return(chr_dat %>%
        mutate(X3 = as.numeric(X3)) %>% filter(!(is.na(X3))) %>% filter(X1 != X2))
}

compute_zscore <- function(dat_file, ncluster) {

    # Load desired HiC data
    chr_dat <- load_sparse_data(dat_file) %>%
        # Compute genomic ditance between interacting bins and add it as a new column (d)
        mutate(d = abs(X1 - X2)) %>%
        # Perform log transform on distance and HiC score (better behaved input for GAM-computation)
        mutate(lw = log10(X3), ld = log10(d))

    # Compute GAM
    hic_gam <- bam(lw ~ s(ld, bs = "ad"), data = chr_dat,)
    pred_vec <- predict(hic_gam, newdata = chr_dat)
    
    # Compute zscore and predicted HiC magnitude
    chr_dat <- chr_dat %>%
    
        # add predicted HiC values from GAM
        mutate(pred = pred_vec) %>%
    
        # Compute zscore using estimation of variance (hic_gam$sig2) and predicted HiC score (pred) from GAM
        mutate(zscore = (lw - pred) / hic_gam$sig2)
    
    # drop the X3, d, lw, ld, and pred columns
    chr_dat <- chr_dat %>%
        mutate(X3 = NULL, d = NULL, lw = NULL, ld = NULL, pred = NULL)
    return(chr_dat)
}

# dat_file = "/storage/mathelierarea/processed/vipin/group/HiC_data/HMEC/HMEC/100kb/chr22.txt"
dat_file = "./data/hi-C/100Kb/chr22.txt"
data = load_sparse_data(dat_file)
data = compute_zscore(dat_file,)

# save the table
write.table(data, "./chr22_100kb.tsv", sep = "\t", row.names = F, col.names = F)
