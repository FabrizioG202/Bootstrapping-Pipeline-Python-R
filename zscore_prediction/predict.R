library(mgcv)
library(tidyverse)

# Load Data
inside = read_table("./results/hubs_data/100kb/inside.tsv", col_names =FALSE)
outside = read_table("./results/hubs_data/100kb/outside.tsv", col_names =FALSE)

# Reading the outside and inside data.
outside = outside %>% mutate(X1 = NULL, X2 = NULL) %>% rename("zscore" = X3, "pvalue" = X4)
inside = inside %>% mutate(X1 = NULL, X2 = NULL) %>% rename("zscore" = X3, "pvalue" = X4)

# Converting the data to dataframe
outside = as.data.frame(outside)
inside = as.data.frame(inside)

# Creating the model
inside_model = gam(zscore ~ s(pvalue), data = inside)
outside_model = gam(zscore ~ s(pvalue), data = outside)

library(MASS)
mcycle <- MASS::mcycle

# Load mgcv
library(mgcv)