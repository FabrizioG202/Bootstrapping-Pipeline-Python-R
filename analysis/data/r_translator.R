data = load("./analysis/data/FANTOM5_cage_peak_type_tbl.Rda")

# cage_peak_type_tbl is a tbl_df, save it to a tsv file
write.table(cage_peak_type_tbl, "./analysis/data/FANTOM5_cage_peak_type_tbl.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


data = load("./analysis/data/FANTOM5_CAGE_peak_entrez_gene_tbl.Rda")

# cage_peak_type_tbl is a tbl_df, save it to a tsv file
write.table(FANTOM5_entrez_tbl, "./analysis/data/FANTOM5_CAGE_peak_entrez_gene_tbl.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
