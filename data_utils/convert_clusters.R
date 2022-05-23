install.packages("rjson")
library(rjson)

directory <- "./data/rda_clusters/H1/"
directory_out <- "./data/clusters/H1/"

for (i in 1:22) {
    path <- paste0(directory, "chr", i, "_spec_res.Rda")
    # save the content of the file in a variable
    file_content <- get(base::load(path))

    da <- rjson::toJSON(file_content, indent = 0, method = "C")

    # Save the data in a file
    write(da, paste0(directory_out, "chr", i, "_spec_res.json"))
}
