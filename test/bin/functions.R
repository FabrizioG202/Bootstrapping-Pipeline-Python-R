################################
## Export the object as RData ##
################################
library(data.table)

create_output_path <- function(path) {
  parts <- unlist(strsplit(path, ".", fixed = T))
  parts[length(parts)] <- "RData"
  output_file <- paste(parts, collapse = ".")
  return(output_file)
}

convert_file <- function(inputfile) {
  if (file.exists(inputfile) & file.size(inputfile) > 0) {
    table = fread(input.file)
    out_path <- create_output_path(inputfile)
    save(table, file = out_path)
  } else {
    message("Input file does not exist or is empty!!")
    stop("")
  }
}