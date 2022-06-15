### Contains the common functions which are shared between the different scripts.

# Loads the data from the given path and returns the object contained in the file.
load_data <- function(filepath) {
  dat.out <- get(load(filepath))
  tmp_obj <- names(mget(load(filepath)))
  rm(list = tmp_obj)
  rm(tmp_obj)
  return(dat.out)
}
