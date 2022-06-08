library(testthat)

source("functions.R")

test_that("output extension", {
    test_input <- "bed.1.tsv"
    expected_output <- "bed.1.RData"
    actual_output <- create_output_path(test_input)
    expect_equal(actual_output, expected_output)
})

# Check if the file exist
test_that("file exists", {
	test_input <- "bed.1.tsv"
	expect_error(convert_file(test_input))
})
