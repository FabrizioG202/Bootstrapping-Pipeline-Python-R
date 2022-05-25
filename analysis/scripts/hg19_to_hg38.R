### TRANSLATES IN PLACE ###
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(valr))

options <- commandArgs(trailingOnly = TRUE)
input_path <- options[[1]]
output_path <- options[[2]]

# the chain
hg19to38_chain <- import.chain("./data/hg19ToHg38.over.chain")

# Read the file.
bed_file <- valr::read_bed(input_path)

# Convert the bed file to a GenomicRanges object
bed_ranges <- GRanges(
    seqnames = bed_file$chrom,
    ranges = IRanges(start = bed_file$start, end = bed_file$end)
)
cat("Size Before Conversion: ", (length(bed_ranges)), "\n")

# Convert the GenomicRanges object to a GenomicRanges object in hg19
bed_file_hg38 <- unlist(liftOver(bed_ranges, hg19to38_chain))

# Save the new file with a hg38 suffix
bed_df <- data.frame(
    chrom = seqnames(bed_file_hg38),
    start = start(bed_file_hg38),
    end = end(bed_file_hg38)
)

# Save the new file to a new file with a hg38 suffix
write.table(bed_df, output_path, sep = "\t", col.names = F, row.names = F, quote = F)
cat("Coordinates translated, new size: ", (length(bed_file_hg38)), "\n")