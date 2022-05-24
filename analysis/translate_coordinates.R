library(rtracklayer)
library(GenomicRanges)
library(valr)
path <- "./data/annotations/chr10/chr10_3utr.BED"

hg38to19_chain <- import.chain("./analysis/data/hg38ToHg19.over.chain")
hg19to38_chain <- import.chain("./analysis/data/hg19ToHg38.over.chain")

bed_file <- valr::read_bed(path)

# Convert the bed file to a GenomicRanges object
bed_ranges <- GRanges(
    seqnames = bed_file$chrom,
    ranges = IRanges(start = bed_file$start, end = bed_file$end)
    )

bed_file_hg38 <- (liftOver(bed_ranges, hg38to19_chain))
