library(rtracklayer)
library(GenomicRanges)
library(valr)
path <- "./data/enhancers_in_clusters.tsv"

hg38to19_chain <- import.chain("./analysis/data/hg38ToHg19.over.chain")
hg19to38_chain <- import.chain("./analysis/data/hg19ToHg38.over.chain")

# Read the file.
bed_file <- valr::read_bed(path)

# Convert the bed file to a GenomicRanges object
bed_ranges <- GRanges(
    seqnames = bed_file$chrom,
    ranges = IRanges(start = bed_file$start, end = bed_file$end)
    )

# Convert the GenomicRanges object to a GenomicRanges object in hg19
bed_file_hg38 <- unlist(liftOver(bed_ranges, hg19to38_chain))
bed_file_hg38

#Save the new file with a hg38 suffix
bed_df = data.frame(
    chrom = seqnames(bed_file_hg38),
    start = start(bed_file_hg38),
    end = end(bed_file_hg38)
)

# Save the new file to a new file with a hg38 suffix
write.table(bed_df, "./data/enhancers_in_clusters_hg38.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
