# This script will extract sequence from Reference fasta file.
# Author: Tian

library("seqinr")

getSeqFromReference <- function(refPath, chr, start, end) {

    # Read fastq into R session.
    ref <- read.fasta(file = refPath)

    # Select matched chromosome.
    idx <- which(names(ref) == chr)

    # Select start/end position.
    toupper(paste(ref[[idx]][start:end], collapse=""))
}

RefSeq <- getSeqFromReference("../Data/genome.fa.gz", "chr7", 55181321, 55181390)

