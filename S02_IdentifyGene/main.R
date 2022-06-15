# This script identify genes each amplicon file belongs to.
# Author: Tian

library("GenomicRanges")

geneBed <- "../Data/gene1.bed"
coordinateBed <- "../Data/amplicon_coordinates.bed"

IdentifyGene <- function(geneBed, coordinateBed) {
    
    # Read into Gene bed file and Coordiniates bed file.
    Genes <- read.csv(geneBed, sep="\t", header=F, row.names=1)
    Coordiniates <- read.csv(coordinateBed, sep="\t", header=F, row.names=4)
    
    # Organize data to fit GRange transfer.
    Genes$V2 <- paste0("chr", Genes$V2)
    colnames(Genes) <- c("seqnames", "start", "end", "width")
    colnames(Coordiniates) <- c("seqnames", "start", "end")
    
    # Use findOverlaps to find the gene.
    ov <- as.data.frame(findOverlaps(makeGRangesFromDataFrame(Coordiniates), makeGRangesFromDataFrame(Genes)))
    
    return(Genes[unique(ov[,2]),])
}

geneName <- IdentifyGene(geneBed, coordinateBed)

