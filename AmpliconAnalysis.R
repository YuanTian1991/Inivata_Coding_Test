# This is a function to read into amplicon sequence file, reference sequence, gene bed file, to return:
# 1. gene name
# 2. nature of the mutation (e.g. C->T etc)
# 3. frequency of the mutation.
# 4. number of amplicons supporting the mutation
# 5. reference sequence for the amplicon
# 6. mutated sequence for the amplicon

# Author: Tian

library("GenomicRanges")
library("seqinr")

# Example Input:
# ampliconFile <- "./Data/amplicon_3.test_2.fa.gz"
# referenceFile <- "./Data/genome.fa.gz"
# geneBed <- "./Data/gene2.bed"
# coordinateBed <- "./Data/amplicon_coordinates.bed"
# ampliconID <- "amplicon_3"

AmpliconAnalysis <- function(ampliconFile, referenceFile, geneBed, coordinateBed, ampliconID) {
    Genes <- read.csv(geneBed, sep="\t", header=F, row.names=1)
    Coordiniates <- read.csv(coordinateBed, sep="\t", header=F, row.names=4)
    
    IdentifyGene <- function(genes, coordinate) {
        # Organize data to fit GRange transfer.
        genes$V2 <- paste0("chr", genes$V2)
        colnames(genes) <- c("seqnames", "start", "end", "width")
        colnames(coordinate) <- c("seqnames", "start", "end")
        # Use findOverlaps to find the gene.
        ov <- as.data.frame(findOverlaps(makeGRangesFromDataFrame(coordinate), makeGRangesFromDataFrame(genes)))
        return(rownames(genes)[unique(ov[,2])])
    }
    
    getSeqFromReference <- function(refPath, chr, start, end) {
        # Read fastq into R session.
        ref <- read.fasta(file = refPath)
        # Select matched chromosome.
        idx <- which(names(ref) == chr)
        # Select start/end position.
        toupper(paste(ref[[idx]][start:end], collapse=""))
    }
    
    MutationDetection <- function(csv) {
        # Read amplicons into R session.
        amplicons <- read.csv(csv, header=FALSE)[,1]
        amplicons <- amplicons[seq(2, length(amplicons), by=2)]
        
        # Break strings into matrix. In this case I directly use data.frame,
        # if the size is too big, data.table should be considered.
        ampliconMatrix <- t(sapply(amplicons, function(x) strsplit(x, "")[[1]]))
        rownames(ampliconMatrix) <- NULL
        
        # Check each column (base) for mutation detection.
        uniqueBase <- apply(ampliconMatrix, 2, function(x) length(table(x)))
        mutationBase <- which(uniqueBase != 1)
        
        # Count frequency for each base.
        Frequency <- list()
        for(i in 1:length(mutationBase)) {
            Frequency[[i]] <- table(ampliconMatrix[, mutationBase[i]])
        }
        return(list(mutationBase=mutationBase, Frequency=Frequency))
    }
    
    geneName <- IdentifyGene(Genes, Coordiniates[ampliconID, ])
    referenceSequence <- getSeqFromReference(referenceFile, 
                                             Coordiniates[ampliconID, "V1"],
                                             Coordiniates[ampliconID, "V2"],
                                             Coordiniates[ampliconID, "V3"])

    mutationResult <- MutationDetection(ampliconFile)
    natureMutation <- c()
    mutationAmpliconNumber <- c()
    mutationFrequency <- c()
    mutationSequence <- c()
    mBases <- c()
    
    for(i in 1:length(mutationResult$mutationBase)) {
        mBase <- mutationResult$mutationBase[i]
        mBases <- c(mBases, mBase)

        natureBase <- substr(referenceSequence, mBase, mBase)
        mutationBase <- setdiff(names(mutationResult$Frequency[[i]]), natureBase)
        
        natureMutation <- c(natureMutation, paste0(natureBase, "->", mutationBase))
        mutationAmpliconNumber <- c(mutationAmpliconNumber, mutationResult$Frequency[[i]][mutationBase])
        mutationFrequency <- c(mutationFrequency, mutationResult$Frequency[[i]][mutationBase] / sum(mutationResult$Frequency[[i]]))
        
        tmpMutationSequence <- referenceSequence
        substr(tmpMutationSequence, mBase, mBase) <- mutationBase
        mutationSequence <- c(mutationSequence, tmpMutationSequence)
    }

    return(list(geneName=geneName,
                natureMutation=natureMutation,
                mutationFrequency=mutationFrequency,
                mutationAmpliconNumber=mutationAmpliconNumber,
                referenceSequence=referenceSequence,
                mutationSequence=mutationSequence,
                # Below are extra output.
                mBases=mBases
                ))
}


