# This script identify mutation from a list of amplicon sequence.
# Author: Tian

amplicon_csv <- "../Data/amplicon_1.fa.gz"

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
    Frequency <- sapply(mutationBase, function(x) table(ampliconMatrix[,x]))

    return(list(mutationBase=mutationBase, Frequency=Frequency))
}


result <- MutationDetection(amplicon_csv)

