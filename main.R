# This script finish the coding test.
# Author: Tian

source("./AmpliconAnalysis.R")

amplicon1Result <- AmpliconAnalysis("./Data/amplicon_1.fa.gz", 
                                    "./Data/genome.fa.gz", 
                                    "./Data/gene1.bed", 
                                    "./Data/amplicon_coordinates.bed", 
                                    "amplicon_1")

amplicon2Result <- AmpliconAnalysis("./Data/amplicon_2.fa.gz", 
                                    "./Data/genome.fa.gz", 
                                    "./Data/gene1.bed", 
                                    "./Data/amplicon_coordinates.bed", 
                                    "amplicon_2")

amplicon3Result <- AmpliconAnalysis("./Data/amplicon_3.fa.gz",
                                    "./Data/genome.fa.gz",
                                    "./Data/gene2.bed",
                                    "./Data/amplicon_coordinates.bed",
                                    "amplicon_3")

# Export Table for amplicon_1, amplicon_2, and amplicon_3
resultTable <- cbind(amplicon1Result, amplicon2Result, amplicon3Result)
write.table(resultTable, "./ampliconMutation.csv", quote=F, row.names=T)

# Run amplicon_3.test_2.fa
amplicon3_2 <- AmpliconAnalysis("./Data/amplicon_3.test_2.fa.gz",
                                "./Data/genome.fa.gz",
                                "./Data/gene2.bed",
                                "./Data/amplicon_coordinates.bed",
                                "amplicon_3")


print(knitr::kable(do.call("rbind", amplicon3_2[c(2,3,4,6,7)])))

# Draw barplot for mutation frequency
df <- unlist(resultTable["mutationFrequency",])
names(df) <- paste0("amplicon_", 1:3)
xx <- barplot(df, border=F, width=0.7, col="steelblue", xlab="Amplicons", ylab="Mutation Frequency")
text(x = xx, y = df * 0.9, label = resultTable["natureMutation",], pos = 3, cex = 1.2, col = "white")


