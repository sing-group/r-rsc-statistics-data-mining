source("load-cancer.R")
source("distance-measures.R")

imagesDirectory <- "images/important-variables/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	1. Ranking features using the sda.ranking function from the SDA package
#
#		Homepage: http://strimmerlab.org/software/sda/
#		CRAN: http://cran.r-project.org/package=sda
#
#	---------------------------------------------------------------------------

library("sda")

ddar <- sda.ranking(Xtrain=consensusBinnedPeaksMatrix, L=consensusData$spectraConditions, fdr=FALSE, diagonal=TRUE)
png(paste0(imagesDirectory, "peak-ranking-t-scores.png"), width = 1200, height = 800)
plot(ddar, top=20, arrow.col="blue", zeroaxis.col="red", ylab="Peaks")
dev.off()

#	---------------------------------------------------------------------------
#
#	2. Using the 'X' top peaks to perform a Hierarchical Clustering Analysis 
#       (using Jaccard distance)
#
#	---------------------------------------------------------------------------

consensusBinnedPeaksMatrix.topPeaks <- consensusBinnedPeaksMatrix[, sort(ddar[1:20,"idx"])]

distances.jaccard.consensus <- jaccard.dist(rownames(consensusBinnedPeaksMatrix.topPeaks), toSpectraList(consensusBinnedPeaksMatrix.topPeaks))

clustering.jaccard.consensus <- hclust(distances.jaccard.consensus)

png(paste0(imagesDirectory, "peak-ranking-t-scores-hierarchical-clustering.png"), width = 1200, height = 800)
plot(clustering.jaccard.consensus, main = "Samples HCA (complete linkage)")
rect.hclust(clustering.jaccard.consensus, k = 3, border = c("blue", "red", "green"))
dev.off()
