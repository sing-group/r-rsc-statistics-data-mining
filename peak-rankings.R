imagesDirectory <- "images/important-variables/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	1 Fold-changes
#
#	---------------------------------------------------------------------------

source("load-cancer.R")

source("peak-rankings-functions.R")

avgPeaks <- meanConditionSpectra(consensusData$datasetConditions, consensusBinnedPeaksMatrix, consensusData$spectraConditions)

# log2(avgPeaks["HEALTHY", "2984.26172166667"] / avgPeaks["LYMPHOMA", "2984.26172166667"])
# log2(avgPeaks["HEALTHY", "2984.26172166667"] / avgPeaks["MYELOMA", "2984.26172166667"])

comparison <- compareConditions(avgPeaks)

comparison.log2 <- log2(comparison)

comparison.index <- 1

data <- sort(comparison.log2[,comparison.index])
data.plot <- data[!is.infinite(data) & abs(data) > 1.5]

png(paste0(imagesDirectory, "peak-ranking-fold-changes.png"), width = 1200, height = 800)
barplot(data.plot, horiz=TRUE, names=substr(names(data.plot), 1, 8), las=1, space=1.5, cex.names=0.8, xlim=c(min(data.plot)-1, max(data.plot)+1), main = colnames(comparison)[comparison.index])
dev.off()

#	---------------------------------------------------------------------------
#
#	2 t-scores
#
#	---------------------------------------------------------------------------

source("load-cancer.R")

source("distance-measures.R")

#	---------------------------------------------------------------------------
#
#	2.1 Ranking features using the sda.ranking function from the SDA package
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
#	2.2 Using the 'X' top peaks to perform a Hierarchical Clustering Analysis 
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
