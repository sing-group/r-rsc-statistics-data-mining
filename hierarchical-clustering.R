source("distance-measures.R")

imagesDirectory <- "images/clustering/hierarchical/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	1. Hierarchical Clustering
#
#	---------------------------------------------------------------------------
#
#	1.1.1 Jaccard using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

distances.jaccard <- jaccard.dist(rownames(binnedPeaksMatrix), toSpectraList(binnedPeaksMatrix))

clustering.jaccard <- hclust(distances.jaccard)

clustering.jaccard.cut <- cutree(clustering.jaccard, k=3)

#	Run the following commands to see the names of the samples that belong to
#	cluster 1, 2 or 3
	#	which(cutree(clustering, k=3)==1)
	#	which(cutree(clustering, k=3)==2)
	#	which(cutree(clustering, k=3)==3)

png(paste0(imagesDirectory, "clustering-jaccard-replicates.png"), width=1200, height=1200)
plot(clustering.jaccard, main = "Replicates HCA (complete linkage)")
rect.hclust(clustering.jaccard, k = 3, border = c("blue", "red", "green"))
dev.off()

#	---------------------------------------------------------------------------
#
#	1.1.2 Jaccard using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

distances.jaccard.consensus <- jaccard.dist(rownames(consensusBinnedPeaksMatrix), toSpectraList(consensusBinnedPeaksMatrix))

clustering.jaccard.consensus <- hclust(distances.jaccard.consensus)

png(paste0(imagesDirectory, "clustering-jaccard-samples.png"), width=1200, height=1200)
plot(clustering.jaccard.consensus, main = "Samples HCA (complete linkage)")
rect.hclust(clustering.jaccard.consensus, k = 3, border = c("blue", "red", "green"))
dev.off()

#	---------------------------------------------------------------------------
#
#	1.2.1 Manhattan using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

presenceBinnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix)

distances.manhattan <- dist(presenceBinnedPeaksMatrix, method="manhattan")

clustering.manhattan <- hclust(distances.manhattan)

clustering.manhattan.cut <- cutree(clustering.manhattan, k=3)

png(paste0(imagesDirectory, "clustering-manhattan-replicates.png"), width=1200, height=1200)
plot(clustering.manhattan, main = "Replicates HCA (complete linkage)")
rect.hclust(clustering.manhattan, k = 3, border = c("blue", "red", "green"))
dev.off()

#	---------------------------------------------------------------------------
#
#	1.2.2 Manhattan using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

presenceConsensusBinnedPeaksMatrix <- asPresenceMatrix(consensusBinnedPeaksMatrix)

distances.manhattan.consensus <- dist(presenceConsensusBinnedPeaksMatrix, method="manhattan")

clustering.manhattan.consensus <- hclust(distances.manhattan.consensus)

png(paste0(imagesDirectory, "clustering-manhattan-samples.png"), width=1200, height=1200)
plot(clustering.manhattan.consensus, main = "Samples HCA (complete linkage)")
rect.hclust(clustering.manhattan.consensus, k = 3, border = c("blue", "red", "green"))
dev.off()

#	---------------------------------------------------------------------------
#
#	1.3.1 Euclidean using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

distances.euclidean <- dist(binnedPeaksMatrix, method="euclidean")

clustering.euclidean <- hclust(distances.euclidean)

clustering.euclidean.cut <- cutree(clustering.euclidean, k=3)

png(paste0(imagesDirectory, "clustering-euclidean-replicates.png"), width=1200, height=1200)
plot(clustering.euclidean, main = "Replicates HCA (complete linkage)")
dev.off()

#	---------------------------------------------------------------------------
#
#	1.3.2 Euclidean using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

distances.euclidean.consensus <- dist(consensusBinnedPeaksMatrix, method="euclidean")

clustering.euclidean.consensus <- hclust(distances.euclidean.consensus)

png(paste0(imagesDirectory, "clustering-euclidean-samples.png"), width=1200, height=1200)
plot(clustering.euclidean.consensus, main = "Samples HCA (complete linkage)")
dev.off()

#	---------------------------------------------------------------------------
#
#	1.4.1 Correlation using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

distances.correlation <- pearson.dist(binnedPeaksMatrix)

clustering.correlation <- hclust(distances.correlation)

clustering.correlation.cut <- cutree(clustering.correlation, k=3)

png(paste0(imagesDirectory, "clustering-correlation-replicates.png"), width=1200, height=1200)
plot(clustering.correlation, main = "Replicates HCA (complete linkage)")
dev.off()

#	---------------------------------------------------------------------------
#
#	1.4.2 Correlation using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

distances.correlation.consensus <- pearson.dist(consensusBinnedPeaksMatrix)

clustering.correlation.consensus <- hclust(distances.correlation.consensus)

png(paste0(imagesDirectory, "clustering-correlation-samples.png"), width=1200, height=1200)
plot(clustering.correlation.consensus, main = "Samples HCA (complete linkage)")
dev.off()

#	---------------------------------------------------------------------------
#
#	2. Hierarchical Clustering Visualization With Heatmaps
#
#	---------------------------------------------------------------------------
#
#	2.1 Reusable "plotHeatmap" Function
#
#	---------------------------------------------------------------------------

library("gplots")

plotHeatmap <- function(peaksMatrix, clustering, spectraColors, breaksCount = 24, rowsep = c(), datasetConditions = c(), datasetConditionsColors = c()) {
	breaks <- seq(
		quantile(peaksMatrix, na.rm=TRUE, probs=0.001),
		quantile(peaksMatrix, na.rm=TRUE, probs=0.999),
		length.out = breaksCount
	)

	rowSideColors <- spectraColors
	
	
	heatmap.2(
		peaksMatrix,
		Colv=FALSE,
		Rowv=as.dendrogram(clustering),
		dendrogram="row",
		breaks=breaks,
		col=bluered(length(breaks)-1),
		scale="none",
		RowSideColors=rowSideColors,
		na.color = "gray",
		rowsep=rowsep,
		lwid=c(3,6),
		trace="none",
		margins =c(12,9),
	)
	if(length(datasetConditions)>0 && length(datasetConditionsColors)>0){
		par(lend = 1)           # square line ends for the color legend
		legend("topright",      # location of the legend on the heatmap plot
			legend = datasetConditions, # category labels
			col = datasetConditionsColors,  # color key
			lty= 1,             # line style
			lwd = 10            # line width
		)
	}
}

#	---------------------------------------------------------------------------
#
#	2.2.1 Jaccard using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-jaccard-replicates-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(binnedPeaksMatrix, clustering.jaccard, data$spectraColors, rowsep = c(25, 35), datasetConditions = data$datasetConditions, datasetConditionsColors = data$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.2.2 Jaccard using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-jaccard-samples-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(consensusBinnedPeaksMatrix, clustering.jaccard.consensus, consensusData$spectraColors, rowsep = c(5, 7), datasetConditions = consensusData$datasetConditions, datasetConditionsColors = consensusData$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.3.1 Manhattan using the replicates matrix (presenceBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-manhattan-replicates-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(presenceBinnedPeaksMatrix, clustering.manhattan, data$spectraColors, breaksCount = 4, rowsep = c(25, 35), datasetConditions = data$datasetConditions, datasetConditionsColors = data$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.3.1 Manhattan using the samples matrix (presenceConsensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-manhattan-samples-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(presenceConsensusBinnedPeaksMatrix, clustering.manhattan.consensus, consensusData$spectraColors, breaksCount = 4, rowsep = c(5, 7), datasetConditions = consensusData$datasetConditions, datasetConditionsColors = consensusData$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.4.1 Euclidean using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-euclidean-replicates-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(binnedPeaksMatrix, clustering.euclidean, data$spectraColors, datasetConditions = data$datasetConditions, datasetConditionsColors = data$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.4.2 Euclidean using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-euclidean-samples-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(consensusBinnedPeaksMatrix, clustering.euclidean.consensus, consensusData$spectraColors, rowsep = c(5, 10), datasetConditions = consensusData$datasetConditions, datasetConditionsColors = consensusData$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.5.1 Correlation using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-correlation-replicates-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(binnedPeaksMatrix, clustering.correlation, data$spectraColors, datasetConditions = data$datasetConditions, datasetConditionsColors = data$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.5.2 Correlation using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "clustering-correlation-samples-with-heatmap.png"), width=1200, height=1200)
plotHeatmap(consensusBinnedPeaksMatrix, clustering.correlation.consensus, consensusData$spectraColors, rowsep = c(5, 7), datasetConditions = consensusData$datasetConditions, datasetConditionsColors = data$datasetConditionsColors)
dev.off()

#	---------------------------------------------------------------------------
#
#	3.  Assessing the uncertainty in hierarchical clustering with pvclust
#
#		Homepage: http://stat.sys.i.kyoto-u.ac.jp/prog/pvclust/
#		CRAN: https://cran.r-project.org/web/packages/pvclust/index.html
#		Paper: https://doi.org/10.1093/bioinformatics/btl117
#
#	---------------------------------------------------------------------------

library("pvclust")

# Set a seed so that the pvclust always gives the same results
set.seed(2019)

#	---------------------------------------------------------------------------
#
#	3.1.1 Jaccard using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

jaccard.dist.2 <- function(binnedPeaksMatrix) {
	jaccard.dist(rownames(t(binnedPeaksMatrix)), toSpectraList(t(binnedPeaksMatrix)))
}

res.pv.jaccard <- pvclust(t(binnedPeaksMatrix), method.hclust = "complete", method.dist = jaccard.dist.2, nboot = 1000, parallel=F)

png(paste0(imagesDirectory, "clustering-jaccard-replicates-pvclust.png"), width=1200, height=1200)
plot(res.pv.jaccard, hang = -1, cex = 0.5)
pvrect(res.pv.jaccard)
dev.off()

#	---------------------------------------------------------------------------
#
#	3.1.2 Jaccard using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.jaccard.consensus <- pvclust(t(consensusBinnedPeaksMatrix), method.hclust = "complete", method.dist = jaccard.dist.2, nboot = 1000, parallel=F)

png(paste0(imagesDirectory, "clustering-jaccard-samples-pvclust.png"), width=1200, height=1200)
plot(res.pv.jaccard.consensus, hang = -1, cex = 0.5)
pvrect(res.pv.jaccard.consensus)
dev.off()

#	---------------------------------------------------------------------------
#
#	3.2.1 Manhattan using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.manhattan <- pvclust(t(presenceBinnedPeaksMatrix), method.hclust = "complete", method.dist = "manhattan", nboot = 1000, parallel=T)

png(paste0(imagesDirectory, "clustering-manhattan-replicates-pvclust.png"), width=1200, height=1200)
plot(res.pv.manhattan, hang = -1, cex = 0.5)
pvrect(res.pv.manhattan)
dev.off()

#	---------------------------------------------------------------------------
#
#	3.2.2 Manhattan using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.manhattan.consensus <- pvclust(t(presenceConsensusBinnedPeaksMatrix), method.hclust = "complete", method.dist = "manhattan", nboot = 1000, parallel=T)

png(paste0(imagesDirectory, "clustering-manhattan-samples-pvclust.png"), width=1200, height=1200)
plot(res.pv.manhattan.consensus, hang = -1, cex = 0.5)
pvrect(res.pv.manhattan.consensus)
dev.off()

#	---------------------------------------------------------------------------
#
#	3.3.1 Euclidean using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.euclidean <- pvclust(t(binnedPeaksMatrix), method.hclust = "complete", method.dist = "euclidean", nboot = 1000, parallel=T)

png(paste0(imagesDirectory, "clustering-euclidean-replicates-pvclust.png"), width=1200, height=1200)
plot(res.pv.euclidean, hang = -1, cex = 0.5)
pvrect(res.pv.euclidean)
dev.off()


#	---------------------------------------------------------------------------
#
#	3.3.2 Euclidean using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.euclidean.consensus <- pvclust(t(consensusBinnedPeaksMatrix), method.hclust = "complete", method.dist = "euclidean", nboot = 1000, parallel=T)

png(paste0(imagesDirectory, "clustering-euclidean-samples-pvclust.png"), width=1200, height=1200)
plot(res.pv.euclidean.consensus, hang = -1, cex = 0.5)
pvrect(res.pv.euclidean.consensus)
dev.off()


#	---------------------------------------------------------------------------
#
#	3.4.1 Correlation using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.correlation <- pvclust(t(binnedPeaksMatrix), method.hclust = "complete", method.dist = "correlation", nboot = 1000, parallel=T)

png(paste0(imagesDirectory, "clustering-correlation-replicates-pvclust.png"), width=1200, height=1200)
plot(res.pv.correlation, hang = -1, cex = 0.5)
pvrect(res.pv.correlation)
dev.off()

#	---------------------------------------------------------------------------
#
#	3.4.2 Correlation using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

res.pv.correlation.consensus <- pvclust(t(consensusBinnedPeaksMatrix), method.hclust = "complete", method.dist = "correlation", nboot = 1000, parallel=T)

png(paste0(imagesDirectory, "clustering-correlation-samples-pvclust.png"), width=1200, height=1200)
plot(res.pv.correlation.consensus, hang = -1, cex = 0.5)
pvrect(res.pv.correlation.consensus)
dev.off()
