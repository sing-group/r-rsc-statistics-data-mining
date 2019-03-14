source("distance-measures.R")

#	---------------------------------------------------------------------------
#
#	1. Hierarchical Clustering
#
#	---------------------------------------------------------------------------

distances <- jaccard.dist(data$names, toSpectraList(binnedPeaksMatrix))

clustering <- hclust(distances)

clustering.cut <- cutree(clustering, k=3)

#	Run the following commands to see the names of the samples that belong to
#	cluster 1, 2 or 3
	#	which(cutree(clustering, k=3)==1)
	#	which(cutree(clustering, k=3)==2)
	#	which(cutree(clustering, k=3)==3)

png("clustering-jaccard.png", width=1200, height=1200)
plot(clustering, main = "Samples HCA (complete linkage)")
rect.hclust(clustering, k = 3, border = c("blue", "red", "green"))
dev.off()

#	---------------------------------------------------------------------------
#
#	2. Hierarchical Clustering Visualization in Heatmaps
#
#	---------------------------------------------------------------------------

library("gplots")

breaks <- seq(quantile(binnedPeaksMatrix,na.rm=TRUE,probs=0.001),quantile(binnedPeaksMatrix,probs=0.999,na.rm=TRUE),length.out = 24)
spectraColors <- data$spectraColors

png("clustering-jaccard-with-heatmap.png", width=1200, height=1200)
heatmap.2(
	binnedPeaksMatrix,
	Colv=FALSE,
	Rowv=as.dendrogram(clustering),
	dendrogram="row",
	breaks=breaks,
	col=bluered(length(breaks)-1),
	scale="none",
	RowSideColors=spectraColors,
	na.color = 'gray',
	#distfun=hamming.dist,
	#hclustfun=function(d) hclust(d, method="average"),
	rowsep=c(25, 35),
	lwid=c(3,6),
	trace="none",
	margins =c(12,9),	# widens margins around plot
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = c("Healthy", "Lymphoma", "Myeloma"), # category labels
    col = unique(data$spectraColors),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()
