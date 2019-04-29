imagesDirectory <- "images/clustering/kmeans/"

dir.create(imagesDirectory,
           recursive = TRUE,
           showWarnings = FALSE)

source("load-cancer.R")

#	---------------------------------------------------------------------------
#
#	1 K-means clustering
#
#	---------------------------------------------------------------------------

# Since k-means starts its centroids randomnly, we set a random seed in order
# to reproduce the same results

set.seed(2019)

#	---------------------------------------------------------------------------
#
#	1.1 K-means using the replicates matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

clustering.kmeans <- kmeans(binnedPeaksMatrix, centers = 3)
clustering.kmeans$cluster

#	---------------------------------------------------------------------------
#
#	1.2 K-means using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

clustering.kmeans.consensus <-
  kmeans(consensusBinnedPeaksMatrix, centers = 3)
clustering.kmeans.consensus$cluster

#	---------------------------------------------------------------------------
#
#	1.3 K-means using the replicates matrix (binnedPeaksMatrix) and using presence/absence
#
#	---------------------------------------------------------------------------

presenceBinnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix)
clustering.kmeans.presence <- kmeans(binnedPeaksMatrix, centers = 3)
clustering.kmeans.presence$cluster

#	---------------------------------------------------------------------------
#
#	1.4 K-means using the samples matrix (consensusBinnedPeaksMatrix) and using presence/abscence
#
#	---------------------------------------------------------------------------

presenceConsensusBinnedPeaksMatrix <-
  asPresenceMatrix(consensusBinnedPeaksMatrix)
clustering.kmeans.presence.consensus <-
  kmeans(consensusBinnedPeaksMatrix, centers = 3)
clustering.kmeans.presence.consensus$cluster


#	---------------------------------------------------------------------------
#
#	1.5 Try different number of clusters and plot the total withinss (calculated
#       with sum(withinss) or using tot.withinss
#
#	---------------------------------------------------------------------------

clustering.kmeans.k_withinss <- vector()
for (c in 1:10) {
  set.seed(2019)
  clust <- kmeans(binnedPeaksMatrix, centers = c)
  clustering.kmeans.k_withinss[c] <- sum(clust$withinss)
}
png(
  paste0(imagesDirectory, "clustering-kmeans-k-withinss.png"),
  width = 640,
  height = 480
)
plot(
  clustering.kmeans.k_withinss,
  type = "b",
  xlab = "Number of clusters",
  ylab = "total within_SS",
  main = "K-means quality at different k (replicates matrix)"
)
dev.off()

clustering.kmeans.k_withinss <- vector()
for (c in 1:10) {
  set.seed(2019)
  clust <- kmeans(consensusBinnedPeaksMatrix, centers = c)
  clustering.kmeans.k_withinss[c] <- sum(clust$withinss)
}
png(
  paste0(imagesDirectory, "clustering-kmeans-samples-k-withinss.png"),
  width = 640,
  height = 480
)
plot(
  clustering.kmeans.k_withinss,
  type = "b",
  xlab = "Number of clusters",
  ylab = "total within_SS",
  main = "K-means quality at different k (samples matrix)"
)
dev.off()

#	---------------------------------------------------------------------------
#
#	1.6 Get a consensus clustering with ConsensusClusterPlus
#
#	---------------------------------------------------------------------------

library("ConsensusClusterPlus")

pdf(paste0(imagesDirectory, "consensus-cluster-plus.pdf"))
ccp <-
  ConsensusClusterPlus(
    t(binnedPeaksMatrix),
    maxK = 3,
    distance = "euclidean",
    clusterAlg = "km",
    reps = 50,
    seed = 2019
  )
dev.off()
# show clustering results for k=3
ccp[[3]][["consensusClass"]]
