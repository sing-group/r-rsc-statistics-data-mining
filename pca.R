library(pca3d)

imagesDirectory <- "images/pca/"

dir.create(imagesDirectory,
           recursive = TRUE,
           showWarnings = FALSE)

source("load-cancer.R")

#	---------------------------------------------------------------------------
#
#	1. Principal Components Analysis
#
#	---------------------------------------------------------------------------
#
#	1.1 Principal Components Analysis using the replicas matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

pca <- prcomp(binnedPeaksMatrix, scale = TRUE)

summary(pca)

# reduce dataset taking the first PCs that accumulate the 95% of variance
pc95 <- min(which(cumsum(pca$sdev ^ 2) / sum(pca$sdev ^ 2) >= 0.95))
binnedPeaksMatrix.pca.95 <- pca$x[, 1:pc95]

png(
  paste0(imagesDirectory, "pca-replicates.png"),
  width = 640,
  height = 480
)
pca2d(
  pca,
  group = data$spectraConditions,
  components = c(1, 2),
  legend = "bottomright"
)
dev.off()

# Uncomment the following line to show the 3D visualiation of the PCA
# pca3d(pca, group=data$spectraConditions, components=c(1,2,3), legend="bottomright")

#	---------------------------------------------------------------------------
#
#	1.2 Principal Components Analysis using the samples matrix (consensusBinnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

pca <- prcomp(consensusBinnedPeaksMatrix, scale = TRUE)
summary(pca)

# reduce dataset taking the first PCs that accumulate the 95% of variance
pc95 <- min(which(cumsum(pca$sdev ^ 2) / sum(pca$sdev ^ 2) >= 0.95))
binnedPeaksMatrix.pca.95 <- pca$x[, 1:pc95]

png(paste0(imagesDirectory, "pca-samples.png"),
    width = 640,
    height = 480)
pca2d(
  pca,
  group = data$samplesConditions,
  components = c(1, 2),
  legend = "bottomright"
)
dev.off()

#with biplot
png(
  paste0(imagesDirectory, "pca-samples-biplot.png"),
  width = 640,
  height = 480
)
pca2d(
  pca,
  group = data$samplesConditions,
  components = c(1, 2),
  legend = "bottomright",
  biplot = TRUE
)
dev.off()

# Uncomment the following line to show the 3D visualiation of the PCA
# pca3d(pca, group=data$samplesConditions, components=c(1,2,3), legend="bottomright")
