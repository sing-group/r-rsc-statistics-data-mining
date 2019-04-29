source("load-cancer.R")

imagesDirectory <- "images/distance-measures/"

dir.create(imagesDirectory,
           recursive = TRUE,
           showWarnings = FALSE)

#	---------------------------------------------------------------------------
#
#	1. Pairwise Comparisons: calculating the distances between two spectra
#
#	---------------------------------------------------------------------------
#
#	1.1 Euclidean Distance
#
#	---------------------------------------------------------------------------

spectrum_1 <- binnedPeaksMatrix[1, ]
spectrum_2 <- binnedPeaksMatrix[2, ]
dist(rbind(spectrum_1, spectrum_2), method = "euclidean")

#	---------------------------------------------------------------------------
#
#	1.2 Manhattan Distance
#
#	---------------------------------------------------------------------------

presenceBinnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix)
spectrum_1 <- presenceBinnedPeaksMatrix[1, ]
spectrum_2 <- presenceBinnedPeaksMatrix[2, ]
dist(rbind(spectrum_1, spectrum_2), method = "manhattan")

#	---------------------------------------------------------------------------
#
#	1.3 Jaccard Distance
#
#	---------------------------------------------------------------------------

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

all_masses <- colnames(binnedPeaksMatrix)
masses_1 <- all_masses[binnedPeaksMatrix[1, ] != 0]
masses_2 <- all_masses[binnedPeaksMatrix[2, ] != 0]
jaccard(masses_1, masses_2)

#	---------------------------------------------------------------------------
#
#	1.4 Pearson Correlation Distance
#
#	---------------------------------------------------------------------------

spectrum_1 <- binnedPeaksMatrix[1, ]
spectrum_2 <- binnedPeaksMatrix[2, ]
1 - cor(spectrum_1, spectrum_2, method = "pearson")

#	---------------------------------------------------------------------------
#
#	2. Creating Distance Matrices
#
#	---------------------------------------------------------------------------
#
#	2.1 Definition of a function using ggplot2 to visualize the distances
#       matrices.
#
#	---------------------------------------------------------------------------

library("ggplot2")

plotDistances <- function(distances, binnedPeaksMatrix) {
  distanceMatrix <- as.matrix(distances)
  sampleNames <- rownames(binnedPeaksMatrix)
  dataPlot <-
    data.frame(row.names = seq(length(sampleNames) * length(sampleNames)))
  dataPlot$sampleA <- rep(sampleNames, each = length(sampleNames))
  dataPlot$sampleB <- rep(sampleNames, length(sampleNames))
  dataPlot$value <- as.vector(t(distanceMatrix))
  
  
  ggplot(data = dataPlot,
         aes(x = sampleA, y = sampleB, fill = value)) + geom_tile() +
    theme_bw() + xlab("sample") + ylab("sample") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

#	---------------------------------------------------------------------------
#
#	2.2 Euclidean Distance
#
#	---------------------------------------------------------------------------

distances <- dist(binnedPeaksMatrix, method = "euclidean")

png(paste0(imagesDirectory, "euclidean.png"),
    width = 1200,
    height = 1200)
plotDistances(distances, binnedPeaksMatrix)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.3 Manhattan Distance
#
#	---------------------------------------------------------------------------

presenceBinnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix)
distances <- dist(presenceBinnedPeaksMatrix, method = "manhattan")

png(paste0(imagesDirectory, "manhattan.png"),
    width = 1200,
    height = 1200)
plotDistances(distances, presenceBinnedPeaksMatrix)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.4 Pearson Correlation Distance
#
#	---------------------------------------------------------------------------

pearson.dist <- function(binnedPeaksMatrix) {
  result <-
    1 - cor(t(binnedPeaksMatrix), t(binnedPeaksMatrix), method = "pearson")
  
  as.dist(result)
}

distances <- pearson.dist(binnedPeaksMatrix)

png(paste0(imagesDirectory, "pearson.png"),
    width = 1200,
    height = 1200)
plotDistances(distances, binnedPeaksMatrix)
dev.off()

#	---------------------------------------------------------------------------
#
#	2.5 Jaccard Distance
#
#	---------------------------------------------------------------------------

toSpectraList <- function(binnedPeaksMatrix) {
  spectraData <- list()
  
  for (i in 1:nrow(binnedPeaksMatrix)) {
    intensities <- as.numeric(binnedPeaksMatrix[i, ])
    masses <-
      as.numeric(as.character(colnames(binnedPeaksMatrix)))
    
    masses <- masses[intensities > 0]
    intensities <- intensities[intensities > 0]
    
    spectrum <-
      createMassPeaks(
        mass = masses,
        intensity = intensities,
        metaData = list(name = rownames(binnedPeaksMatrix)[i])
      )
    spectraData <- c(spectraData, spectrum)
  }
  spectraData
}

jaccard.dist <- function(dataNames, spectraData) {
  result <- matrix(nrow = length(dataNames), ncol = length(dataNames))
  colnames(result) <- dataNames
  rownames(result) <- dataNames
  
  for (i in 1:length(spectraData)) {
    currentValues <- rep(NA, length(dataNames))
    for (j in 1:length(spectraData)) {
      massesA <- mass(spectraData[[i]])
      massesB <- mass(spectraData[[j]])
      currentValues[j] <- 1 - jaccard(massesA, massesB)
    }
    result[i, ] <- currentValues
  }
  
  as.dist(result)
}

distances <-
  jaccard.dist(rownames(binnedPeaksMatrix),
               toSpectraList(binnedPeaksMatrix))

png(paste0(imagesDirectory, "jaccard.png"),
    width = 1200,
    height = 1200)
plotDistances(distances, binnedPeaksMatrix)
dev.off()
