library("MALDIquant")

toConsensusSpectraData <- function(data,
                                   tolerance = 0.002,
                                   POP = 0.5) {
  consensusSpectra <- list()
  sampleNames <- data$sampleNames
  firstReplicateIndex <- 1
  
  for (i in 1:length(sampleNames)) {
    currentSample <- sampleNames[i]
    currentSampleIndexes <-
      which(data$spectraSampleNames == currentSample)
    numSamples <- length(currentSampleIndexes)
    sampleSpectra <- data$spectra[currentSampleIndexes]
    
    binnedPeaks <- binPeaks(sampleSpectra, tolerance = tolerance)
    
    binnedPeaksMatrix <- intensityMatrix(binnedPeaks)
    
    
    consensusMasses <- vector()
    consensusIntensities <- vector()
    
    for (j in 1:ncol(binnedPeaksMatrix)) {
      consensusPeaks <- which(binnedPeaksMatrix[, j] != 0)
      currentPresence <- length(consensusPeaks) / numSamples
      if (currentPresence >= POP) {
        consensusIntensities <-
          c(consensusIntensities, mean(binnedPeaksMatrix[consensusPeaks, j]))
        consensusMasses <-
          c(consensusMasses, colnames(binnedPeaksMatrix)[j])
      }
    }
    
    currentConsensusSpectrum <- createMassPeaks(
      mass = as.numeric(consensusMasses),
      intensity = as.numeric(consensusIntensities),
      metaData = list(name = currentSample)
    )
    
    consensusSpectra <-
      c(consensusSpectra, currentConsensusSpectrum)
  }
  
  list(
    spectraNames = sampleNames,
    spectra = consensusSpectra,
    spectraColors = data$samplesColors,
    spectraConditions = data$samplesConditions,
    datasetConditions = data$datasetConditions,
    datasetConditionsColors = data$datasetConditionsColors
  )
}

asPresenceMatrix <- function(data,
                             present = 1,
                             absent = 0) {
  toret <- data
  toret[data > 0] <- present
  
  toret[is.na(data)] <- absent
  
  toret
}
