library("MALDIquant")
library("MALDIquantForeign")
library("MALDIquantExamples")

spectra <- import(getPathSpecies(), verbose = FALSE)

spectra <- trim(spectra)

set.seed(2019)

spectra <- transformIntensity(spectra, method = "sqrt")

spectra <-
  smoothIntensity(spectra, method = "SavitzkyGolay", halfWindowSize = 10)

spectra <- removeBaseline(spectra, method = "SNIP", iterations = 25)

spectra <- calibrateIntensity(spectra, method = "TIC")

spectra <- alignSpectra(spectra)

metaData(spectra[[1]])$spot

spots <- sapply(spectra, function(x)
  metaData(x)$spot)
species <- sapply(spectra, function(x)
  metaData(x)$sampleName)

avgSpectra <-
  averageMassSpectra(spectra, labels = paste0(species, spots))

peaks <- detectPeaks(avgSpectra, SNR = 2, halfWindowSize = 10)

peaks.binned <- binPeaks(peaks)

peaks.binned.filtered <-
  filterPeaks(peaks.binned, minFrequency = 0.25)

spots <- sapply(avgSpectra, function(x)
  metaData(x)$spot)
species <- sapply(avgSpectra, function(x)
  metaData(x)$sampleName)
species <- factor(species)

#binnedPeaksMatrix <- intensityMatrix(peaks.binned.filtered, avgSpectra)
binnedPeaksMatrix <- intensityMatrix(peaks.binned.filtered)
rownames(binnedPeaksMatrix) <- paste(species, spots, sep = ".")

binnedPeaksMatrix.conditions <- species
