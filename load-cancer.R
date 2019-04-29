library("MALDIquant")

source("data-functions.R")

# Load all files inside the dataDir directory and returns a list of MassPeaks
# objects.
#
# Files are expected to be comma-separated files with two columns: mass (m/z)
# and intensity.
loadDirectory <- function(dataDir) {
  spectra <- list.files(dataDir)
  spectraData <- list()
  
  for (spectrumFile in spectra) {
    spectrumFile <- paste(dataDir, "/", spectrumFile, sep = "")
    data <- read.csv(spectrumFile)
    spectrum <-
      createMassPeaks(
        mass = data[, 1],
        intensity = data[, 2],
        metaData = list(name = spectrumFile)
      )
    spectraData <- c(spectraData, spectrum)
  }
  
  spectraData
}

loadSamples <- function(dataDir) {
  samples <- list.dirs(dataDir, full.names = FALSE, recursive = FALSE)
  spectra <- list()
  spectraNames <- list()
  spectraConditions <-  list()
  spectraSampleNames <- list()
  
  for (sampleDir in samples) {
    sampleDirectory <- paste(dataDir, "/", sampleDir, sep = "")
    sampleSpectra <- loadDirectory(sampleDirectory)
    spectra <- c(spectra, sampleSpectra)
    spectraNames <-
      c(spectraNames, sapply(
        1:length(sampleSpectra),
        FUN = function(x)
          paste(sampleDir, "_R", x, sep = "")
      ))
    spectraSampleNames <-
      c(spectraSampleNames, rep(sampleDir, length(sampleSpectra)))
    spectraConditions <-
      c(spectraConditions, rep(tail(strsplit(dataDir, "/")[[1]], n = 1), length(sampleSpectra)))
  }
  
  list(
    spectraNames = spectraNames,
    spectraSampleNames = spectraSampleNames,
    spectraConditions = spectraConditions,
    spectra = spectra
  )
}

loadDirectories <- function(dataDirs, col) {
  if (missing(col)) {
    palette <- sample(colors(TRUE))[1:length(dataDirs)]
  } else {
    palette <- col[1:length(dataDirs)]
  }
  
  spectraColors <- list()
  spectraNames <- list()
  spectraSampleNames <- list()
  spectraConditions <- list()
  spectra <- list()
  
  datasetConditions <- vector()
  datasetConditionsColors <- vector()
  samplesColors <- vector()
  samplesConditions <- vector()
  
  i <- 1
  for (dataDir in dataDirs) {
    datasetCondition <- tail(strsplit(dataDir, "/")[[1]], n = 1)
    datasetConditionsColors <-
      c(datasetConditionsColors, palette[i])
    datasetConditions <- c(datasetConditions, datasetCondition)
    data <- loadSamples(dataDir)
    spectraNames <- unlist(c(spectraNames, data$spectraNames))
    spectraSampleNames <-
      unlist(c(spectraSampleNames, data$spectraSampleNames))
    spectra <- c(spectra, data$spectra)
    spectraConditions <-
      unlist(c(spectraConditions, data$spectraConditions))
    spectraColors <-
      unlist(c(spectraColors, rep(palette[i], length(data$spectra))))
    samplesColors <-
      c(samplesColors, rep(palette[i], length(unique(
        data$spectraSampleNames
      ))))
    samplesConditions <-
      c(samplesConditions, rep(data$spectraConditions[[1]], length(unique(
        data$spectraSampleNames
      ))))
    
    i <- i + 1
  }
  
  sampleNames <- unique(spectraSampleNames)
  
  list(
    spectraNames = spectraNames,
    spectraSampleNames = spectraSampleNames,
    spectraConditions = spectraConditions,
    spectra = spectra,
    spectraColors = spectraColors,
    sampleNames = sampleNames,
    samplesColors = samplesColors,
    samplesConditions = samplesConditions,
    datasetConditions = datasetConditions,
    datasetConditionsColors = datasetConditionsColors
  )
}

getBinnedPeaksMatrix <-
  function(data,
           tolerance = 0.002,
           peakIntensityThreshold = 0) {
    binnedPeaks <- binPeaks(data$spectra, tolerance = tolerance)
    binnedPeaksMatrix <- intensityMatrix(binnedPeaks)
    binnedPeaksMatrix[is.na(binnedPeaksMatrix)] <- 0
    binnedPeaksMatrix <-
      binnedPeaksMatrix[, apply(binnedPeaksMatrix, 2, max) >= peakIntensityThreshold]
    rownames(binnedPeaksMatrix) <- data$spectraNames
    
    binnedPeaksMatrix
  }

dataDirs <- c(
  "cancer-dataset-supernatant/HEALTHY/",
  "cancer-dataset-supernatant/LYMPHOMA/",
  "cancer-dataset-supernatant/MYELOMA/"
)

colors <- c("red", "blue", "green")

data <- loadDirectories(dataDirs, colors)

binnedPeaksMatrix <- getBinnedPeaksMatrix(data)

consensusData <- toConsensusSpectraData(data, POP = 0.6)

consensusBinnedPeaksMatrix <- getBinnedPeaksMatrix(consensusData)
