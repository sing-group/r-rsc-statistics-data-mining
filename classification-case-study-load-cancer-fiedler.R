library("MALDIquant")
library("MALDIquantForeign")
library("MALDIquantExamples")
library("caret")

# Global preprocessing parameters
mq.params.SNR <- 2
mq.params.halfWindowSize <- 20
mq.params.tolerance <- 0.002

# Dataset loading
spectra <- import(getPathFiedler2009()["spectra"], verbose = FALSE)
spectra.info <-
  read.table(getPathFiedler2009()["info"], sep = ",", header = TRUE)

# Dataset filtering to leave only Leipzig samples
spectra <- spectra[spectra.info$location == "leipzig"]
spectra.info <- spectra.info[spectra.info$location == "leipzig", ]

# Common preprocessing
set.seed(2019)

spectra <- trim(spectra)
spectra <- transformIntensity(spectra, method = "sqrt")
spectra <-
  smoothIntensity(spectra, method = "SavitzkyGolay", halfWindowSize = mq.params.halfWindowSize)
spectra <- removeBaseline(spectra, method = "SNIP", iterations = 150)
spectra <- calibrateIntensity(spectra, method = "TIC")

patientConditions <-
  spectra.info[seq(1, nrow(spectra.info), 4), c("patientID", "health")]

# Train-test partition (75% training vs 25% test)
inTrain <-
  createDataPartition(y = patientConditions$health, p = 0.75, list = FALSE)

# Training subset preprocessing
inTrain.patients <- patientConditions[inTrain, "patientID"]

inTrain.spectraIndexes <-
  which(spectra.info$patientID %in% inTrain.patients)
inTrain.spectra <- spectra[inTrain.spectraIndexes]
inTrain.spectra.info <- spectra.info[inTrain.spectraIndexes, ]

inTrain.alignedSpectra <-
  alignSpectra(
    inTrain.spectra,
    SNR = mq.params.SNR,
    halfWindowSize = mq.params.halfWindowSize,
    tolerance = mq.params.tolerance
  )
inTrain.avgSpectra <-
  averageMassSpectra(inTrain.alignedSpectra, inTrain.spectra.info$patientID)
inTrain.avgSpectra.info <-
  inTrain.spectra.info[!duplicated(inTrain.spectra.info$patientID),]

inTrain.peaks <-
  detectPeaks(inTrain.avgSpectra,
              SNR = mq.params.SNR,
              halfWindowSize = mq.params.halfWindowSize)
inTrain.referencePeaks <-
  referencePeaks(inTrain.peaks, tolerance = mq.params.tolerance)

inTrain.peaks.binned <- binPeaks(inTrain.peaks)
inTrain.peaks.binned.filtered <-
  filterPeaks(
    inTrain.peaks.binned,
    minFrequency = c(0.5, 0.5),
    labels = inTrain.avgSpectra.info$health,
    mergeWhitelists = TRUE
  )

inTrain.binnedPeaksMatrix <-
  intensityMatrix(inTrain.peaks.binned.filtered)
inTrain.binnedPeaksMatrix.conditions <-
  inTrain.avgSpectra.info$health
rownames(inTrain.binnedPeaksMatrix) <-
  inTrain.avgSpectra.info$patientID

## The final training set is created in the format required by the caret library
trainingSet <- as.data.frame(inTrain.binnedPeaksMatrix)
trainingSet[is.na(trainingSet)] <- 0
trainingSet$condition <- inTrain.binnedPeaksMatrix.conditions

## Features are renamed to avoid problems with some classification models
featureNames <-
  c(sapply(1:(ncol(trainingSet) - 1), function(x)
    paste0("V", x)), "condition")

colnames(trainingSet) <- featureNames

# Testing subset preprocessing
inTest.patients <- patientConditions[-inTrain, "patientID"]
inTest.spectraIndexes <-
  which(spectra.info$patientID %in% inTest.patients)

inTest.spectra <- spectra[inTest.spectraIndexes]
inTest.spectra.info <- spectra.info[inTest.spectraIndexes, ]

## Peak alignment is done using the reference peaks calculated for the training set
inTest.alignedSpectra <-
  alignSpectra(
    inTest.spectra,
    reference = inTrain.referencePeaks,
    SNR = mq.params.SNR,
    halfWindowSize = mq.params.halfWindowSize,
    tolerance = mq.params.tolerance
  )
inTest.avgSpectra <-
  averageMassSpectra(inTest.alignedSpectra, inTest.spectra.info$patientID)
inTest.avgSpectra.info <-
  inTest.spectra.info[!duplicated(inTest.spectra.info$patientID),]

inTest.peaks <-
  detectPeaks(inTest.avgSpectra, SNR = mq.params.SNR, halfWindowSize = mq.params.halfWindowSize)

## Peaks are binned to the nearest peak of the training set
binPeaksToReference <-
  function (spectraPeaks,
            referenceMasses,
            tolerance = 0.002,
            mapNATo = 0) {
    closestMasses <-
      match.closest(mass(spectraPeaks),
                    referenceMasses,
                    tolerance = tolerance * mass(spectraPeaks))
    # If two peaks are binned to the same reference peak, only the first one is kept
    closestMasses[duplicated(closestMasses)] <- NA
    
    mapToReference <- function(m, mapTo) {
      index <- match(match(m, referenceMasses), closestMasses)
      ifelse(is.na(index), mapNATo, mapTo[index])
    }
    
    intensities <-
      sapply(referenceMasses, mapToReference, intensity(spectraPeaks))
    snrs <-
      sapply(referenceMasses, mapToReference, snr(spectraPeaks))
    
    createMassPeaks(
      mass = referenceMasses,
      intensity = intensities,
      snr = snrs,
      metaData = metaData(spectraPeaks)
    )
  }

inTest.peaks.binned <- lapply(inTest.peaks,
                              function(peaks, referenceMasses)
                                binPeaksToReference(peaks, referenceMasses, tolerance = mq.params.tolerance),
                              attr(inTrain.binnedPeaksMatrix, "mass"))

inTest.binnedPeaksMatrix <- intensityMatrix(inTest.peaks.binned)
inTest.binnedPeaksMatrix[inTest.binnedPeaksMatrix == 0] <- NA
rownames(inTest.binnedPeaksMatrix) <-
  inTest.avgSpectra.info$patientID

inTest.binnedPeaksMatrix.conditions <- inTest.avgSpectra.info$health

## The final testing set is created in the format required by the caret library
testingSet <- as.data.frame(inTest.binnedPeaksMatrix)
testingSet[is.na(testingSet)] <- 0
testingSet$condition <- inTest.binnedPeaksMatrix.conditions

colnames(testingSet) <- featureNames
