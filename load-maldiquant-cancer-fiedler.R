library("MALDIquant")
library("MALDIquantForeign")
library("MALDIquantExamples")

spectra <- import(getPathFiedler2009()["spectra"], verbose = FALSE)

spectra.info <-
  read.table(getPathFiedler2009()["info"], sep = ",", header = TRUE)

spectra <- trim(spectra)

set.seed(2019)

spectra <- transformIntensity(spectra, method = "sqrt")

spectra <-
  smoothIntensity(spectra, method = "SavitzkyGolay", halfWindowSize = 20)

spectra <- removeBaseline(spectra, method = "SNIP", iterations = 150)

spectra <- calibrateIntensity(spectra, method = "TIC")

spectra <- alignSpectra(spectra)

avgSpectra <-
  averageMassSpectra(spectra, labels = spectra.info$patientID)
avgSpectra.info <-
  spectra.info[!duplicated(spectra.info$patientID),]

peaks <- detectPeaks(avgSpectra, SNR = 2, halfWindowSize = 20)

peaks.binned <- binPeaks(peaks)

peaks.binned.filtered <-
  filterPeaks(
    peaks.binned,
    minFrequency = c(0.5, 0.5),
    labels = avgSpectra.info$health,
    mergeWhitelists = TRUE
  )

# binnedPeaksMatrix <- intensityMatrix(peaks.binned.filtered, avgSpectra)
binnedPeaksMatrix <- intensityMatrix(peaks.binned.filtered)
rownames(binnedPeaksMatrix) <- avgSpectra.info$patientID

binnedPeaksMatrix.conditions <- avgSpectra.info$health
