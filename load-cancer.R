library("MALDIquant")

# Load all files inside the dataDir directory and returns a list of MassPeaks 
# objects.
#
# Files are expected to be comma-separated files with two columns: mass (m/z)
# and intensity.
loadDirectory <- function(dataDir) {
    spectra <- list.files(dataDir)
    spectraData <- list()
    for (spectrumFile in spectra){
   	 spectrumFile <- paste(dataDir,"/", spectrumFile, sep='')
   	 data <- read.csv(spectrumFile)
   	 spectrum <- createMassPeaks(mass=data[,1], intensity=data[,2], metaData=list(name=spectrumFile))
   	 spectraData <- c(spectraData, spectrum)
    }
    
    spectraData
}

loadSamples <- function(dataDir) {
    samples <- list.dirs(dataDir, full.names=FALSE, recursive=FALSE)
    spectra <- list()
    names <- list()
    sampleNames <- list()
    for (sampleDir in samples){
   	 sampleDirectory <- paste(dataDir,"/", sampleDir, sep='')
   	 sampleSpectra <- loadDirectory(sampleDirectory)
   	 spectra <- c(spectra, sampleSpectra)
   	 names <- c(names, sapply(1:length(sampleSpectra), FUN=function(x) paste(sampleDir,"_R",x,sep='')))
   	 sampleNames <- c(sampleNames, rep(sampleDir, length(sampleSpectra)))
    }
    
    list(names=names, sampleNames=sampleNames, spectra=spectra)
}

loadDirectories <- function(dataDirs, col) {
    if(missing(col)) {
   	 palette <- sample(colors(TRUE))[1:length(dataDirs)]
    } else {
   	 palette <- col[1:length(dataDirs)]
    }
    spectraColors <- list()
    names <- list()
    sampleNames <- list()
    spectra <- list()

    i <- 1
    for (dataDir in dataDirs){
   	 data <- loadSamples(dataDir)
   	 names <- unlist(c(names, data$names))
   	 sampleNames <- unlist(c(sampleNames, data$sampleNames))
   	 spectra <- c(spectra, data$spectra)
   	 spectraColors <- unlist(c(spectraColors, rep(palette[i], length(data$spectra))))
   	 i <- i+1
    }
    
    list(names=names, sampleNames=sampleNames, spectra=spectra, spectraColors=spectraColors)
}

getBinnedPeaksMatrix <- function(data, tolerance=0.002, peakIntensityThreshold=0) {
    binnedPeaks <- binPeaks(data$spectra, tolerance=tolerance)
    binnedPeaksMatrix <- intensityMatrix(binnedPeaks)
    binnedPeaksMatrix[is.na(binnedPeaksMatrix)] <- 0
    binnedPeaksMatrix <- binnedPeaksMatrix[,apply(binnedPeaksMatrix, 2, max) >= peakIntensityThreshold]
    rownames(binnedPeaksMatrix) <- data$names

    binnedPeaksMatrix
}

toConsensusSpectraData <- function(data, tolerance=0.002, POP=0.5) {
	consensusSpectra <- list()
	consensusColors <- vector()
	sampleNames <- unique(data$sampleNames)
	firstReplicateIndex <- 1

	for(i in 1:length(sampleNames)) {
		currentSample <- sampleNames[i]
		currentSampleIndexes <- which(data$sampleNames == currentSample)
		numSamples <- length(currentSampleIndexes)
		sampleSpectra <- data$spectra[currentSampleIndexes]

		binnedPeaks <- binPeaks(sampleSpectra, tolerance=tolerance);
		binnedPeaksMatrix <- intensityMatrix(binnedPeaks);

		consensusMasses <- vector()
		consensusIntensities <- vector()

		for(j in 1:ncol(binnedPeaksMatrix)) {
			consensusPeaks <- which(binnedPeaksMatrix[,j] != 0)
			currentPresence <- length(consensusPeaks) / numSamples
			if(currentPresence >= POP) {
				consensusIntensities <- c(consensusIntensities, mean(binnedPeaksMatrix[consensusPeaks, j]))
				consensusMasses <- c(consensusMasses, colnames(binnedPeaksMatrix)[j])
			}
		}

		currentConsensusSpectrum <- createMassPeaks(
			mass=as.numeric(consensusMasses),
			intensity=as.numeric(consensusIntensities),
			metaData=list(name=currentSample)
		)

		consensusSpectra <- c(consensusSpectra, currentConsensusSpectrum)
		consensusColors <- c(consensusColors, data$spectraColors[firstReplicateIndex])
		firstReplicateIndex <- firstReplicateIndex + numSamples
	}

	list(names=sampleNames, spectra=consensusSpectra, spectraColors=consensusColors)
}

asPresenceMatrix <- function(data) {
	toret <- data
	toret[data > 0] <- 1;
	toret[data == NA] <- 0;
	toret
}

dataDirs <- c(
    "cancer-dataset-supernatant/HEALTHY/",
    "cancer-dataset-supernatant/LYMPHOMA/",
    "cancer-dataset-supernatant/MYELOMA/"
)

colors <- c("red", "blue", "green")

data <- loadDirectories(dataDirs, colors)

binnedPeaksMatrix <- getBinnedPeaksMatrix(data)

consensusData <- toConsensusSpectraData(data, POP=0.6)

consensusBinnedPeaksMatrix <- getBinnedPeaksMatrix(consensusData)
