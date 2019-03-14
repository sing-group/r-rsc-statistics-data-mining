source("load-cancer.R")

#	---------------------------------------------------------------------------
#
#	1. Outlier detection
#
#	---------------------------------------------------------------------------

# The spectra length is the number of peaks that it holds.
spectraLengths <- sapply(data$spectra, length);

# We also calculate the mean intensity of the peaks in the spectra.
intensityMeans <- sapply(data$spectra, function(spectra) mean(intensity(spectra)));

# Data is grouped in a data frame, so that we can use it in the boxplot function.
dataFrame <- data.frame(
  "spectra" = data$spectraNames,
  "samples" = data$sampleNames,
  # "colors" = c("red", "red", "blue", "blue", "blue", "blue", "blue", "green", "green", "green", "green", "green"),#data$spectraColors,
  "spectraLengths" = spectraLengths,
  "intensityMeans" = intensityMeans,
  stringsAsFactors = FALSE
)

# This graphical configuration allows including two plots in the same window.
par(mfrow = c(1,2))

# Generation of the boxplot for the spectra lengths.
bpLength = boxplot(
  spectraLengths ~ samples, data = dataFrame,
  main = "Spectra Size", ylab = "Number of peaks", xlab = "Sample"
  # col = dataFrame$colors
)


# Generation of the boxplot for the mean peak intensity.
bpIntensity = boxplot(
  intensityMeans ~ samples, data = dataFrame,
  main = "Intensities", ylab = "Mean Intensity", xlab = "Sample"
  # , col = dataFrame$colors
)

# The boxplot result includes the "out" attribute, a vector with the outlier values.
# Here we extract the spectra with outlier values in the characteristics analyzed.
dataFrame[dataFrame$spectraLengths %in% bpLength$out,]
dataFrame[dataFrame$intensityMeans %in% bpIntensity$out,]

#	---------------------------------------------------------------------------
#
#	2. Spectra comparison
#
#	---------------------------------------------------------------------------

compareSpectra <- function(positive, negative) {
  plot(NULL, xlim=c(min(mass(positive), mass(negative)), max(mass(positive), mass(negative))), ylim=c(-1.1,1.1), xlab="m/z", ylab = "Intensity")
  for(i in 1:length(positive)) {
      segments(mass(positive)[i], 0, y1 = intensity(positive)[i])
  }
  for(i in 1:length(negative)) {
      segments(mass(negative)[i], 0, y1 = -intensity(negative)[i])
  }
  abline(h=0)
}

par(mfrow = c(2,2))
compareSpectra(data$spectra[[7]], data$spectra[[6]])
compareSpectra(data$spectra[[7]], data$spectra[[8]])
compareSpectra(data$spectra[[7]], data$spectra[[9]])
compareSpectra(data$spectra[[7]], data$spectra[[10]])
