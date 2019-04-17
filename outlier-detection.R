source("load-cancer.R")

#	---------------------------------------------------------------------------
#
#	1. Outlier detection
#
#	---------------------------------------------------------------------------

# The spectra length is the number of peaks that it holds.
spectraLengths <- sapply(data$spectra, length);

# Calculate the mean intensity of the peaks in the spectra.
intensityMeans <- sapply(data$spectra, function(spectra) mean(intensity(spectra)));

# Data is grouped in a data frame, so that it can be used in the boxplot function.
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
# The spectra with outlier values in the characteristics analyzed are extracted here.
dataFrame[dataFrame$spectraLengths %in% bpLength$out,]
dataFrame[dataFrame$intensityMeans %in% bpIntensity$out,]

#	---------------------------------------------------------------------------
#
#	2. Spectra comparison
#
#	---------------------------------------------------------------------------

source("multiple-sample-visualization-functions.R")

par(mfrow = c(2,2))
compareSpectra(data$spectra[[7]], data$spectra[[6]])
compareSpectra(data$spectra[[7]], data$spectra[[8]])
compareSpectra(data$spectra[[7]], data$spectra[[9]])
compareSpectra(data$spectra[[7]], data$spectra[[10]])
