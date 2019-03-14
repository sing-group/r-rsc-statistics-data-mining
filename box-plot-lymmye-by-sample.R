source("load-cancer.R")

# Each spectra belongs to a sample with 5 technical replicates.
# The spectra name is prefixed with the corresponding sample name, and here
# we extract the sample name for each spectra.
sampleNames <- sapply(data$names, function(sample) substr(sample, 1, 2));

# The spectra length is the number of peaks that it holds.
spectraLengths <- sapply(data$spectra, length);

# We also calculate the mean intensity of the peaks in the spectra.
intensityMeans <- sapply(data$spectra, function(spectra) mean(intensity(spectra)));

# Data is grouped in a data frame, so that we can use it in the boxplot function.
dataFrame <- data.frame(
  "spectra" = data$names,
  "samples" = sampleNames,
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
