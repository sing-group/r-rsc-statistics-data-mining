source("load-cancer.R")

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