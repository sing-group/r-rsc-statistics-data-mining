source("load-cancer.R")
source("multiple-sample-visualization-functions.R")

imagesDirectory <- "images/multiple-sample-visualization/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	Intensity Matrix Visualization
#
#	---------------------------------------------------------------------------

library("gplots")

#
# The plotIntensityMatrix function is defined in the "multiple-sample-visualization-functions.R" 
# source file
#
png(paste0(imagesDirectory, "intensity-matrix.png"), width = 1200, height = 800)
plotIntensityMatrix(binnedPeaksMatrix)
dev.off()

#	---------------------------------------------------------------------------
#
#	Pairwise Spectra Comparison: visualization of shared and unshared peaks
#
#	---------------------------------------------------------------------------

a <- binnedPeaksMatrix[1,]
b <- binnedPeaksMatrix[2,]

col.onlyA <- "#7285a5"
col.onlyB <- "#4b560e"
col.both <- "#ff7369"

comparison.parts <- 6
comparison.lwd <- 2
masses <- as.numeric(colnames(binnedPeaksMatrix))
min.mass <- min(as.numeric(colnames(binnedPeaksMatrix)))
max.mass <- max(as.numeric(colnames(binnedPeaksMatrix)))
comparison.range <- (max.mass - min.mass) / comparison.parts

png(paste0(imagesDirectory, "spectra-comparison-multi-", comparison.parts, ".png"), width = 1200, height = 800)
par(mfrow=c(comparison.parts,1), mar=c(2, 2, 2, 2), oma = c(4, 4, 4, 4))

for(part in 1:comparison.parts) {
	peak.start <- min.mass + (part-1) * comparison.range
	peak.end <- min(max.mass, min.mass + part * comparison.range)
	
	indexes <- intersect(which(masses > peak.start), which(masses < peak.end))
	
	plot(NULL, xlim=c(peak.start, peak.end), ylim=c(0,1), xlab="m/z", ylab = NA, yaxt='n')
	
	for(index in 1:length(indexes)) {
		i <- indexes[index]
		if(a[i] > 0 && b[i] == 0) {
			abline(v=as.numeric(colnames(binnedPeaksMatrix)[i]), col=col.onlyA, lwd=comparison.lwd)
		} else if(a[i] == 0 && b[i] > 0) {
			abline(v=as.numeric(colnames(binnedPeaksMatrix)[i]), col=col.onlyB, lwd=comparison.lwd)
		} else if(a[i] > 0 && b[i] > 0) {
			abline(v=as.numeric(colnames(binnedPeaksMatrix)[i]), col=col.both, lwd=comparison.lwd)
		}
	}
}
dev.off()

#	---------------------------------------------------------------------------
#
#	Inverse Spectra Comparison
#
#	---------------------------------------------------------------------------

#
# The compareSpectra function is defined in the "multiple-sample-visualization-functions.R" 
# source file
#

positive <- data$spectra[[1]]
negative <- data$spectra[[2]]

png(paste0(imagesDirectory, "inverse-comparison.png"), width = 1200, height = 800)
compareSpectra(positive, negative)
dev.off()
