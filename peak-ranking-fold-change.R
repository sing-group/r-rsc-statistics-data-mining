source("load-cancer.R")
source("peak-ranking-fold-change-functions.R")

imagesDirectory <- "images/important-variables/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

avgPeaks <- meanConditionSpectra(consensusData$datasetConditions, consensusBinnedPeaksMatrix, consensusData$spectraConditions)

# log2(avgPeaks["HEALTHY", "2984.26172166667"] / avgPeaks["LYMPHOMA", "2984.26172166667"])
# log2(avgPeaks["HEALTHY", "2984.26172166667"] / avgPeaks["MYELOMA", "2984.26172166667"])

comparison <- compareConditions(avgPeaks)

comparison.log2 <- log2(comparison)

comparison.index <- 1

data <- sort(comparison.log2[,comparison.index])
data.plot <- data[!is.infinite(data) & abs(data) > 1.5]

png(paste0(imagesDirectory, "peak-ranking-fold-changes.png"), width = 1200, height = 800)
barplot(data.plot, horiz=TRUE, names=substr(names(data.plot), 1, 8), las=1, space=1.5, cex.names=0.8, xlim=c(min(data.plot)-1, max(data.plot)+1), main = colnames(comparison)[comparison.index])
dev.off()
