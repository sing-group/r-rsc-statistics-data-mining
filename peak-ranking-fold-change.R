source("load-cancer.R")

imagesDirectory <- "images/important-variables/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

meanConditionSpectra <- function(conditions, peaksMatrix, rowsConditions) {
    avgPeaks <- c()

    for(condition in conditions) {
        spectraIndexes <- which(rowsConditions == condition)
        avgPeaks <- rbind(avgPeaks, apply(peaksMatrix[spectraIndexes,], 2, mean))
    }

    rownames(avgPeaks) <- conditions
    
    avgPeaks
}

avgPeaks <- meanConditionSpectra(consensusData$datasetConditions, consensusBinnedPeaksMatrix, consensusData$spectraConditions)

# log2(avgPeaks["HEALTHY", "2984.26172166667"] / avgPeaks["LYMPHOMA", "2984.26172166667"])
# log2(avgPeaks["HEALTHY", "2984.26172166667"] / avgPeaks["MYELOMA", "2984.26172166667"])

compareConditions <- function(averagePeaksMatrix, zeroIncrease=0) {
    avgPeaks <- averagePeaksMatrix
    comparisons <- combn(rownames(avgPeaks), 2)
    result <- matrix(nrow = ncol(avgPeaks), ncol = ncol(comparisons))
    rownames(result) <- colnames(avgPeaks)
    colnames(result) <- apply(comparisons,2, function(x) paste0(x[1], " / ", x[2]))

    for(peak in colnames(avgPeaks)) {
        for(cmpIndex in 1:ncol(comparisons)) {
            a <- comparisons[1, cmpIndex]
            b <- comparisons[2, cmpIndex]

            avg.a <- avgPeaks[a, peak]
            if(avg.a == 0) {
                avg.a <- zeroIncrease
            }
            
            avg.b <- avgPeaks[b, peak]
            if(avg.b == 0) {
                avg.b <- zeroIncrease
            }

            avg <- avg.a / avg.b
            if(is.infinite(avg) | is.nan(avg)) {
                avg <- NA
            }
            result[peak, cmpIndex] <- avg
        }
    }
    
    result
}

comparison <- compareConditions(avgPeaks)

comparison.log2 <- log2(comparison)

comparison.index <- 1

data <- sort(comparison.log2[,comparison.index])
data.plot <- data[!is.infinite(data) & abs(data) > 1.5]

png(paste0(imagesDirectory, "peak-ranking-fold-changes.png"), width = 1200, height = 800)
barplot(data.plot, horiz=TRUE, names=substr(names(data.plot), 1, 8), las=1, space=1.5, cex.names=0.8, xlim=c(min(data.plot)-1, max(data.plot)+1), main = colnames(comparison)[comparison.index])
dev.off()
