meanConditionSpectra <-
  function(conditions, peaksMatrix, rowsConditions) {
    avgPeaks <- c()
    
    for (condition in conditions) {
      spectraIndexes <- which(rowsConditions == condition)
      avgPeaks <-
        rbind(avgPeaks, apply(peaksMatrix[spectraIndexes, ], 2, mean))
    }
    
    rownames(avgPeaks) <- conditions
    
    avgPeaks
  }

compareConditions <- function(averagePeaksMatrix, zeroIncrease = 0) {
  avgPeaks <- averagePeaksMatrix
  comparisons <- combn(rownames(avgPeaks), 2)
  result <-
    matrix(nrow = ncol(avgPeaks), ncol = ncol(comparisons))
  rownames(result) <- colnames(avgPeaks)
  colnames(result) <-
    apply(comparisons, 2, function(x)
      paste0(x[1], " / ", x[2]))
  
  for (peak in colnames(avgPeaks)) {
    for (cmpIndex in 1:ncol(comparisons)) {
      a <- comparisons[1, cmpIndex]
      b <- comparisons[2, cmpIndex]
      
      avg.a <- avgPeaks[a, peak]
      if (avg.a == 0) {
        avg.a <- zeroIncrease
      }
      
      avg.b <- avgPeaks[b, peak]
      if (avg.b == 0) {
        avg.b <- zeroIncrease
      }
      
      avg <- avg.a / avg.b
      if (is.infinite(avg) | is.nan(avg)) {
        avg <- NA
      }
      result[peak, cmpIndex] <- avg
    }
  }
  
  result
}
