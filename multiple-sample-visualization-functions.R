#	---------------------------------------------------------------------------
#
#	Intensity Matrix Visualization
#
#	---------------------------------------------------------------------------

plotIntensityMatrix <- function(peaksMatrix) {
  colorPalette <- colorpanel(n = 10, low = "white", high = "black")
  
  
  heatmap.2(
    peaksMatrix,
    # peaks matrix
    main = "",
    # heat map title
    notecol = "black",
    # change font color of cell labels to black
    density.info = "none",
    # turns off density plot inside color legend
    key = FALSE,
    # turns off the color key
    trace = "none",
    # turns off trace lines inside the heat map
    margins = c(12, 9),
    # widens margins around plot
    col = colorPalette,
    # use on color palette defined earlier
    dendrogram = "none",
    # disable dendogram
    Colv = "NA",
    # disable columns clustering
    Rowv = "NA",
    # disable rows clustering,
    lwid = c(0.5, 5),
    # adjust plot margins
    lhei = c(0.5, 5) # adjust plot margins
  )
}

#	---------------------------------------------------------------------------
#
#	Inverse Spectra Comparison
#
#	---------------------------------------------------------------------------

compareSpectra <- function(positive, negative) {
  plot(
    NULL,
    xlim = c(min(mass(positive), mass(negative)), max(mass(positive), mass(negative))),
    ylim = c(-1.1, 1.1),
    xlab = "m/z",
    ylab = "Intensity"
  )
  for (i in 1:length(positive)) {
    segments(mass(positive)[i], 0, y1 = intensity(positive)[i])
  }
  for (i in 1:length(negative)) {
    segments(mass(negative)[i], 0, y1 = -intensity(negative)[i])
  }
  abline(h = 0)
}
