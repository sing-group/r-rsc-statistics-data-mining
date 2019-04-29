library(kohonen)

imagesDirectory <- "images/som/"

dir.create(imagesDirectory,
           recursive = TRUE,
           showWarnings = FALSE)

source("load-cancer.R")

#	---------------------------------------------------------------------------
#
#	1. Self-Organized Maps
#
#	---------------------------------------------------------------------------
#
#	1.1 SOM without conditions using the replicas matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------

set.seed(2019)
somresult <- som(
  scale(binnedPeaksMatrix),
  grid = somgrid(xdim = 13, ydim = 3, topo = "hexagonal"),
  alpha = c(0.5, 0.01),
  rlen = 4000
)

plot(somresult, type = "changes")
plot(somresult, type = "counts")
plot(somresult, type = "quality")
plot(somresult, type = "dist.neighbours")

classification <- somresult$unit.classif
names(classification) <- data$spectraNames
classification

#	---------------------------------------------------------------------------
#
#	1.2 SOM with conditions using the replicas matrix (binnedPeaksMatrix)
#
#	---------------------------------------------------------------------------
set.seed(2019)
somresult <- xyf(
  scale(binnedPeaksMatrix),
  factor(data$spectraConditions),
  grid = somgrid(xdim = 13, ydim = 3, topo = "hexagonal"),
  alpha = c(0.5, 0.01),
  rlen = 4000
)

plot(somresult, type = "codes")
plot(somresult, type = "dist.neighbours")
