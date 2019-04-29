source("data-functions.R")

imagesDirectory <- "images/biomarker/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	1. Biomarker discovery
#
#	---------------------------------------------------------------------------
#
#	1.1 Intensities
#
#	---------------------------------------------------------------------------
#
#	1.1.1 Two conditions
#
#	---------------------------------------------------------------------------

source("load-maldiquant-cancer-fiedler.R")

binnedPeaksMatrix[is.na(binnedPeaksMatrix)] <- 0

# Test one peak (the first one)
wilcox.test(binnedPeaksMatrix[,1] ~ binnedPeaksMatrix.conditions)

# Test all peaks
pvals <- c()
for (i in 1:ncol(binnedPeaksMatrix)) {
  pvals[colnames(binnedPeaksMatrix)[i]] <-
    wilcox.test(binnedPeaksMatrix[,i] ~ binnedPeaksMatrix.conditions)$p.value
}
pvals <- p.adjust(pvals, method="fdr")
pvals[pvals < 0.01]

# Show the boxplot for a peak with p-val < 0.01
png(paste0(imagesDirectory, "biomarker-fiedler-intensities.png"), width = 1200, height = 800)

boxplot(binnedPeaksMatrix[,"1450.33683095267"] ~ binnedPeaksMatrix.conditions, ylab="intensity", xlab="condition", main="Peak 1450.33683095267 by condition")

dev.off()

#	---------------------------------------------------------------------------
#
#	1.1.2 Multiple conditions
#
#	---------------------------------------------------------------------------

source("load-maldiquant-species.R")

binnedPeaksMatrix[is.na(binnedPeaksMatrix)] <- 0
pvals <- c()
for (i in 1:ncol(binnedPeaksMatrix)) {
  pvals[colnames(binnedPeaksMatrix)[i]] <-
    kruskal.test(binnedPeaksMatrix[,i] ~ binnedPeaksMatrix.conditions)$p.value
}
pvals <- p.adjust(pvals, method="fdr")

# Show the boxplot for a peak with p-val < 0.01
png(paste0(imagesDirectory, "biomarker-species-intensities.png"), width = 1200, height = 800)

boxplot(binnedPeaksMatrix[,"1974.3176290975"] ~ binnedPeaksMatrix.conditions, ylab="intensity", xlab="species", main="Peak 1974.3176290975 by species")

dev.off()

#	---------------------------------------------------------------------------
#
#	1.2 Presence/absence
#
#	---------------------------------------------------------------------------
#
#	1.2.1 Two conditions
#
#	---------------------------------------------------------------------------

source("load-maldiquant-cancer-fiedler.R")
source("data-functions.R")

binnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix, "present", "absent")

# Test the first peak
contingencyTable <- table(factor(binnedPeaksMatrix[,1], c("present","absent")), binnedPeaksMatrix.conditions)
contingencyTable
fisher.test(contingencyTable)

pvals <- c()
for (i in 1:ncol(binnedPeaksMatrix)) {
  pvals[colnames(binnedPeaksMatrix)[i]] <-
    fisher.test(table(factor(binnedPeaksMatrix[,i], c("present","absent")), binnedPeaksMatrix.conditions))$p.value
}
pvals <- p.adjust(pvals, method="fdr")

# Show the contingency for a significant peak
table(factor(binnedPeaksMatrix[,"2093.05826713668"], c("present","absent")), binnedPeaksMatrix.conditions)

#	---------------------------------------------------------------------------
#
#	1.2.2 Multiple conditions
#
#	---------------------------------------------------------------------------

source("load-maldiquant-species.R")
source("data-functions.R")

binnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix, "present", "absent")

# Test the first peak
contingencyTable <- table(factor(binnedPeaksMatrix[,1], c("present","absent")), binnedPeaksMatrix.conditions)
contingencyTable
chisq.test(contingencyTable)

pvals <- c()
for (i in 1:ncol(binnedPeaksMatrix)) {
  pvals[colnames(binnedPeaksMatrix)[i]] <-
    chisq.test(table(factor(binnedPeaksMatrix[,i], c("present","absent")), binnedPeaksMatrix.conditions))$p.value
}
pvals <- p.adjust(pvals, method="fdr")

# Show the contingency for a significant peak
table(factor(binnedPeaksMatrix[,"2008.24525043867"], c("present","absent")), binnedPeaksMatrix.conditions)

#	---------------------------------------------------------------------------
#
#	1.3 Visualization
#
#	---------------------------------------------------------------------------
#
#	1.3.1 Filtered heatmaps
#
#	---------------------------------------------------------------------------

source("load-maldiquant-cancer-fiedler.R")
source("data-functions.R")

binnedPeaksMatrix <- asPresenceMatrix(binnedPeaksMatrix)
pvals <- c()
for (i in 1:ncol(binnedPeaksMatrix)) {
  pvals[colnames(binnedPeaksMatrix)[i]] <-
    fisher.test(table(factor(binnedPeaksMatrix[,i], c(1,0)), binnedPeaksMatrix.conditions))$p.value
}
pvals <- p.adjust(pvals, method="fdr")

mapConditionToColor <- function(conditions) {
    colorsVector = ifelse(conditions=="cancer", "#D95DA5", "#82D3A5")
    return(colorsVector)
}

# Presence/absence heatmap with the top-20 most significative peaks
png(paste0(imagesDirectory, "biomarker-cancer-top20-heatmap.png"), width = 1200, height = 800)

topPeaksMatrixWithPvals <- t(binnedPeaksMatrix[,names(pvals[pvals<0.05])])
rownames(topPeaksMatrixWithPvals) <- paste(rownames(topPeaksMatrixWithPvals)," (pval: ", as.character(format(pvals[pvals<0.05])), ")")

heatmap(
    topPeaksMatrixWithPvals, 
    scale="none", 
    ColSideColors=mapConditionToColor(binnedPeaksMatrix.conditions), 
    col=rev(terrain.colors(2))
)

legend(
    "left", 
    inset=.02, 
    title="Peak presence", 
    c("present","absent"), 
    fill=terrain.colors(2), 
    horiz=TRUE, 
    cex=0.8
)

legend(
    "bottomleft", 
    inset=.02, 
    title="Group", 
    c("cancer","control"), 
    fill=c("#D95DA5", "#82D3A5"), 
    horiz=TRUE, 
    cex=0.8
)

dev.off()

#	---------------------------------------------------------------------------
#
#	1.3.2 Volcano plot    
#
#	---------------------------------------------------------------------------

#
# Data and test from subsection 1.1.1 Biomarker Discovery / Intensities / Two conditions
#

source("load-maldiquant-cancer-fiedler.R")

binnedPeaksMatrix[is.na(binnedPeaksMatrix)] <- 0

# Test all peaks
pvals <- c()
for (i in 1:ncol(binnedPeaksMatrix)) {
  pvals[colnames(binnedPeaksMatrix)[i]] <-
    wilcox.test(binnedPeaksMatrix[,i] ~ binnedPeaksMatrix.conditions)$p.value
}
pvals <- p.adjust(pvals, method="fdr")

#
# Load the data.functions.R script in order to compute the fold-changes of the
# Fiedler et al. 2009 dataset
#

source("peak-ranking-fold-change-functions.R")

avgPeaks <- meanConditionSpectra(unique(binnedPeaksMatrix.conditions), binnedPeaksMatrix, binnedPeaksMatrix.conditions)

comparison <- compareConditions(avgPeaks)

comparison.log2 <- log2(comparison)

logFC <- comparison.log2[,1]

# 
# Create the Volcano plot data
#

library("calibrate")

abs_fcs <- cbind(as.integer(abs(logFC))) 
o <- order(abs_fcs, decreasing=TRUE)
abs_fcs_order <- abs_fcs[o,]
limit <- abs_fcs_order[1] + 1

names.logFC.filter <- names(which(abs(logFC)>1))
names.pValue.filter <- names(which(pvals<0.05))
names.both.filter <- intersect(names.pValue.filter, names.logFC.filter)

#
# Make the scatter plot
#

png(paste0(imagesDirectory, "biomarker-cancer-volcano-plot.png"), width = 900, height = 600)

plot(
    logFC, -log10(pvals), 
    col=colors[1], pch=20, 
    xlim=c(-limit, limit), 
    xlab="log2(fold change)", 
    ylab="-log10(adjusted p-value)"
)  

#
# Highlight peaks with different colors for those with an |log2FC| > 1 (red), with an 
# adjusted p-value < 0.05 (orange) and both (green)
#

colors <- c("black", "red", "orange", "green")

points(logFC[names.logFC.filter],-log10(pvals[names.logFC.filter]), col=colors[2], pch=20) 
points(logFC[names.pValue.filter],-log10(pvals[names.pValue.filter]), col=colors[3], pch=20) 
points(logFC[names.both.filter],-log10(pvals[names.both.filter]), col=colors[4], pch=20) 

#
# Label the significant peaks
# 
textxy(
    logFC[names.both.filter],
    -log10(pvals[names.both.filter]),
    labs=names.both.filter, 
    cex=0.9
)

legend(
    "bottomright", 
    xjust=1, 
    yjust=1, 
    legend=c("|log2(fold change)|>1", "adjusted p-value < 0.05", "both"), 
    pch=20, 
    cex=1,
    col=c(colors[2],colors[3],colors[4])
)

dev.off()
