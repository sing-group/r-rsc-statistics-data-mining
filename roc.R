library("pROC")

imagesDirectory <- "images/roc/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

source("load-maldiquant-cancer-fiedler.R")

binnedPeaksMatrix[is.na(binnedPeaksMatrix)] <- 0

#	---------------------------------------------------------------------------
#
#	1. ROC analysis
#
#	---------------------------------------------------------------------------
#
#	1.1 ROC analysis for biomarkers
#
#	---------------------------------------------------------------------------

png(paste0(imagesDirectory, "roc-biomarker.png"), width = 800, height = 800)
# plot the biomarker

rocobj <- plot.roc(binnedPeaksMatrix.conditions, binnedPeaksMatrix[,"1450.33683095267"])
dev.off()

# show area under roc curve
rocobj$auc
# show best threshold
coords(rocobj, "best", ret=c("threshold", "sensitivity", "specificity"))

# compare 2 biomarkers
png(paste0(imagesDirectory, "roc-2biomarkers.png"), width = 800, height = 800)
rocobj <- plot.roc(binnedPeaksMatrix.conditions, binnedPeaksMatrix[,"1450.33683095267"], col="red")
legend("bottomright", col=c("#ff0000", "#0000ff"), legend=c("peak 1450.33683095267","peak 1546.52890935556"), lwd=2)
rocobj2 <- plot.roc(binnedPeaksMatrix.conditions, binnedPeaksMatrix[,"1546.52890935556"], col="blue", add=TRUE)
dev.off()

# show area under roc curve of the second biomarker
rocobj2$auc

#	---------------------------------------------------------------------------
#
#	1.2 ROC analysis for machine learning models
#
#	---------------------------------------------------------------------------

library("caret")

#	---------------------------------------------------------------------------
#
#	1.2.1 Data preparation
#
#	---------------------------------------------------------------------------

data <- as.data.frame(binnedPeaksMatrix)
data <- cbind(data, binnedPeaksMatrix.conditions)
colnames(data)[ncol(data)] <- "condition"

set.seed(2019)

trainSamplesIndexes <- createDataPartition(y = data$condition, p = 0.7, list = FALSE)

train <- data[trainSamplesIndexes,]
test <- data[-trainSamplesIndexes,]

# Replace column names with V1, V2, ..., Vn to use them in some models that do
# not accept variable names starting with numbers

data.colnames <- c(sapply(1:(ncol(data)-1), function(x) paste0("V", x)), "condition")

train.fixedColnames <- train
colnames(train.fixedColnames) <- data.colnames

test.fixedColnames <- test
colnames(test.fixedColnames) <- data.colnames

#	---------------------------------------------------------------------------
#
#	1.2.2 ROC for a Artificial Neural Network
#
#	---------------------------------------------------------------------------

nn <- train(condition~., data = train.fixedColnames, method = "mlp", trControl = trainControl(method = "none"), preProcess = c("center", "scale"))
nn.test.result <- predict(nn, type="prob", newdata=test)

testSamples.cancerProb <- nn.test.result[,1] 
testSamples.class <- binnedPeaksMatrix.conditions[match(rownames(nn.test.result), rownames(binnedPeaksMatrix))]

png(paste0(imagesDirectory, "roc-ann.png"), width = 800, height = 800)
rocobj <- plot.roc(testSamples.class, testSamples.cancerProb)
rocobj$auc
dev.off()

