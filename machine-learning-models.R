library("caret")

imagesDirectory <- "images/machine-learning/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	1. Data preparation: use 70% for train and keep 30% for test
#
#	---------------------------------------------------------------------------

source("load-maldiquant-cancer-fiedler.R")

data <- as.data.frame(binnedPeaksMatrix)
data[is.na(data)] <- 0
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
#	2. Logistic Regression
#
#	---------------------------------------------------------------------------

lr <- train(condition ~ ., data = train.fixedColnames, method = "glm", trControl = trainControl(method = "none"))

# Train data predictions
lr.train.result <- predict(lr, type="raw")
caret::confusionMatrix(data = lr.train.result, train.fixedColnames$condition)

# Test data predictions
lr.test.result <- predict(lr, type="raw", newdata=test.fixedColnames)
caret::confusionMatrix(data = lr.test.result, test.fixedColnames$condition)

#	---------------------------------------------------------------------------
#
#	3. Decision Tree
#
#	---------------------------------------------------------------------------

library("rpart.plot")

tree <- train(condition ~ ., data = train.fixedColnames, method = "rpart1SE", trControl = trainControl(method = "none"))

# Look at the caret available models for other rpart methods: 
#   https://rdrr.io/cran/caret/man/models.html
#
# For instance, the following method also fits a tree by controlling them
# complexity parameter (cp):
#   tree <- train(condition ~ ., data = train.fixedColnames, method = "rpart", 
# trControl = trainControl(method = "none"), tuneGrid = expand.grid(cp=0.003))

# Train data predictions
lr.train.result <- predict(tree, type="raw")
caret::confusionMatrix(data = lr.train.result, train.fixedColnames$condition)

# Test data predictions
lr.test.result <- predict(tree, type="raw", newdata=test.fixedColnames)
caret::confusionMatrix(data = lr.test.result, test.fixedColnames$condition)

# Simple plot: https://stat.ethz.ch/R-manual/R-devel/library/rpart/html/plot.rpart.html
png(paste0(imagesDirectory, "decision-tree-1.png"), width=1200, height=1200)
plot(tree$finalModel, uniform=TRUE, margin=.05)
text(tree$finalModel, use.n = TRUE)
dev.off()

png(paste0(imagesDirectory, "decision-tree-2.png"), width=1200, height=1200)
rpart.plot(tree$finalModel, box.palette="RdBu", shadow.col="gray", nn=TRUE)
dev.off()

#	---------------------------------------------------------------------------
#
#	4. Random Forest
#
#	---------------------------------------------------------------------------

randomforest <- train(condition~., data = train.fixedColnames, method = "rf", trControl = trainControl(method = "none"))

# Example of how to fit a random forest specifying the number randomly selected 
# predictors at each split:
#   randomforest <- train(condition~., data = train.fixedColnames, method = "rf", 
# trControl = trainControl(method = "none"), tuneGrid = expand.grid(mtry=40))

# Train data predictions
lr.train.result <- predict(randomforest, type="raw")
caret::confusionMatrix(data = lr.train.result, train.fixedColnames$condition)

# Test data predictions
lr.test.result <- predict(randomforest, type="raw", newdata=test.fixedColnames)
caret::confusionMatrix(data = lr.test.result, test.fixedColnames$condition)

# Variable importance
library("randomForest")

randomforest$finalModel$importance
rownames(randomforest$finalModel$importance) <- colnames(train)[-ncol(train)]

varImpPlot(randomforest$finalModel, type=2)
dev.print(png, paste0(imagesDirectory, "random-forest-variable-importance.png"), width=600, height=600)

#	---------------------------------------------------------------------------
#
#	5. Artificial Neural Networks
#
#	---------------------------------------------------------------------------

nn <- train(condition~., data = train.fixedColnames, method = "mlp", trControl = trainControl(method = "none"), preProcess = c("center", "scale"))

nn.train.result <- predict(nn, type="raw")
caret::confusionMatrix(data = nn.train.result, train.fixedColnames$condition)

# Test data predictions
nn.test.result <- predict(nn, type="raw", newdata=test.fixedColnames)
caret::confusionMatrix(data = nn.test.result, test.fixedColnames$condition)

#	---------------------------------------------------------------------------
#
#	6. Support Vector Machines
#
#	---------------------------------------------------------------------------

svm.model <- train(condition~., data = train.fixedColnames, method = "svmLinear2", trControl = trainControl(method = "none"))

# Train data predictions
svm.pred.train <- predict(svm.model)
caret::confusionMatrix(data = svm.pred.train, train.fixedColnames$condition)

# Test data predictions
svm.pred.test <- predict(svm.model, newdata=test.fixedColnames)
caret::confusionMatrix(data = svm.pred.test, test.fixedColnames$condition)
