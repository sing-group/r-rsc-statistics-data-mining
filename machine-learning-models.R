imagesDirectory <- "images/machine-learning/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

#	---------------------------------------------------------------------------
#
#	1. Data preparation: use 70% for train and keep 30% for test
#
#	---------------------------------------------------------------------------

source("load-maldiquant-cancer-fiedler.R")

data <- as.data.frame(binnedPeaksMatrix)
data <- cbind(data, binnedPeaksMatrix.conditions)
colnames(data)[ncol(data)] <- "condition"

set.seed(2019)

trainSamples <- round(0.7*nrow(data))
trainSamplesIndexes <- sample(1:nrow(data), trainSamples)

train <- data[trainSamplesIndexes,]
test <- data[-trainSamplesIndexes,]

#	---------------------------------------------------------------------------
#
#	2. Logistic Regression
#
#	---------------------------------------------------------------------------

train.lr <- train
train.lr$condition <- ifelse(train.lr$condition == "cancer", 1, 0)

test.lr <- test
test.lr$condition <- ifelse(test.lr$condition == "cancer", 1, 0)

lr <- glm(condition~., family=binomial(link='logit'), data=train.lr)

lr

pred.prob.train <- predict(lr, type='response')

#
# Alternative ways of obtaining the predicted probabilities:
#
# pred.odds <- predict(lr); pred.prob.train <- 1/(1+exp(-pred.odds))
#
# pred.odds <- predict(lr); pred.prob.train <- plogis(pred.odds)
#

lr.train.result <- ifelse(pred.prob.train > 0.5, 1, 0)
table(lr.train.result, train.lr$condition)

pred.prob.test <- predict(lr, newdata=test.lr, type='response')

lr.test.result <- ifelse(pred.prob.test > 0.5, 1, 0)
table(lr.test.result, test.lr$condition)

#	---------------------------------------------------------------------------
#
#	3. Decision Tree
#
#	---------------------------------------------------------------------------

library("rpart")
library("rpart.plot")

tree <- rpart(condition~., method = 'class', data=train)
tree

# Train data predictions
tree.pred.train <- predict(tree, type = 'class')
table(tree.pred.train, train$condition)

# Test data predictions
tree.pred.test <- predict(tree, newdata=test, type = 'class')
table(tree.pred.test, test$condition)

# Simple plot: https://stat.ethz.ch/R-manual/R-devel/library/rpart/html/plot.rpart.html
png(paste0(imagesDirectory, "decision-tree-1.png"), width=1200, height=1200)
plot(tree, uniform=TRUE, margin=.05)
text(tree, use.n = TRUE)
dev.off()

png(paste0(imagesDirectory, "decision-tree-2.png"), width=1200, height=1200)
rpart.plot(tree, box.palette="RdBu", shadow.col="gray", nn=TRUE)
dev.off()

#	---------------------------------------------------------------------------
#
#	4. Random Forest
#
#	---------------------------------------------------------------------------

library("randomForest")

train.rf <- train
test.rf <- test

nn.colnames <- c(sapply(1:(ncol(data)-1), function(x) paste0("v", x)), "condition")

colnames(train.rf) <- nn.colnames
colnames(test.rf) <- nn.colnames

randomforest <- randomForest(condition~., data=train.rf)
randomforest

# Train data predictions
randomforest.pred.train <- predict(randomforest)
table(randomforest.pred.train, train.rf$condition)

# Test data predictions
randomforest.pred.test <- predict(randomforest, newdata=test.rf)
table(randomforest.pred.test, test.rf$condition)

# Variable importance
randomforest$importance
rownames(randomforest$importance) <- colnames(train)[-ncol(train)]

varImpPlot(randomforest, type=2)
dev.print(png, paste0(imagesDirectory, "random-forest-variable-importance.png"), width=600, height=600)

#	---------------------------------------------------------------------------
#
#	5. Artificial Neural Networks
#
#	---------------------------------------------------------------------------

library("neuralnet")

train.nn <- train
train.nn$condition <- ifelse(train.nn$condition == "cancer", 1, 0)

test.nn <- test
test.nn$condition <- ifelse(test.nn$condition == "cancer", 1, 0)

# The neuralnet package does not accept column names starting with numbers, so it
# is necessary to rename them to v1, v2, ... vN.

nn.colnames <- c(sapply(1:(ncol(data)-1), function(x) paste0("v", x)), "condition")

colnames(train.nn) <- nn.colnames
colnames(test.nn) <- nn.colnames

# The neuralnet package does not accept formulas like condition~., so it is 
# necessary to write all the variables.

n <- names(train.nn)
f <- as.formula(paste("condition ~", paste(n[!n %in% "condition"], collapse = " + ")))

nn <- neuralnet(f, data=train.nn, hidden=3)

plot(nn)
dev.print(png, paste0(imagesDirectory, "neural-network.png"), width=1200, height=1200)

# Train data predictions
nn.pred.train <- compute(nn, train.nn)
nn.pred.train.prob <- nn.pred.train$net.result
result.train <- ifelse(nn.pred.train.prob > 0.5, 1, 0)
table(result.train, train.nn$condition)

# Test data predictions
nn.pred.test <- compute(nn, test.nn)
nn.pred.test.prob <- nn.pred.test$net.result
result.test <- ifelse(nn.pred.test.prob > 0.5, 1, 0)
table(result.test, test.nn$condition)

#	---------------------------------------------------------------------------
#
#	6. Support Vector Machines
#
#	---------------------------------------------------------------------------

library("e1071")

svm.model <- svm(condition ~ ., kernel = "linear", data = train) 
svm.model

# Train data predictions
svm.pred.train <- predict(svm.model)
table(svm.pred.train, train$condition)

# Test data predictions
svm.pred.test <- predict(svm.model, newdata=test)
table(svm.pred.test, test$condition)
