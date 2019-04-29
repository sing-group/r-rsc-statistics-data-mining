source("classification-case-study-load-cancer-fiedler.R")

imagesDirectory <- "images/machine-learning/"

dir.create(imagesDirectory, recursive=TRUE, showWarnings=FALSE)

# 5x10-CV configuration
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5
)

# Preprocessing configuration. Features are standardized.
pp <- c("center", "scale")

# List of models that will be evaluated
models <- c("rf", "glm", "mlp", "rpart1SE", "svmLinear2")

# Model evaluation
modelFits <- lapply(models, function(model) {
  # Random seed is reseted on each iteration to use the same partitions for all
  # the models
  set.seed(2019)
  train(condition ~ ., data = trainingSet, method = model, preProcess = pp, trControl = ctrl)
})
names(modelFits) <- models

# Resamples allow model fit comparision
results <- resamples(modelFits)

summary(results)

bwplot(results)
dev.print(png, paste0(imagesDirectory, "classification-case-study-bwplot.png"), width=800, height=800)

dotplot(results)
dev.print(png, paste0(imagesDirectory, "classification-case-study-dotplot.png"), width=800, height=800)

# Final prediction of the testing set
# The model used should be selected based on the previous model comparison
cancerPrediction <- predict(modelFits$mlp, newdata=testingSet)
confusionMatrix(data = cancerPrediction, testingSet$condition)
