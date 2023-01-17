library(MASS)
library(caret)
library(lmtest)
library(psych)
library(car)

rm(list=ls())
setwd("C:/Users/nikit/OneDrive - The University of Colorado Denver/levothyroxine_dose_project/manuscript/data/")

# Prescribable LT4 doses
LT4doses=c(50, 75, 88, 100, 112, 125, 137, 150, 175, 200, 224, 250, 275, 300)

train_polished <- read.csv("./training_split.csv", stringsAsFactors = F)
test_polished <- read.csv("./test_split.csv", stringsAsFactors = F)

# centering height for interaction term
train_polished$Heightunscaled <- train_polished$Height
train_polished$Height <- scale(train_polished$Height, center = TRUE, scale = FALSE)
test_polished$Heightunscaled <- test_polished$Height
test_polished$Height <- scale(test_polished$Height, center = TRUE, scale = FALSE)

# For Ordinal Regression Round prescribed doses to nearest prescribable dose
test_polished$rounded_LT4dose <- unlist(apply(test_polished[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))
train_polished$rounded_LT4dose <- unlist(apply(train_polished[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

#################################################
# Fitting parsimonious model on training data
################################################

formula <- "rounded_LT4dose ~ logTSH + Weight + Height + Height:Sex + Sex + Age + Calcium"

#Run Ordinal Regression
Ordinal_data <- train_polished

#Define Levels (Questionable since its all the prescribed doses, not the prescribable doses)
Ordinal_data$rounded_LT4dose <- factor(Ordinal_data$rounded_LT4dose)

mdl <- polr(formula , data = Ordinal_data, Hess = TRUE)
model_sum <- summary(mdl)

#store table
ctable <- coef(summary(mdl))

#Calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

#Combined table with p value
(ctable <- cbind(ctable, "p value" = p))

#testing on independent dataset
test_polished$predLT4 <- predict(mdl, newdata = test_polished)
test_polished$predLT4 <- as.numeric(as.character(test_polished$predLT4))

test_polished$diff <- test_polished$LT4 - test_polished$predLT4
test_polished.norm <- test_polished[test_polished$postopTSH > 0.45 & test_polished$postopTSH <5, ]

# calculating fraction of accurate predictions using +/- 12.5 metric
print(paste0("Test dataset accuracy on all data: ", sum(abs(test_polished$diff) <=12.5)/nrow(test_polished)))
print(paste0("Test dataset accuracy on samples with normal TSH: ", sum(abs(test_polished.norm$diff) <=12.5)/nrow(test_polished.norm)))

# rounding prescribed doses to formulary
test_polished$rounded_LT4dose <- unlist(apply(test_polished[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# rounding predicted doses to formulary 
test_polished$predicted_rounded_LT4dose <- unlist(apply(test_polished[, "predLT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# calculating % correct
sum(test_polished$rounded_LT4dose == test_polished$predicted_rounded_LT4dose)/nrow(test_polished)

# CROSSVALIDATION

#dataframe for CV folds predictions 
predictions <- train_polished[0 ,]
predictions$predLT4 <- numeric()

train_polished$rounded_LT4dose <- factor(train_polished$rounded_LT4dose)

set.seed(123)
nfolds = 5

ctrl <- trainControl(method = "cv", number = nfolds)
cv_model <- train(as.formula(formula), data = train_polished, method = "polr", Hess = TRUE, trControl = ctrl)

print(cv_model)

for (k in 1:nfolds){
  
  train.set <- train_polished[unlist(cv_model$control$index[k]), ]
  test.set <- train_polished[unlist(cv_model$control$indexOut[k]), ]
  
  fit.tmp <- polr(formula, data = train.set, Hess = TRUE)
  
  print(paste0("Fold: ", k))
  print(summary(fit.tmp))
  
  test.set$predLT4 <- predict(fit.tmp, newdata = test.set)
  
  predictions <- rbind(predictions, test.set)
}

# estimating accuracy
predictions$predLT4 <- as.numeric(as.character(predictions$predLT4))
predictions$LT4 <- as.numeric(as.character(predictions$LT4))

predictions$diff <- predictions$LT4 - predictions$predLT4
predictions.norm <- predictions[predictions$postopTSH > 0.45 & predictions$postopTSH <5, ]

# calculating fraction of accurate predictions using +/- 12.5 metric
print(paste0("CV accuracy on all data: ", sum(abs(predictions$diff) <=12.5)/nrow(predictions)))
print(paste0("CV accuracy on samples with normal TSH: ", sum(abs(predictions.norm$diff) <=12.5)/nrow(predictions.norm)))

# rounding prescribed doses to formulary
predictions$rounded_LT4dose <- unlist(apply(predictions[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# rounding predicted doses to formulary 
predictions$predicted_rounded_LT4dose <- unlist(apply(predictions[, "predLT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# calculating % correct
sum(predictions$rounded_LT4dose == predictions$predicted_rounded_LT4dose)/nrow(predictions)


##########################################################
# Combining train and test for production model training
##########################################################

combined <- rbind(predictions, test_polished)
combined$rounded_LT4dose <- factor(combined$rounded_LT4dose)

# initial fit (unweighted)
fit <- polr(formula, data = combined, Hess = TRUE)
summary(fit) 

# CROSSVALIDATION

#dataframe for CV folds predictions 
predictions <- combined[0 ,]
predictions$predLT4 <- numeric()

set.seed(123)
nfolds = 5

ctrl <- trainControl(method = "cv", number = nfolds)
cv_model <- train(as.formula(formula), data = combined, method = "polr", Hess = TRUE, trControl = ctrl)

print(cv_model)

for (k in 1:nfolds){
  
  train.set <- combined[unlist(cv_model$control$index[k]), ]
  test.set <- combined[unlist(cv_model$control$indexOut[k]), ]
  
  fit.tmp <- polr(formula, data = train.set, Hess = TRUE)
  
  print(paste0("Fold: ", k))
  print(summary(fit.tmp))
  
  test.set$predLT4 <- predict(fit.tmp, newdata = test.set)
  
  predictions <- rbind(predictions, test.set)
}

# estimating accuracy
predictions$predLT4 <- as.numeric(as.character(predictions$predLT4))
predictions$LT4 <- as.numeric(as.character(predictions$LT4))

predictions$diff <- predictions$LT4 - predictions$predLT4
predictions.norm <- predictions[predictions$postopTSH > 0.45 & predictions$postopTSH <5, ]

# calculating fraction of accurate predictions using +/- 12.5 metric
print(paste0("CV accuracy on all data: ", sum(abs(predictions$diff) <=12.5)/nrow(predictions)))
print(paste0("CV accuracy on samples with normal TSH: ", sum(abs(predictions.norm$diff) <=12.5)/nrow(predictions.norm)))

# rounding prescribed doses to formulary
predictions$rounded_LT4dose <- unlist(apply(predictions[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# rounding predicted doses to formulary 
predictions$predicted_rounded_LT4dose <- unlist(apply(predictions[, "predLT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# calculating % correct
sum(predictions$rounded_LT4dose == predictions$predicted_rounded_LT4dose)/nrow(predictions)
