library(MASS)
library(caret)
library(lmtest)
library(psych)
library(car)

#rm(list=ls())
#setwd("E:/OneDrive - The University of Colorado Denver/levothyroxine_dose_project/")

# Prescribable LT4 doses
LT4doses=c(50, 75, 88, 100, 112, 125, 137, 150, 175, 200, 224, 250, 275, 300)

train_polished <- read.csv("./manuscript/data/training_split.csv", stringsAsFactors = F)
test_polished <- read.csv("./manuscript/data/test_split.csv", stringsAsFactors = F)

# centering height for interaction term
train_polished$Height <- scale(train_polished$Height, center = TRUE, scale = FALSE)
test_polished$Height <- scale(test_polished$Height, center = TRUE, scale = FALSE)

# converting levothyroxine doses to whole numbers
train_polished$LT4 <- round(train_polished$LT4, digits = 0)
test_polished$LT4 <- round(test_polished$LT4, digits = 0)


################################
# Fitting parsimonious model
################################

formula <- "LT4 ~ logTSH + Weight + Height:Sex + Sex + Age + Calcium"

fit <- glm(as.formula(formula), data = train_polished, family = poisson)
summary(fit)

# calculating dispersion parameter
print(paste0("Dispersion Parameter: ", (sum(residuals(fit, type = "pearson")^2))/fit$df.residual))


# testing on independent dataset
test_polished$predLT4 <- predict(fit, newdata = test_polished, type = "response")
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

predictions <- train_polished[0 ,]
predictions$predLT4 <- numeric()

set.seed(123)
nfolds = 5
folds <- createFolds(train_polished$LT4, k = nfolds, list = TRUE, returnTrain = FALSE)

for (k in 1:nfolds){
  
  train.set <- train_polished[-folds[[k]], ]
  test.set <- train_polished[folds[[k]], ]
  
  fit.tmp <- glm(formula, data = train.set, family = poisson)
  
  print(paste0("Fold: ", k))
  print(summary(fit.tmp))
  
  test.set$predLT4 <- predict(fit.tmp, newdata = test.set, type = "response")
  
  predictions <- rbind(predictions, test.set)
}


# estimating accuracy
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











