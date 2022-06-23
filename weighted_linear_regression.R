# Weighted OLS Linear regression
library(MASS)
library(caret)
library(lmtest)
library(psych)
library(car)

rm(list=ls())
setwd("E:/OneDrive - The University of Colorado Denver/levothyroxine_dose_project/")

# Prescribable LT4 doses
LT4doses=c(50, 75, 88, 100, 112, 125, 137, 150, 175, 200, 224, 250, 275, 300)


train_polished <- read.csv("./manuscript/data/test_train_split/training_split_deid_May_8_2022.csv", stringsAsFactors = F)
test_polished <- read.csv("./manuscript/data/test_train_split/test_split_deid_May_8_2022.csv", stringsAsFactors = F)

# centering height for interaction term
train_polished$Height <- scale(train_polished$Height, center = TRUE, scale = FALSE)
test_polished$Height <- scale(test_polished$Height, center = TRUE, scale = FALSE)

#################################
# EXPLORING BIVARIATE STATISTICS
#################################

bivariate <- train_polished[, c("LT4", "logTSH", "Weight", "Height", "Gender", "Age", "Calcium"), ]

pdf("./manuscript/Results/Weighted_LR/bivariate_analysis.pdf")
pairs.panels(bivariate, col="red", lm=TRUE)
dev.off()

#################################################
# Fitting parsimonious model on training data
################################################

formula <- "LT4 ~ logTSH + Weight + Height:Gender + Gender + Age + Calcium"

# initial fit (unweighted)
fit <- lm(formula, data = train_polished)
summary(fit) 

plot(fit, which = 1, cook.levels = cutoff) # LT4 and predictors are linearly dependent. 
plot(fit, which = 2, cook.levels = cutoff) # Residuals are approximately normally distributed. 
plot(fit, which = 3, cook.levels = cutoff) # Heteroscedasticity is present, which is confirmed by Breusch-Pagan Test
bptest(fit)

# To fix heteroscedasticity using weighted linear regression
#calculating weights
wt <- 1 / lm(abs(fit$residuals) ~ fit$fitted.values)$fitted.values^2

#perform weighted least squares regression
wls_fit <- lm(formula, data = train_polished, weights=wt)

#view summary of model
summary(wls_fit)

# Residual plots
pdf("./manuscript/Results/Weighted_LR/Residual_plots.pdf")
plot(wls_fit, which = 1, cook.levels = cutoff) # linear association is maintained in a weighted model
plot(wls_fit, which = 2, cook.levels = cutoff) # residuals are nearly perfectly distributed
plot(wls_fit, which = 3, cook.levels = cutoff) # there is no heteroscedasticity anymore, which is confirmed by
#perform Breusch-Pagan Test
bptest(wls_fit)

### CHECKING INFLATION
vif(wls_fit) # minor inflation is not suprising considering that weight/height and gender are correlated. 
dev.off()


# testing on independent dataset
test_polished$predLT4 <- predict(wls_fit, newdata = test_polished)
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

#dataframe to collect source data and predictions for CV folds. 
predictions <- train_polished[0 ,]
predictions$predLT4 <- numeric()

set.seed(123)
nfolds = 5

ctrl <- trainControl(method = "cv", number = nfolds)
#define weights to use
wt <- 1 / lm(abs(fit$residuals) ~ fit$fitted.values)$fitted.values^2
wls_model <- train(as.formula(formula), data = train_polished, method = "lm", trControl = ctrl, weights = wt)

print(wls_model)


for (k in 1:nfolds){
  
  train.set <- train_polished[unlist(wls_model$control$index[k]), ]
  test.set <- train_polished[unlist(wls_model$control$indexOut[k]), ]
  
  fit.tmp <- lm(formula, data = train.set)
  wt <- 1 / lm(abs(fit.tmp$residuals) ~ fit.tmp$fitted.values)$fitted.values^2
  
  wls_model.tmp <- lm(formula, data = train.set, weights=wt)
  
  
  print(paste0("Fold: ", k))
  print(summary(wls_model.tmp))
  
  test.set$predLT4 <- predict(wls_model.tmp, newdata = test.set)
  
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


##########################################################
# Combining train and test for production model training
##########################################################

combined <- rbind(predictions, test_polished)

# initial fit (unweighted)
fit <- lm(formula, data = combined)
summary(fit) 


