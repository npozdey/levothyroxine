# OLS Linear regression

library(caret)
library(randomForest)
library(caTools)

# rm(list=ls())
# setwd("E:/OneDrive - The University of Colorado Denver/levothyroxine_dose_project/")

# Prescribable LT4 doses
LT4doses=c(50, 75, 88, 100, 112, 125, 137, 150, 175, 200, 224, 250, 275, 300)

train_polished <- read.csv("./manuscript/data/training_split.csv", stringsAsFactors = F)
test_polished <- read.csv("./manuscript/data/test_split.csv", stringsAsFactors = F)

# centering height for interaction term
train_polished$Heightunscaled <- train_polished$Height
train_polished$Height <- scale(train_polished$Height, center = TRUE, scale = FALSE)
test_polished$Heightunscaled <- test_polished$Height
test_polished$Height <- scale(test_polished$Height, center = TRUE, scale = FALSE)

# converting levothyroxine doses to whole numbers
train_polished$LT4 <- round(train_polished$LT4, digits = 0)
test_polished$LT4 <- round(test_polished$LT4, digits = 0)



# rounding prescribed doses to formulary
train_polished$rounded_LT4dose <- unlist(apply(train_polished[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))
test_polished$rounded_LT4dose <- unlist(apply(test_polished[, "LT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))

# rounding predicted doses to formulary 
# test_polished$predicted_rounded_LT4dose <- unlist(apply(test_polished[, "predLT4", drop = FALSE], MARGIN = 1,  function(v) LT4doses[which(abs(LT4doses-v)==min(abs(LT4doses-v)))[[1]]]))


train_polished <- transform(train_polished,
                            rounded_LT4dose=factor(rounded_LT4dose, levels = LT4doses),
                            logTSH = as.numeric(logTSH),
                            Weight=as.numeric(Weight),
                            Height=as.numeric(Height),
                            Sex=as.factor(Sex),
                            Age=as.numeric(Age),
                            Calcium=as.factor(Calcium)
                            )
traindata <- train_polished[, c("rounded_LT4dose", "logTSH", "Weight", "Height", "Sex", "Age", "Calcium")]
traindata$rounded_LT4dose <- droplevels(traindata$rounded_LT4dose)
levels(traindata$rounded_LT4dose)

test_polished <- transform(test_polished,
                           rounded_LT4dose=factor(rounded_LT4dose, levels = levels(traindata$rounded_LT4dose)),
                           logTSH = as.numeric(logTSH),
                           Weight=as.numeric(Weight),
                           Height=as.numeric(Height),
                           Sex=as.factor(Sex),
                           Age=as.numeric(Age),
                           Calcium=as.factor(Calcium)
)
testdata <- test_polished[, c("rounded_LT4dose", "logTSH", "Weight", "Height", "Sex", "Age", "Calcium")]
testdata$rounded_LT4dose

formula <- as.formula("rounded_LT4dose ~ logTSH + Weight + Height + Height:Sex + Sex + Age + Calcium")



for (ntree in c(100,300,500,1000,2000,5000, 10000, 20000)){
  
  print(paste0("Ntree: ", ntree))
  
  err_rate <- numeric(10)
  
  for (i in 1:10) {
    data.rf <- randomForest(formula, data=traindata, ntree = ntree)
    conf <- data.rf$confusion[,-ncol(data.rf$confusion)]
    err_rate[i] <- (1 - sum(diag(conf))/sum(conf)) * 100
    
  }
  
  print(paste0("OOB Error Rate: ", mean(err_rate)))
  
    
}

# It appears that 5000 trees is an optimal number 
#
# [1] "Ntree: 100"
# [1] "OOB Error Rate: 58.9107611548556"
# [1] "Ntree: 300"
# [1] "OOB Error Rate: 57.8608923884514"
# [1] "Ntree: 500"
# [1] "OOB Error Rate: 57.7821522309711"
# [1] "Ntree: 1000"
# [1] "OOB Error Rate: 57.5590551181102"
# [1] "Ntree: 2000"
# [1] "OOB Error Rate: 57.3753280839895"
# [1] "Ntree: 5000"
# [1] "OOB Error Rate: 57.0866141732283"
# [1] "Ntree: 10000"
# [1] "OOB Error Rate: 57.1784776902887"
# [1] "Ntree: 20000"
# [1] "OOB Error Rate: 57.2309711286089"






for (mtry in c(2,3,4,5)){
  
  print(paste0("Mtry: ", mtry))
  
  err_rate <- numeric(10)
  
  for (i in 1:10) {
    data.rf <- randomForest(formula, data=traindata, ntree = 5000, mtry = mtry)
    conf <- data.rf$confusion[,-ncol(data.rf$confusion)]
    err_rate[i] <- (1 - sum(diag(conf))/sum(conf)) * 100
    
  }
  
  print(paste0("OOB Error Rate: ", mean(err_rate)))
  
  
}




accuracy <- numeric(10)

for (i in 1:10){

  # prediction on the test split
  data.rf <- randomForest(formula, data=traindata, ntree = 5000, mtry = 2)
  testdata$prediction <- predict(data.rf, newdata=testdata)
  accuracy[i] <- sum(testdata$rounded_LT4dose == testdata$prediction)/nrow(testdata)
}

paste0("Mean accuracy: ", mean(accuracy))





