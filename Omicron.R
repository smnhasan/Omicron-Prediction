library(neuralnet)
require(foreign)
require(MASS)
require(pROC)
require(survey)
require(ResourceSelection)
require(ROCR)
require(car)
require(ggplot2)
require(maptools)
library(nnet)
library(FSA)
library(caret)
require(mapproj)
require(e1071)
library(randomForest)
library(caret)
library(ROCR)

## Creating index variable 

# Read the Data
data = read.csv("E:\\Bio-informatics\\Prediction - Omicron\\All Mutations - copy.csv", header=T)

data$Type_cat_s <- (data$Type_cat - mean(data$Type_cat)) / sd(data$Type_cat)
data$Gene_Cat_s <- (data$Gene_Cat - mean(data$Gene_Cat)) / sd(data$Gene_Cat)
data$Position_s <- (data$Position - mean(data$Position)) / sd(data$Position)
data$AA1_cat_s <- (data$AA1_cat - mean(data$AA1_cat)) / sd(data$AA1_cat)
data$AA2_cat_s <- (data$AA2_cat - mean(data$AA2_cat)) / sd(data$AA2_cat)
data$Mutation.Score_s <- (data$Mutation.Score - mean(data$Mutation.Score)) / sd(data$Mutation.Score)


traindata <- data[which(data$Variant_cat != 10),]
testdata <- data[which(data$Variant_cat == 10),]

ftable(factor(traindata$Type_cat))


#logis
model <- lm(Type_cat_s ~ Gene_Cat_s + Position_s 
             + AA1_cat_s + AA2_cat_s + Mutation.Score_s,data=traindata)

summary(model)
pred = predict(model, newdata=testdata, type = "response")

pred1 <- ifelse(pred>0, 1, 0)   ## if greatr than median value then 1 otherwise 0
pred1 <- factor(pred1)
confusionMatrix(pred1,factor(testdata$Type_cat))   ## Original Yes=106 No=1926

#auc value
library(ROCR)
detach(package:neuralnet)
pr <- prediction(as.numeric(pred1),as.numeric(testdata$Type_cat))
perf <- performance(pr, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pr,"auc"); 
auc <- as.numeric(auc.tmp@y.values)
auc


#auc curve

prob <- predict(model,type="response")
pred <- prediction(as.numeric(prob),as.numeric(traindata$Type_cat_s))
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
auc

#roc curve

plot(perf, main="ROC Curve ", xlab="specificity",  ylab="sensitivity")  
grid()
abline(0,1, col="blue", lty=2)


#lda
model <- lda(Type_cat_s ~ Gene_Cat_s + Position_s 
             + AA1_cat_s + AA2_cat_s + Mutation.Score_s,data=traindata)

summary(model)
pred = predict(model, newdata=testdata, type = "response")
pred$x = as.numeric(unlist(pred$x))
pred1 <- ifelse(pred$x>0, 1, 0)   ## if greatr than median value then 1 otherwise 0
pred1 <- factor(pred1,levels=c(0,1))
confusionMatrix(pred1,factor(testdata$Type_cat))   ## Original Yes=106 No=1926

#auc value
library(ROCR)
detach(package:neuralnet)
pr <- prediction(as.numeric(pred1),as.numeric(testdata$Type_cat))
perf <- performance(pr, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pr,"auc"); 
auc <- as.numeric(auc.tmp@y.values)
auc


#auc curve

prob <- predict(model,type="response")
pred <- prediction(as.numeric(prob$x),as.numeric(traindata$Type_cat_s))
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
auc

#roc curve

plot(perf, main="ROC Curve ", xlab="specificity",  ylab="sensitivity")  
grid()
abline(0,1, col="blue", lty=2)



#ANN
library(neuralnet)
model <- neuralnet(Type_cat_s ~ Gene_Cat_s + Position_s 
                   + AA1_cat_s + AA2_cat_s + Mutation.Score_s,data=traindata)

summary(model)
pred = predict(model, newdata=testdata, type = "response")
pred1 <- ifelse(pred>0, 1, 0)   ## if greatr than median value then 1 otherwise 0
pred1 <- factor(pred1,levels=c(0,1))
confusionMatrix(pred1,factor(testdata$Type_cat))   ## Original Yes=106 No=1926

#auc value
library(ROCR)
detach(package:neuralnet)
pr <- prediction(as.numeric(pred1),as.numeric(testdata$Type_cat))
perf <- performance(pr, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pr,"auc"); 
auc <- as.numeric(auc.tmp@y.values)
auc


#auc curve

prob <- predict(model,newdata=traindata, type="response")
pred <- prediction(as.numeric(prob),as.numeric(traindata$Type_cat_s))
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
auc

#roc curve

plot(perf, main="ROC Curve ", xlab="specificity",  ylab="sensitivity")  
grid()
abline(0,1, col="blue", lty=2)


#SVM

library(e1071)
model <- svm(Type_cat_s ~ Gene_Cat_s + Position_s 
                   + AA1_cat_s + AA2_cat_s + Mutation.Score_s, data=traindata)

summary(model)
pred = predict(model, newdata=testdata, type = "response")
pred1 <- ifelse(pred>0.44, 1, 0)   ## if greatr than median value then 1 otherwise 0
pred1 <- factor(pred1,levels=c(0,1))
confusionMatrix(pred1,factor(testdata$Type_cat))   ## Original Yes=106 No=1926

#auc value
library(ROCR)
detach(package:neuralnet)
pr <- prediction(as.numeric(pred1),as.numeric(testdata$Type_cat))
perf <- performance(pr, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pr,"auc"); auc <- as.numeric(auc.tmp@y.values)
auc


#auc curve

prob <- predict(model, type="response")
pred <- prediction(as.numeric(prob),as.numeric(traindata$Type_cat_s))
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
auc

#roc curve

plot(perf, main="ROC Curve ", xlab="specificity",  ylab="sensitivity")  
grid()
abline(0,1, col="blue", lty=2)

