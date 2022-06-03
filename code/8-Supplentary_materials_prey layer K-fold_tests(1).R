### Prey layer K-fold tests----
# Author: S. Zabihi-Seissan, C.M. Prokopenko, E. Vander Wal
# K-fold method from Robert et al. 2017


### Load packages ----
library(data.table)
library(lme4)
library(splitstackshape)


### Moose habitat RSF K-fold test (5-folds) ----

# Load data
Data<-fread("data/moose_survey_data.csv")
Data[Count==0,Count:=1]
Data<-expandRows(Data,"Count")
Data$Year<-as.factor(Data$Year)

# Randomly bin data
set.seed(5)
Data$rand.vec<-sample(1:5,nrow(Data),replace=TRUE)

# Fit the model in all folds but one
Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Ruggedness:Year+Stream_D:Year,family=binomial,data=Data)

model_fold_1 <- glm(Model, Data[Data$rand.vec != 1,], family = binomial)
model_fold_2 <- glm(Model, Data[Data$rand.vec != 2,], family = binomial)
model_fold_3 <- glm(Model, Data[Data$rand.vec != 3,], family = binomial)
model_fold_4 <- glm(Model, Data[Data$rand.vec != 4,], family = binomial)
model_fold_5 <- glm(Model, Data[Data$rand.vec != 5,], family = binomial)


# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.
# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(Data)

Data$RSFscores[rand.vec ==1] = predict(model_fold_1,newdata=Data[rand.vec == 1],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==2] = predict(model_fold_2,newdata=Data[rand.vec == 2],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==3] = predict(model_fold_3,newdata=Data[rand.vec == 3],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==4] = predict(model_fold_3,newdata=Data[rand.vec == 4],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==5] = predict(model_fold_3,newdata=Data[rand.vec == 5],type = "response",allow.new.levels = TRUE)

detach(Data)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- Data[complete.cases(Data[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients

# Fold 1
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 2
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[2]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[2] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 3
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[3]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[3] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 4
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[4]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[4] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 5
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[5]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[5] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies


K_fold_scores<-as.data.table(rho_model)

#Final output
mean(K_fold_scores$rho_model)


### Elk habitat RSF K-fold test (5-folds) ----

# Load data
Data<-fread("data/elk_survey_data.csv")
Data[Count==0,Count:=1]
Data<-expandRows(Data,"Count")
Data$Year<-as.factor(Data$Year)

# Randomly bin data
set.seed(5)
Data$rand.vec<-sample(1:5,nrow(Data),replace=TRUE)

# Fit the model in all folds but one
Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Ruggedness:Year+Stream_D:Year,family=binomial,data=Data)

model_fold_1 <- glm(Model, Data[Data$rand.vec != 1,], family = binomial)
model_fold_2 <- glm(Model, Data[Data$rand.vec != 2,], family = binomial)
model_fold_3 <- glm(Model, Data[Data$rand.vec != 3,], family = binomial)
model_fold_4 <- glm(Model, Data[Data$rand.vec != 4,], family = binomial)
model_fold_5 <- glm(Model, Data[Data$rand.vec != 5,], family = binomial)


# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.
# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(Data)

Data$RSFscores[rand.vec ==1] = predict(model_fold_1,newdata=Data[rand.vec == 1],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==2] = predict(model_fold_2,newdata=Data[rand.vec == 2],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==3] = predict(model_fold_3,newdata=Data[rand.vec == 3],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==4] = predict(model_fold_3,newdata=Data[rand.vec == 4],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==5] = predict(model_fold_3,newdata=Data[rand.vec == 5],type = "response",allow.new.levels = TRUE)

detach(Data)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- Data[complete.cases(Data[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients

# Fold 1
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 2
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[2]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[2] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 3
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[3]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[3] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 4
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[4]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[4] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 5
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[5]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[5] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies


K_fold_scores<-as.data.table(rho_model)

#Final output
mean(K_fold_scores$rho_model)


### Moose catchability RSF K-fold test (5-folds) ----

# 2016 and 2017 data with 600m x 600m availability grid
Data<- fread("data/moose_kill_data.csv")

# Set year as factor
Data$Year<-as.factor(Data$Year)

# Set used as integer
Data$Used<-as.integer(Data$Used)

# Randomly bin data
set.seed(5)
Data$rand.vec<-sample(1:5,nrow(Data),replace=TRUE)

# Fit the model in all folds but one
Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Edge_D:Year+Stream_D:Year+Ruggedness:Year,family=binomial(link="logit"),data=Data)

model_fold_1 <- glm(Model, Data[Data$rand.vec != 1,], family = binomial)
model_fold_2 <- glm(Model, Data[Data$rand.vec != 2,], family = binomial)
model_fold_3 <- glm(Model, Data[Data$rand.vec != 3,], family = binomial)
model_fold_4 <- glm(Model, Data[Data$rand.vec != 4,], family = binomial)
model_fold_5 <- glm(Model, Data[Data$rand.vec != 5,], family = binomial)


# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.
# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(Data)

Data$RSFscores[rand.vec ==1] = predict(model_fold_1,newdata=Data[rand.vec == 1],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==2] = predict(model_fold_2,newdata=Data[rand.vec == 2],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==3] = predict(model_fold_3,newdata=Data[rand.vec == 3],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==4] = predict(model_fold_3,newdata=Data[rand.vec == 4],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==5] = predict(model_fold_3,newdata=Data[rand.vec == 5],type = "response",allow.new.levels = TRUE)

detach(Data)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- Data[complete.cases(Data[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients

# Fold 1
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 2
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[2]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[2] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 3
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[3]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[3] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 4
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[4]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[4] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 5
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[5]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[5] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies


K_fold_scores<-as.data.table(rho_model)

#Final output
mean(K_fold_scores$rho_model)


### Elk catchability RSF K-fold test (5-folds) ----

# 2016 and 2017 data with 600m x 600m availability grid
Data<- fread("data/elk_kill_data.csv")

# Set year as factor
Data$Year<-as.factor(Data$Year)

# Set used as integer
Data$Used<-as.integer(Data$Used)

# Randomly bin data
set.seed(5)
Data$rand.vec<-sample(1:5,nrow(Data),replace=TRUE)

# Fit the model in all folds but one
Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Edge_D:Year+Stream_D:Year+Ruggedness:Year,family=binomial(link="logit"),data=Data)

model_fold_1 <- glm(Model, Data[Data$rand.vec != 1,], family = binomial)
model_fold_2 <- glm(Model, Data[Data$rand.vec != 2,], family = binomial)
model_fold_3 <- glm(Model, Data[Data$rand.vec != 3,], family = binomial)
model_fold_4 <- glm(Model, Data[Data$rand.vec != 4,], family = binomial)
model_fold_5 <- glm(Model, Data[Data$rand.vec != 5,], family = binomial)


# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.
# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(Data)

Data$RSFscores[rand.vec ==1] = predict(model_fold_1,newdata=Data[rand.vec == 1],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==2] = predict(model_fold_2,newdata=Data[rand.vec == 2],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==3] = predict(model_fold_3,newdata=Data[rand.vec == 3],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==4] = predict(model_fold_3,newdata=Data[rand.vec == 4],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==5] = predict(model_fold_3,newdata=Data[rand.vec == 5],type = "response",allow.new.levels = TRUE)

detach(Data)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- Data[complete.cases(Data[,"RSFscores"]),]
rho_model <- numeric(5) ## it will store Spearman's coefficients

# Fold 1
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[1]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[1] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 2
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[2]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[2] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 3
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[3]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[3] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 4
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[4]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[4] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies

# Fold 5
fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[5]) # run the procedure for the 1st fold - the for-loop below will work on folds 2 and 3
q.pp <- quantile(fold$RSFscores,probs=seq(0,1,.1))  ## computing quantiles of RSF scores
# --------------------------------------------------------
bin <- rep(NA,length(fold$RSFscores)) ## binning RSF scores (10 bins)
for (j in 1:10){
  bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j
}
used <- fold$Used
# --------------------------------------------------------
a <- table(used,bin) ## Occurrence of presence and available data by bin
a <- t(a) #transpose the table
a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores 
a$areaadjusted <- rep(NA,length(10))
sum0 <- sum(a[,1])
sum1 <- sum(a[,2])
a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
a$bins <- seq(1,10,by=1);a

# --------------------------------------------------------
rho_model[5] <- with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate ## Spearman correlation coefficient between RSF bin ranks and area-adjusted frequencies


K_fold_scores<-as.data.table(rho_model)

#Final output
mean(K_fold_scores$rho_model)
