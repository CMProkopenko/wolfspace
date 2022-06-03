### Wolf RSF K-fold test----
# Author: S. Zabihi-Seissan, C.M. Prokopenko, E. Vander Wal
# K-fold method from Robert et al. 2017


### Load packages ----
library(data.table)
library(lme4)


### LOAD DATA----

# 2016 and 2017 data with 600m x 600m availability grid
Data<- fread("data/wolf_data.csv")

# Set year as factor
Data$Year<-as.factor(Data$Year)


# STANDARDIZE VARIABLES----

# Reduced to 99.9th percentile
# Elk vulnerability
Data[E_V30>quantile(Data$E_V30,0.999),E_V30:=quantile(Data$E_V30,0.999)]

# Function
Standardize <- function(x){(x-min(x))/(max(x)-min(x))}

#Moose Habitat
Data$M_H30<-Standardize(Data$M_H30)
#Moose Density
Data$M_D30<-Standardize(Data$M_D30)
#Moose Vulnerability
Data$M_V30<-Standardize(Data$M_V30)
#Elk Habitat
Data$E_H30<-Standardize(Data$E_H30)
#Elk Density
Data$E_D30<-Standardize(Data$E_D30)
#Elk Vulnerability
Data$E_V30<-Standardize(Data$E_V30)


### RUN K-FOLD (3 folds)----

# Model formula
Model<- Used~M_D30+E_D30+M_H30+E_H30+M_V30+E_V30+(1|PackID/WolfID) 

# Fix ID column
Data$WolfID<-sub("W", "", Data$WolfID)
Data$WolfID<-as.numeric(Data$WolfID)
Data$WolfID<-as.factor(Data$WolfID)

#K-fold by pack (Spatially independent)
## 3-fold cross-validation with data split by spatially independent individuals. Data from each individual contribute only to one fold and individuals closer than 20 km are never allocated to the same fold. This results in 3 folds with 8 individuals each, 1 fold with 7 individuals, and 1 fold with 12 individuals. Home ranges of individuals assigned to different folds do not overlap

# Select spatially independent individuals (see PART 9 for how these animals have been selected)
Data$rand.vec = ifelse(Data$WolfID == "15" | Data$WolfID == "11"| Data$WolfID == "26"| Data$WolfID == "1"| Data$WolfID == "2"| Data$WolfID == "10",
                       1,ifelse(Data$WolfID == "25"| Data$WolfID == "9"| Data$WolfID == "19"| Data$WolfID == "6"| Data$WolfID == "3", 
                                2, ifelse(Data$WolfID == "24"| Data$WolfID == "22"| Data$WolfID == "12"| Data$WolfID == "4"| Data$WolfID == "5", 3, "ERROR")))
Data$rand.vec<-as.factor(Data$rand.vec)

## visualize sample size split by folds ##
with(Data, tapply(WolfID, rand.vec, unique))

# Fit the model in all folds but one.
#################################  SLOW STEP  ##################################
model_fold_1 <- glmer(Model, Data[Data$rand.vec != 1,], family = binomial)
model_fold_2 <- glmer(Model, Data[Data$rand.vec != 2,], family = binomial)
model_fold_3 <- glmer(Model, Data[Data$rand.vec != 3,], family = binomial)

# Calculation of the contribution to RSF scores by the levels of the categorical predictor (reference category: coniferous forest) and calculation of RSF scores assuming the exponential form. Repeating it for each withhold fold.
# NOTE THE ATTACH COMMAND HERE - DETACH BELOW
attach(Data)

Data$RSFscores[rand.vec ==1] = predict(model_fold_1,newdata=Data[rand.vec == 1],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==2] = predict(model_fold_2,newdata=Data[rand.vec == 2],type = "response",allow.new.levels = TRUE)
Data$RSFscores[rand.vec ==3] = predict(model_fold_3,newdata=Data[rand.vec == 3],type = "response",allow.new.levels = TRUE)

detach(Data)

# Run the k-fold CV evaluation sensu Boyce et al. 2002
dataset <- Data[complete.cases(Data[,"RSFscores"]),]
rho_model <- numeric(3) ## it will store Spearman's coefficients

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

K_fold_scores<-as.data.table(rho_model)

#Final output
mean(K_fold_scores$rho_model)
