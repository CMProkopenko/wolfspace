# Moose kill site RSF for prey (moose) catchability layers
# Author: S. Zabihi-Seissan, C.M. Prokopenko and E. Vander Wal


### PACKAGES----

library(data.table)
library(car)
library(piecewiseSEM)
library(lme4)


### LOAD DATA----

# 2016 and 2017 data with 600m x 600m availability grid
Data<- fread("data/moose_kill_data.csv")

# Set year as factor
Data$Year<-as.factor(Data$Year)

# Set used as integer
Data$Used<-as.integer(Data$Used)


### MODEL----

Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Edge_D:Year+Stream_D:Year+Ruggedness:Year,family=binomial(link="logit"),data=Data)
summary(Model) # GLOBAL


### Model diagnostics----

# Check for collinearity inflation
vif(Model)

# Calculate R-squared
rsquared(Model)


### Calculate confidence intervals for table----

Moose_Coefs<-(cbind(Odds=coef(Model), confint(Model,level=0.95)))
