# Moose RSF for prey habitat layers
# Author: S. Zabihi-Seissan, C.M. Prokopenko, E. Vander Wal


### PACKAGES----
library(data.table)
library(splitstackshape)
library(car)
library(piecewiseSEM)


### LOAD DATA----

# DATA 100 percent cover 600m
Moose_Data<- fread("data/moose_survey_data.csv")
Moose_Data[Count==0,Count:=1]
Moose_Data<-expandRows(Moose_Data,"Count")
Moose_Data$Year<-as.factor(Moose_Data$Year)


### MODELS----

# All
Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Ruggedness:Year+Stream_D:Year,family=binomial,data=Moose_Data)
summary(Model)


### Model diagnostics----

# Check for collinearity inflation
vif(Model)

# Calculate R-squared
rsquared(Model)


### Calculate confidence intervals for table----

Moose_Coefs<-(cbind(Odds=coef(Model), confint(Model,level=0.95)))