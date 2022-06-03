# Elk RSF for prey habitat layers
# Authors: S. Zabihi-Seissan, C.M. Prokopenko and E. Vander Wal


### PACKAGES----

library(data.table)
library(splitstackshape)
library(car)
library(piecewiseSEM)


### LOAD DATA----

# DATA 100 percent cover 600m
Elk_Data<- fread("data/elk_survey_data.csv")
Elk_Data[Count==0,Count:=1]
Elk_Data<-expandRows(Elk_Data,"Count")
Elk_Data$Year<-as.factor(Elk_Data$Year)


### MODELS----

# All
Model<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Ruggedness:Year+Stream_D:Year,family=binomial,data=Elk_Data)
summary(Model)


### Model diagnostics----

# Check for collinearity inflation
vif(Model)

# Calculate R-squared
rsquared(Model)


### Calculate confidence intervals for table----

Elk_Coefs<-(cbind(Odds=coef(Model), confint(Model,level=0.95)))
