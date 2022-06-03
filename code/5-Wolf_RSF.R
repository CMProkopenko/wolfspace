# Wolf 2016-2017 RSFs
# Author: S. Zabihi-Seissan, C.M. Prokopenko, E. Vander Wal


### PACKAGES----

library(data.table)
library(lme4)
library(piecewiseSEM)
library(ggplot2)
library(ggeffects)
library(corrplot)
library(ggpubr)
library(AICcmodavg)


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


# MODELS----

Model1<-glmer(Used~M_D30+E_D30+M_H30+E_H30+M_V30+E_V30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)
summary(Model1)
Model2<-glmer(Used~M_D30*E_D30+M_H30+E_H30+M_V30+E_V30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)
summary(Model2)
Model3<-glmer(Used~M_D30+E_D30+M_H30*E_H30+M_V30+E_V30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)
summary(Model3)
Model4<-glmer(Used~M_D30+E_D30+M_H30+E_H30+M_V30*E_V30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)
summary(Model4)
Model5<-glmer(Used~M_D30*M_V30+E_D30*E_V30+M_H30+E_H30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)
summary(Model5)
Model6<-glmer(Used~M_D30+E_D30+M_H30*M_V30+E_H30*E_V30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)
summary(Model6)
Models<-list(Model1,Model2,Model3,Model4,Model5,Model6)
aictab(Models)

# Model 6 wins


# Calculate r-squared and variance inflation factor (change model name for different model outputs)
rsquared(Model6)
vif(Model6)

# Calculate confidence intervals
Wolf_Coefs<-(confint(Model6,level=0.95,method="boot"))

(t1 <- system.time(
  gm1 <- glmer(Used~M_D30+E_D30+M_H30*M_V30+E_H30*E_V30+(1|PackID/WolfID),family=binomial(link="logit"),data=Data)))
##    user  system elapsed 
##   0.188   0.000   0.186

nranpars <- length(getME(gm1,"theta"))
nfixpars <- length(fixef(gm1))

(t2 <- system.time(c1 <- confint(Model6,method="boot", nsim=200,
                                 parm=(nranpars+1):(nranpars+nfixpars),
                                 .progress="txt")))


### Code for generating figures ----

#Predict probability
MooseA<-ggpredict(Model6,terms=c("M_H30 [all]"),ci.lvl=0.95)
ElkA<-ggpredict(Model6,terms=c("E_H30 [all]"),ci.lvl=0.95)
MooseA$Prey<-"Moose_A"
ElkA$Prey<- "Elk_A"
PreyA <- rbind(MooseA,ElkA)
PreyA$group<-NULL
PreyA$conf.high<-NULL
PreyA$conf.low<-NULL
AvailabilityPlot<-ggplot()+
  geom_ribbon(data=MooseA, aes(x,ymin=conf.low,ymax=conf.high),fill="grey35",alpha=0.25)+
  geom_ribbon(data=ElkA, aes(x,ymin=conf.low,ymax=conf.high),fill="grey69",alpha=0.25)+
  geom_line(data=PreyA, aes(x,predicted,linetype=Prey),color="black")+
  ylab("Predicted probability of use")+xlab("Prey habitat selection")+
  annotate("text", x = 1.05, y = 0.15, label = "a",fontface="bold",size=7)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"))+
  scale_linetype_manual(guide = guide_legend(override.aes = list(alpha = 0) ),values = c("Moose_A" = "solid","Elk_A" = "dashed"),labels = c("Elk", "Moose"))+
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white")))
AvailabilityPlot

MooseD<-ggpredict(Model6,terms=c("M_D30 [all]"),ci.lvl=0.95)
ElkD<-ggpredict(Model6,terms=c("E_D30 [all]"),ci.lvl=0.95)
MooseD$Prey<-"Moose_D"
ElkD$Prey<- "Elk_D"
PreyD <- rbind(MooseD,ElkD)
PreyD$group<-NULL
PreyD$conf.high<-NULL
PreyD$conf.low<-NULL
DensityPlot<-ggplot()+
  geom_ribbon(data=MooseD, aes(x,ymin=conf.low,ymax=conf.high),fill="grey35",alpha=0.25)+
  geom_ribbon(data=ElkD, aes(x,ymin=conf.low,ymax=conf.high),fill="grey69",alpha=0.25)+
  geom_line(data=PreyD, aes(x,predicted,linetype=Prey),color="black")+
  ylab("Predicted probability of use")+xlab("Prey density")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  annotate("text", x = 1.05, y = 0.14, label = "b",fontface="bold",size=7)+
  scale_linetype_manual(values = c("Moose_D" = "solid","Elk_D" = "dashed"),labels = c("Elk", "Moose"))
DensityPlot

MooseV<-ggpredict(Model6,terms=c("M_V30 [all]"),ci.lvl=0.95)
ElkV<-ggpredict(Model6,terms=c("E_V30 [all]"),ci.lvl=0.95)
MooseV$Prey<-"Moose_V"
ElkV$Prey<- "Elk_V"
PreyV <- rbind(MooseV,ElkV)
PreyV$group<-NULL
PreyV$conf.high<-NULL
PreyV$conf.low<-NULL
VulnerabilityPlot<-ggplot()+
  geom_ribbon(data=MooseV, aes(x,ymin=conf.low,ymax=conf.high),fill="grey35",alpha=0.25)+
  geom_ribbon(data=ElkV, aes(x,ymin=conf.low,ymax=conf.high),fill="grey69",alpha=0.25)+
  geom_line(data=PreyV, aes(x,predicted,linetype=Prey),color="black")+
  ylab("Predicted probability of use")+xlab("Prey catchability")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white"),
        panel.background = element_blank())+
  annotate("text", x = 1.05, y = 0.60, label = "c",fontface="bold",size=7)+
  scale_linetype_manual(guide = guide_legend(override.aes = list(alpha = 0) ),values = c("Moose_V" = "solid","Elk_V" = "dashed"),labels = c("Elk", "Moose"))
VulnerabilityPlot

# Combine plots

effect_plot<-ggarrange(AvailabilityPlot,DensityPlot,VulnerabilityPlot,nrow=3,ncol=1)
ggsave("effect_plot.jpeg",plot=effect_plot,dpi=450, width = 4.5, height = 7.5, units = "in")

### Interaction plot----

MooseX<-ggpredict(Model6,terms=c("M_H30 [all]","M_V30 [all]"),ci.lvl=NA,)
MooseX$group<-as.factor(MooseX$group)
MooseX<-as.data.table(MooseX)
MooseX<-MooseX[group==0|group==1|group==0.474|group==0.257]
InteractionPlot_moose<-ggplot()+
  geom_line(data=MooseX, aes(x,predicted,linetype=group),color="black")+
  ylab("Predicted probability of use")+xlab("Moose catchability")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_linetype_manual(values = c("0" = "solid","0.257"="dotted","0.474"="dashed","1" = "longdash"),labels = c("0","0.257","0.474", "1"),name="Moose habitat selection")
InteractionPlot_moose
ggsave("interaction_plot.jpeg",plot=InteractionPlot_moose,dpi=450, width = 7.5, height = 4.5, units = "in")

#Full model matrix
CorrData<-Data
colnames(CorrData)[colnames(CorrData)=="M_D30"] <- "Moose density"
colnames(CorrData)[colnames(CorrData)=="M_H30"] <- "Moose habitat"
colnames(CorrData)[colnames(CorrData)=="M_V30"] <- "Moose catchability"
colnames(CorrData)[colnames(CorrData)=="E_D30"] <- "Elk density"
colnames(CorrData)[colnames(CorrData)=="E_H30"] <- "Elk habitat"
colnames(CorrData)[colnames(CorrData)=="E_V30"] <- "Elk catchability"
FullModelcorr<-corrplot(cor(CorrData[,c("Moose density","Moose habitat","Moose catchability","Elk density","Elk habitat","Elk catchability")]),col=c("#F0F0F0", "White", "White"),addCoef.col = 'black',type = 'lower',cl.pos = "n",tl.col = "black")
