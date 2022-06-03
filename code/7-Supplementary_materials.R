### Supplementary Materials----
# Author: S. Zabihi-Seissan, C.M. Prokopenko, E. Vander Wal


### Load packages ----
library(data.table)
library(ggplot2)
library(lme4)
library(corrplot)
library(ggpubr)


### Load data ----
Data_kills<- fread("data/supplementary/kill_site_data.csv")
Data_snow<- fread("data/supplementary/snow_depth_data.csv")
Data_elk<- fread("data/supplementary/elk_gps_data.csv")
Data_elk_coef<-fread("data/supplementary/elk_coef_comp.csv")
Data_prey_catchability<-fread("data/supplementary/prey_catchability_comp.csv")


### Kill site bar plots ----
Data_kills<-Data_kills[Pack_ID=="Gunn Lake"|Pack_ID=="Whitewater Lake"|(Pack_ID=="Baldy Lake"&Year==2016)|Pack_ID=="Ranch Creek"|Pack_ID=="Birdtail Valley"|Pack_ID=="Block"]
diet_plot<-ggplot(Data_kills,aes(Prey_spp,fill=as.factor(Year))) +
  geom_bar(position=position_dodge(preserve = 'single'))+
  scale_fill_manual(name = "Year",values = c("grey", "grey15"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ylab("Number of kills")+xlab("Prey species")
diet_plot
ggsave("diet_plot.jpeg",plot=diet_plot,dpi=450, width = 7, height = 5.5, units = "in")

### Snow accumulation plot ----
snow_plot<-ggplot(Data_snow)+
  geom_line(aes(julian_day,snow_depth_cm,colour=as.factor(year)),size=0.8)+
  scale_color_manual(name="Year",values=c("black","grey"))+
  ylab("Snow depth (cm)")+xlab("Julian day")+
  geom_vline(xintercept=32,linetype="dotted")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
snow_plot
ggsave("snow_plot.jpeg",plot=snow_plot,dpi=450, width = 9, height = 5.5, units = "in")

### GPS Elk data RSF (2016) ----
Model<-glmer(Used~ConBog+MarshGrass+Mixedwood+Opendec+BTrail_D+Road_D+Trail_D+Water_D+Ruggedness+Stream_D+(1|elk_id),family=binomial,data=Data_elk)
summary(Model)



Elk_Coefs<-(confint(Model,level=0.95,method="boot"))

(t1 <- system.time(
  gm1 <- glmer(Used~ConBog+MarshGrass+Mixedwood+Opendec+BTrail_D+Road_D+Trail_D+Water_D+Ruggedness+Stream_D+(1|elk_id),family=binomial,data=Data_elk)))
##    user  system elapsed 
##   0.188   0.000   0.186

nranpars <- length(getME(gm1,"theta"))
nfixpars <- length(fixef(gm1))

(t2 <- system.time(c1 <- confint(gm1,method="boot", nsim=200,
                                 parm=(nranpars+1):(nranpars+nfixpars),
                                 .progress="txt")))

# Plot comparing the aerial survey elk RSF coefficients (2016) to those of the 2016 GPS RSF coefficients.
elk_coef_plot<-ggplot(Data_elk_coef,aes(variable,coefficient,shape=model))+
  geom_point(position=position_dodge(width=.5),size=1.5)+
  geom_errorbar(aes(variable,ymin=ci_low,ymax=ci_high,group=model),position=position_dodge(width=.5))+
  scale_fill_manual(name="Model",values=c("grey45","grey"), labels = c("Aerial survey RSF", "GPS collar RSF"))+
  ylab("Coefficient value")+xlab("Variable")+
  geom_hline(yintercept=0,linetype="dashed")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10),
        legend.text=element_text(size=10))
elk_coef_plot
ggsave("elk_coef_plot.jpeg",plot=elk_coef_plot,dpi=450, width = 9, height = 5.5, units = "in")

### Comparing prey catchability layers with prey catchability including prey density ----
FullModelcorr<-corrplot(cor(Data_prey_catchability[,c("Moose catchability 2016","Moose catchability+density 2016","Moose density 2016","Moose catchability 2017","Moose catchability+density 2017","Moose density 2017","Elk catchability 2016","Elk catchability+density 2016","Elk density 2016","Elk catchability 2017","Elk catchability+density 2017","Elk density 2017")]),col=c("#F0F0F0", "White", "White"),addCoef.col = 'black',type = 'lower',cl.pos = "n",tl.col = "black")

### Elk catchability RSF with density ----
Data<- fread("data/elk_kill_data.csv")
Data$Year<-as.factor(Data$Year)
Data$Used<-as.integer(Data$Used)
Model_elk<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Edge_D:Year+Stream_D:Year+Ruggedness:Year+E_Density:Year,family=binomial(link="logit"),data=Data)
summary(Model_elk)


### Moose catchability RSF with density ----
Data<- fread("data/moose_kill_data.csv")
Data$Year<-as.factor(Data$Year)
Data$Used<-as.integer(Data$Used)
Model_moose<-glm(Used~ConBog:Year+MarshGrass:Year+Mixedwood:Year+Opendec:Year+BTrail_D:Year+Road_D:Year+Trail_D:Year+Water_D:Year+Edge_D:Year+Stream_D:Year+Ruggedness:Year+M_Density:Year,family=binomial(link="logit"),data=Data)
summary(Model_moose) 