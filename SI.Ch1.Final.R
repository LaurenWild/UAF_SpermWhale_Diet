# This code was developed to show evidence for dietary time series in cetacean skin using stable isotope analysis
# Species include sperm whale, humpback whale, and fin whale
# Author: Lauren Wild, lawild@alaska.edu
# Date published: August 2018
 
############################################################
# Load necessary libraries:
library(ggplot2)
library(Hmisc)
library(psych)
library(devtools)
library(plotrix)  #library for the std.error function; can use mean_se() as well, get same answer
library(MASS)
library(MuMIn) #for dredge
library(lme4)
library(car)
library(nlme)
library(multcomp) #for glht to compare lme parameters
library(AICcmodavg)
library(stats)
library(ggpubr)
library(plyr)

#################################################################################
#### Data for Part 1 - all 3 species - goes until line 727; 
#### Data for Part 2 - just Pm skin (n=28) - starts line 727;
#################################################################################

#One sample, 10 layers:
Pm1<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmIsotopes1_forR.csv',sep=",",header=TRUE)
View(Pm1)
Pm1
str(Pm1) # gives structure of each column (what type of variable)

Mn1<-read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/MnIsotopes1_forR.csv',sep=",",header=TRUE)
View(Mn1) 
CpreciseMn<-0.22
NpreciseMn<-0.13

Bp1<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/BpIsotopes1_forR.csv',sep=",",header=TRUE)
View(Bp1)
Bp2<- Bp1[-c(2,3,4),]
View(Bp2)
Bp2<-Bp1
View(Bp2)
Bp2$d13C2<-NA
Bp2$d13C2<- c(-20.19, -20.22, -20.25, -20.41, -20.28, -20.14, -20.01, -20.17, -20.16, -19.95)

### Exploratory Plots of each species:
######
ggplot(Pm1,aes(d13C,d15N,label=SampleNumber)) + geom_point(size=5) +
  #geom_errorbarh(aes(xmax=Pm1$d13C+0.06,xmin=Pm1$d13C-0.06, height = 0.01)) +
  #geom_errorbar(aes(ymax=Pm1$d15N+0.24,ymin=Pm1$d15N-0.24, width = 0.01))+
  geom_text(hjust=-0.2,vjust = -0.6, size = 6)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-17.4,-16.1,0.2), limits=c(-17.4,-16.1))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(16.5,18.0,0.2), limits=c(16.5,18.0))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

ggplot(Mn1,aes(d13C,d15N,label=SampleNumber)) + geom_point(size=5) +
  #geom_errorbarh(aes(xmax=Pm1$d13C+0.06,xmin=Pm1$d13C-0.06, height = 0.01)) +
  #geom_errorbar(aes(ymax=Pm1$d15N+0.24,ymin=Pm1$d15N-0.24, width = 0.01))+
  geom_text(hjust=-0.2,vjust = -0.6, size = 6)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") 

ggplot(Bp1,aes(d13C,d15N,label=SampleNumber)) + geom_point(size=5)+
  geom_text(hjust=-0.2,vjust = -0.6, size = 6)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

ggplot(Pm1,aes(C.N,d13C,label=SampleNumber, fill=Layer, shape=Vertical, stroke=2)) + geom_point(size=4) +
  geom_text(hjust=-0.2,vjust = -0.6, size = 4)+
  #scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-17.4,-16.1,0.2), limits=c(-17.4,-16.1))+
  #scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(16.5,18.0,0.2), limits=c(16.5,18.0))+
  xlab("C:N")+
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(legend.position='none')+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)
#######

#New 10-samples manuscript figures:
Palette1<- c('black','white', 'grey', 'black')
Shape1<- c(3,21,22,24)

setwd('/Users/laurenwild/Desktop')
tiff(filename="Pm10samp.tif", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
p3<-ggplot(Pm1,aes(d13C,d15N,label=SampleNumber, fill=Layer, shape=Vertical, stroke=1)) + geom_point(size=4) +
  geom_text(hjust=-0.2,vjust = -0.6, size = 4)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-17.4,-16.1,0.2), limits=c(-17.4,-16.1))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(16.5,18.0,0.5), limits=c(16.5,18.0))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(legend.position='none')+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)
dev.off()

setwd('/Users/laurenwild/Desktop')
tiff(filename="Mn10samp.tif", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
p2<- ggplot(Mn1,aes(d13C,d15N,label=SampleNumber, fill=Layer, shape=Vertical, stroke=1)) + geom_point(size=4) +
  geom_text(hjust=-1,vjust = -0.3, size = 4)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-18.0,-16.8,0.2), limits=c(-18.0,-16.8))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(11.5,13.8,0.5), limits=c(11.5,13.8))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.position='none')+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)+
  theme(axis.title.x=element_blank())
dev.off()

setwd('/Users/laurenwild/Desktop')
tiff(filename="Bp10samp.tif", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
p1<-ggplot(Bp1,aes(d13C,d15N,label=SampleNumber, fill=Layer, shape=Vertical, stroke=1)) + geom_point(size=4) +
  geom_text(hjust=-0.8,vjust = -0.8, size = 4)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-20.8,-19.8,0.2), limits=c(-20.8,-19.8))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(11.0,13.0,0.5), limits=c(11.0,13.0))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.position='none')+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)+
  theme(axis.title.x=element_blank())
dev.off()

setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig3_Revised.tif", height = 24, width = 14, units = 'cm', 
     compression = "lzw", res = 300)
ggarrange(p1,p2,p3, labels = c("A", "B", "C"), ncol=1, nrow=3)
dev.off()

ggplot(Bp2,aes(d13C2,C.N,label=SampleNumber, fill=Layer, shape=Vertical, stroke=1)) + geom_point(size=4) +
  geom_text(hjust=-0.8,vjust = -0.8, size = 4)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-20.8,-19.8,0.2), limits=c(-20.8,-19.8))+
  #scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(11.0,13.0,0.5), limits=c(11.0,13.0))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab("C:N ratio")+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(legend.position='none')+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)+
  theme(axis.title.x=element_blank())

##########   ------------------------------  ##############

library(plotrix)  #library for the std.error function; can use mean_se() as well, get same answer

### Summary statistics on each species: 
mean(Pm1$d13C) #-16.788
std.error(Pm1$d13C) #0.124176
sd(Pm1$d13C) #0.39
mean(Pm1$d15N) #17.193
std.error(Pm1$d15N) #0.11201
sd(Pm1$d15N) #0.35
mean(Pm1$C.N) #3.35
std.error(Pm1$C.N) #0.08
sd(Pm1$C.N) #0.26

mean(Mn1$d13C) #-17.472
std.error(Mn1$d13C) #0.08
sd(Mn1$d13C) #0.25
mean(Mn1$d15N) #12.354
std.error(Mn1$d15N) #0.2393
sd(Mn1$d15N) #0.76
mean(Mn1$C.N) #3.28
std.error(Mn1$C.N) #0.03
sd(Mn1$C.N) #0.08

mean(Bp1$d13C) #-20.227
std.error(Bp1$d13C) #0.065
sd(Bp1$d13C) #0.21
mean(Bp1$d15N) #12.081
std.error(Bp1$d15N) #0.11
sd(Bp1$d15N) #0.35
mean(Bp1$C.N) #3.53
std.error(Bp1$C.N) #0.02
sd(Bp1$C.N) #0.07

#subset each layer to calculate mean and std. error of layers
OuterPm<-subset(Pm1[1:3,])
MiddlePm<-subset(Pm1[4:6,])
InnerPm<-subset(Pm1[7:9,])
mean(OuterPm$d13C) #-17.09
mean(OuterPm$d15N) #17.24333
std.error(OuterPm$d13C) #0.140119
mean_se(OuterPm$d13C) #0.14012 as well! Doesn't require plotrix package 
std.error(OuterPm$d15N)  #0.23255
range(OuterPm$d13C)
range(OuterPm$d15N)
sd(OuterPm$d13C) #0.243
sd(OuterPm$d15N) #0.403
range(MiddlePm$d13C)
range(MiddlePm$d15N)
sd(MiddlePm$d13C) #0.268
sd(MiddlePm$d15N) #0.197
range(InnerPm$d13C)
range(InnerPm$d15N)
sd(InnerPm$d13C) #0.093
sd(InnerPm$d15N) #0.235

OuterMn<-subset(Mn1[1:3,])
MiddleMn<-subset(Mn1[4:6,])
InnerMn<-subset(Mn1[7:9,])
mean(OuterMn$d13C) #-17.35
std.error(OuterMn$d13C) #0.25
mean_se(OuterMn$d13C) # 
mean(OuterMn$d15N) #13.41667
std.error(OuterMn$d15N)  #0.069
mean(MiddleMn$d15N) #12.08
std.error(MiddleMn$d15N) #0.059
mean(InnerMn$d15N) #11.72
std.error(InnerMn$d15N) #0.107
range(OuterMn$d13C)
range(OuterMn$d15N)
sd(OuterMn$d13C) #0.425
sd(OuterMn$d15N) #0.119
range(MiddleMn$d13C)
range(MiddleMn$d15N)
sd(MiddleMn$d13C) #0.133
sd(MiddleMn$d15N) #0.102
range(InnerMn$d13C)
range(InnerMn$d15N)
sd(InnerMn$d13C) #0.173
sd(InnerMn$d15N) #0.185

OuterBp<-subset(Bp1[1:3,])
MiddleBp<-subset(Bp1[4:6,])
InnerBp<-subset(Bp1[7:9,])
mean(OuterBp$d13C) # -20.3833
std.error(OuterBp$d13C) # 0.1486
mean_se(OuterBp$d13C) # 
mean(OuterBp$d15N) #12.3133
std.error(OuterBp$d15N)  #0.214
mean(MiddleBp$d15N) #11.813
std.error(MiddleBp$d15N) #0.262
mean(InnerBp$d15N) #12.09
std.error(InnerBp$d15N) #0.0436
range(OuterBp$d13C)
range(OuterBp$d15N)
sd(OuterBp$d13C) #0.257
sd(OuterBp$d15N) #0.37
range(MiddleBp$d13C)
range(MiddleBp$d15N)
sd(MiddleBp$d13C) #0.135
sd(MiddleBp$d15N) #0.454
range(InnerBp$d13C)
range(InnerBp$d15N)
sd(InnerBp$d13C) #0.089
sd(InnerBp$d15N) #0.075

# Avg (+/- SD) range of layers, cores

#First Try on ANOVAs to determine if sperm whale Layer or Column varies significantly by d13C, d15N, or C:N
# Not ideal to use these because they only measure one variable at a time... 
#DO NOT use in final, have moved on... 
#######
fit2<-aov(C.N~Layer, data=Pm1[1:9,])
summary(fit2) #significant p=0.00507
TukeyHSD(fit2) #Check p-values for each layer against each other; O-I most sig, but all are <0.1
fit3<-aov(d15N~Layer, data=Pm1[1:9,])
summary(fit3) #not signficant p=0.101
TukeyHSD(fit3)
fit3.2<-lm(d15N~Layer, data=Pm1[1:9,])
summary(fit3.2) #layer M is sig different from Intercept (Inner); intercept is sig too. 
#overall p value is 0.101 so I say d15N is not sig by layer
fit4<-aov(d13C~Layer, data=Pm1[1:9,])
summary(fit4) #significant p=0.00815
TukeyHSD(fit4) #O-I layer seems to have the most significance for d13C.
fit4.2<-lm(d13C~Layer, data=Pm1[1:9,])
summary(fit4.2) #all layers M & O, plus intercept are sig different; 
#again overal pvalue is also sig at 0.00815

#ANOVAs to determine if sperm whale column varies significantly by d13C, d15N, or C:N
fit8<-aov(d15N~Vertical, data=Pm1[1:9,])
summary(fit8) #not significant p=0.524
fit9<-aov(d13C~Vertical, data=Pm1[1:9,])
summary(fit9) #not significant p=0.99
fit10<-aov(C.N~Vertical, data=Pm1[1:9,])
summary(fit10) #not significant
########

#Two-way ANCOVAs as suggested by Franz, Sep 2017:
# These work well but show C:N being important, and it really shouldn't be in model
#C:N really just used for extraction of lipids, not really structure. 
########
PmN<-lm(d15N~Layer+Vertical+C.N, data=Pm1[1:9,])
anova(PmN) #No significance for any variable
summary(PmN)
AICc(PmN)
PmN2<-lm(d15N~Layer+Vertical, data=Pm1[1:9,])
summary(PmN2)
AICc(PmN2)
PmN3<-lm(d15N~Layer+C.N, data=Pm1[1:9,])
summary(PmN3)
AICc(PmN3)
PmN4<-lm(d15N~Vertical+C.N, data=Pm1[1:9,])
summary(PmN4)
AIC(PmN4)
PmN5<-lm(d15N~Layer, data=Pm1[1:9,])
summary(PmN5)
AICc(PmN5)
PmN6<-lm(d15N~Vertical, data=Pm1[1:9,])
summary(PmN6)
AICc(PmN6)
PmN7<-lm(d15N~C.N, data=Pm1[1:9,])
summary(PmN7)
AICc(PmN7)

PmC<-lm(d13C~Layer+Vertical+C.N+Layer*C.N, data=Pm1[1:9,])
anova(PmC) #No significance for any variable #gives p-values you find in summary of the model
summary(PmC)
AICc(PmC)
PmC2<-lm(d13C~Layer+Vertical, data=Pm1[1:9,])
summary(PmC2)
AICc(PmC2)
PmC3<-lm(d13C~Layer+C.N, data=Pm1[1:9,])
summary(PmC3)
AICc(PmC3)
PmC4<-lm(d13C~Vertical+C.N, data=Pm1[1:9,])
summary(PmC4)
AICc(PmC4)
PmC5<-lm(d13C~Layer, data=Pm1[1:9,])
summary(PmC5)
AICc(PmC5)
PmC6<-lm(d13C~Vertical, data=Pm1[1:9,])
summary(PmC6)
AIC(PmC6)
PmC7<-lm(d13C~C.N, data=Pm1[1:9,])
summary(PmC7)
AIC(PmC7)
anova(PmC7) #gives same p-value as overal summary of this model.. 
PmC8<-lm(d13C~C.N+Layer*C.N, data=Pm1[1:9,])
summary(PmC8)
AIC(PmC8)
anova(PmC8)
PmC9<-lm(d13C~Layer+Layer*C.N, data=Pm1[1:9,])
summary(PmC9)
AIC(PmC9)
anova(PmC9)
############

#Third try Two-Way ANOVAs as suggested by Geraldine after discussion of C:N
#Determined C:N does not belong in model, 
# USE THESE for final! 
Pm.d15N1<-lm(d15N~Layer+Vertical, data=Pm1[1:9,])
options(na.action = "na.fail")
anova(Pm.d15N1) #no significance
summary(Pm.d15N1)
AICc(Pm.d15N1) #48.9
dredge(Pm.d15N1)
Pm.d15N2<-lm(d15N~Layer, data=Pm1[1:9,])
anova(Pm.d15N2) #None significant
summary(Pm.d15N2)
AICc(Pm.d15N2) #17.75
Pm.d15N3<-lm(d15N~Vertical, data=Pm1[1:9,])
anova(Pm.d15N3) #None significant
summary(Pm.d15N3)
AICc(Pm.d15N3) #22.68
dredge(Pm.d15N1) #Best model is that of intercept only
Pm.d15N0<-lm(d15N~1, data=Pm1[1:9,])
summary(Pm.d15N0)
anova(Pm.d15N0)

Pm.d13C1<-lm(d13C~Layer+Vertical, data=Pm1[1:9,])
anova(Pm.d13C1) #Layer significant
summary(Pm.d13C1)
AICc(Pm.d13C1) #48.13
dredge(Pm.d13C1)
Pm.d13C2<-lm(d13C~Layer, data=Pm1[1:9,])
anova(Pm.d13C2) #Layer very significant
summary(Pm.d13C2)
AICc(Pm.d13C2) #12.29
summary(glht(Pm.d13C2, linfct=mcp(Layer="Tukey")))
Pm.d13C3<-lm(d13C~Vertical, data=Pm1[1:9,])
anova(Pm.d13C3) # No significance
summary(Pm.d13C3)
AICc(Pm.d13C3) #26.68
dredge(Pm.d13C1) #Layer is best model. 

#ANOVAs to determine if humpback whale layer or column varies significantly by d13C, d15N, or C:N
# These were first try, and not the right way to do it individually... 
#######
fit5<-aov(C.N~Layer, data=Mn1[1:9,])
summary(fit5) #almost significant at 5% level p=0.054
TukeyHSD(fit5) #O-I is most significant, p=0.05 ... 
fit6<-aov(d15N~Layer, data=Mn1[1:9,])
summary(fit6) #significant p=1.36e-05
TukeyHSD(fit6) #O-I most sig, also O-M... M-I slightly sig p=0.04
fit6.2<-lm(d15N~Layer, data=Mn1[1:9,])
summary(fit6.2) #significant, same p-value as above
fit7<-aov(d13C~Layer, data=Mn1[1:9,])
summary(fit7) #not significant p =0.498
fit7.2<-lm(d13C~Layer, data=Mn1[1:9,])
summary(fit7.2) #not significant same p-value as above

#ANOVAs to determine if humpback whale column varies significantly by d13C, d15N, or C:N
fit11<-aov(d15N~Vertical, data=Mn1[1:9,])
summary(fit11) #Not significant, p=0.978
fit12<-aov(d13C~Vertical, data=Mn1[1:9,])
summary(fit12)  #Not significant, p=0.466
fit13<-aov(C.N~Vertical, data=Mn1[1:9,])
summary(fit13) # Not significant, p=0.795
###########

#Two-way ANCOVAs as suggested by Franz, Sep 2017
##These were better but shouldn't include C:N, as discussed with Geraldine 28 Sep 2017
##########
anova(lm(d15N~Layer+Vertical+C.N, data=Mn1[1:9,]))
anova(lm(d13C~Layer+Vertical+C.N, data=Mn1[1:9,]))
anova(lm(d15N~Layer, data=Mn1[1:9,]))
anova(lm(d13C~Layer+Vertical+C.N, data=Mn1[1:9,]))
MnN<-lm(d15N~Layer+Vertical+C.N, data=Mn1[1:9,])
anova(MnN) #Layer is significant
summary(MnN)
summary.aov(MnN)
AICc(MnN)
MnN2<-lm(d15N~Layer+Vertical, data=Mn1[1:9,])
summary(MnN2)
AICc(MnN2)
MnN3<-lm(d15N~Layer+C.N, data=Mn1[1:9,])
summary(MnN3)
AICc(MnN3)
MnN4<-lm(d15N~Vertical+C.N, data=Mn1[1:9,])
summary(MnN4)
AICc(MnN4)
MnN5<-lm(d15N~Layer, data=Mn1[1:9,])
summary(MnN5)
AICc(MnN5)
MnN6<-lm(d15N~Vertical, data=Mn1[1:9,])
summary(MnN6)
AICc(MnN6)
MnN7<-lm(d15N~C.N, data=Mn1[1:9,])
summary(MnN7)
AICc(MnN7)
summary(step(MnN))

MnC<-lm(d13C~Layer+Vertical+C.N, data=Mn1[1:9,])
anova(MnC) #Layer is significant
summary(MnC)
summary.aov(MnC)
AICc(MnC)
MnC2<-lm(d13C~Layer+Vertical, data=Mn1[1:9,])
summary(MnC2)
AICc(MnC2)
MnN3<-lm(d15N~Layer+C.N, data=Mn1[1:9,])
summary(MnN3)
AICc(MnN3)
MnC4<-lm(d13C~Vertical+C.N, data=Mn1[1:9,])
summary(MnC4)
AICc(MnC4)
MnC5<-lm(d13C~Layer, data=Mn1[1:9,])
summary(MnC5)
AICc(MnC5)
MnC6<-lm(d13C~Vertical, data=Mn1[1:9,])
summary(MnC6)
AICc(MnC6)
MnC7<-lm(d13C~C.N, data=Mn1[1:9,])
summary(MnC7)
AICc(MnC7)
##################

Mn.d15N1<-lm(d15N~Layer+Vertical, data=Mn1[1:9,])
anova(Mn.d15N1) #layer significant
summary(Mn.d15N1)
AICc(Mn.d15N1) #37.099
dredge(Mn.d15N1)
Mn.d15N2<-lm(d15N~Layer, data=Mn1[1:9,])
anova(Mn.d15N2) #Layer even more significant
summary(Mn.d15N2)
AICc(Mn.d15N2) #4.497
Mn.d15N3<-lm(d15N~Vertical, data=Mn1[1:9,])
anova(Mn.d15N3) #no significance
summary(Mn.d15N3)
AICc(Mn.d15N3) #38.05
dredge(Mn.d15N1)

Mn.d13C1<-lm(d13C~Layer+Vertical, data=Mn1[1:9,])
options(na.action = "na.fail")
anova(Mn.d13C1) #None significant
summary(Mn.d13C1)
AICc(Mn.d13C1) #49.72
dredge(Mn.d13C1)
Mn.d13C2<-lm(d13C~Layer, data=Mn1[1:9,])
anova(Mn.d13C2) #No significance
summary(Mn.d13C2)
AICc(Mn.d13C2) #16.72
Mn.d13C3<-lm(d13C~Vertical, data=Mn1[1:9,])
anova(Mn.d13C3) #No significance
summary(Mn.d13C3)
AICc(Mn.d13C3) #16.52
dredge(Mn.d13C1)


#ANOVAs to determine if Fin whale layer varies sig by d13C, d15N, or C:N
#### Not using these, first version... 
########
fit14<-aov(d13C~Layer, data=Bp1[1:9,])
summary(fit14) #Not significant, p=0.244
fit15<-aov(d15N~Layer, data=Bp1[1:9,])
summary(fit15)  #Not significant, p=0.274
fit16<-aov(C.N~Layer, data=Bp1[1:9,])
summary(fit16)  #Not significant, p=0.125

#ANOVAs to determine if Fin whale column vaires sig by d13C, d15N, or C:N
fit17<-aov(d13C~Vertical, data=Bp1)
summary(fit17)  #Not significant, p=0.395
fit18<-aov(d15N~Vertical, data=Bp1)
summary(fit18)   #Not significant, p=0.247
fit19<-aov(C.N~Vertical, data=Bp1)
summary(fit19)   #Not significant, p=0.72
#########

#Final two-way ANCOVAs as suggested by Franz, Sep 2017
# These not great because don't need C:N, 
########
BpN<-lm(d15N~Layer+Vertical+C.N, data=Bp1[1:9,])
anova(BpN) #No significance for any variable
summary(BpN)
AICc(BpN)
BpN2<-lm(d15N~Layer+Vertical, data=Bp1[1:9,])
summary(BpN2)
AICc(BpN2)
BpN3<-lm(d15N~Layer+C.N, data=Bp1[1:9,])
summary(BpN3)
AICc(BpN3)
BpN4<-lm(d15N~Vertical+C.N, data=Bp1[1:9,])
summary(BpN4)
AICc(BpN4)
BpN5<-lm(d15N~Layer, data=Bp1[1:9,])
summary(BpN5)
AICc(BpN5)
BpN6<-lm(d15N~Vertical, data=Bp1[1:9,])
summary(BpN6)
AICc(BpN6)
BpN7<-lm(d15N~C.N, data=Bp1[1:9,])
summary(BpN7)
AICc(BpN7)

BpC<-lm(d13C~Layer+Vertical+C.N, data=Bp1[1:9,])
anova(BpC) #No significance for any variable
summary(BpC)
AICc(BpC)
BpC2<-lm(d13C~Layer+Vertical, data=Bp1[1:9,])
summary(BpC2)
AICc(BpC2)
BpC3<-lm(d13C~Layer+C.N, data=Bp1[1:9,])
summary(BpC3)
AICc(BpC3)
BpC4<-lm(d13C~Vertical+C.N, data=Bp1[1:9,])
summary(BpC4)
AICc(BpC4)
BpC5<-lm(d13C~Layer, data=Bp1[1:9,])
summary(BpC5)
AICc(BpC5)
BpC6<-lm(d13C~Vertical, data=Bp1[1:9,])
summary(BpC6)
AICc(BpC6)
BpC7<-lm(d13C~C.N, data=Bp1[1:9,])
summary(BpC7)
AICc(BpC7)
#########

#### Two-way ANOVA's advised by geraldine to cut C:N on 28 Sep 2017
Bp.d15N1<-lm(d15N~Layer+Vertical, data=Bp1[1:9,])
anova(Bp.d15N1) #core almost significant, but not
summary(Bp.d15N1)
AICc(Bp.d15N1) #44.92
dredge(Bp.d15N1)
Bp.d15N2<-lm(d15N~Layer, data=Bp1[1:9,])
anova(Bp.d15N2) #No significance
summary(Bp.d15N2)
AICc(Bp.d15N2) #20.524
Bp.d15N3<-lm(d15N~Vertical, data=Bp1[1:9,])
anova(Bp.d15N3) #No significance 
summary(Bp.d15N3)
AICc(Bp.d15N3) #18.683

Bp.d13C1<-lm(d13C~Layer+Vertical, data=Bp1[1:9,])
anova(Bp.d13C1) #None significant
summary(Bp.d13C1)
AICc(Bp.d13C1) #41.35
dredge(Bp.d13C1)
Bp.d13C2<-lm(d13C~Layer, data=Bp1[1:9,])
anova(Bp.d13C2) #No significance
summary(Bp.d13C2)
AICc(Bp.d13C2) #8.58
Bp.d13C3<-lm(d13C~Vertical, data=Bp1[1:9,])
anova(Bp.d13C3) #No significance
summary(Bp.d13C3)
AICc(Bp.d13C3) #10.929

### ----------------------------------

#Concentration of Carbon & Nitrogen in Pm b/c C:N was significantly different - WHY?
plot(Pm1$Layer, Pm1$X.N, xlab="Layers", ylab="Nitrogen Concentration") 
fitXN<-aov(X.N~Layer, data=Pm1)
summary(fitXN) #significant p=0.00561
TukeyHSD(fitXN) #mostly just the O-I difference that is sig, also the M-I slightly

plot(Pm1$Layer, Pm1$X.C, xlab="Layers", ylab="Carbon Concentration")
fitXC<-aov(X.C~Layer, data=Pm1)
summary(fitXC) #Not significant p =0.82

plot(Mn1$Layer, Mn1$X.N, xlab="Layers", ylab="Nitrogen Concentration")
fitXN.Mn<-aov(X.N~Layer, data=Mn1[1:9,])
summary(fitXN.Mn) #Not significant p=0.341
TukeyHSD(fitXN.Mn)
plot(Mn1$Layer, Mn1$X.C, xlab="Layers", ylab="Carbon Concentration")
fitXC.Mn<-aov(X.C~Layer, data=Mn1[1:9,])
summary(fitXC.Mn) #Not significant p=0.708
TukeyHSD(fitXC.Mn)

#### PLOTS ########
plot(Pm1$Layer, Pm1$C.N, xlab="Layers", ylab="C:N ratio")
plot(Pm1$Layer, Pm1$d13C, xlab="Layers", ylab="d13C value")
plot(Pm1$C.N, Pm1$d13C, xlab="C:N", ylab="d13C value")

plot(Mn1$Layer, Mn1$C.N, xlab="Layers", ylab="C:N ratio")
plot(Mn1$Layer, Mn1$d15N, xlab="Layers", ylab="d15N value")
plot(Mn1$C.N, Mn1$d15N, xlab="C:N", ylab="d15N value")

ggplot(Pm1, aes(C.N, d13C)) +geom_point() +
  xlab("C:N")+
  ylab(expression(paste(delta^13, "C \u2030",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

ggplot(Pm1, aes(C.N, d15N)) +geom_point() +
  xlab("C:N")+
  ylab(expression(paste(delta^15, "N \u2030",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

ggplot(Pm1,aes(C.N,d13C,label=SampleNumber, fill=Layer, shape=Layer, stroke=2)) + geom_point(size=4) +
  geom_text(hjust=-0.8,vjust = -0.6, size = 4)+
  #scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-17.4,-16.1,0.2), limits=c(-17.4,-16.1))+
  #scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(16.5,18.0,0.2), limits=c(16.5,18.0))+
  xlab("C:N")+
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(legend.position='none')+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)

ggplot(Pm1,aes(C.N,d15N,label=SampleNumber, fill=Layer, shape=Layer, stroke=2)) + geom_point(size=4) +
  geom_text(hjust=-0.8,vjust = -0.6, size = 4)+
  #scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-17.4,-16.1,0.2), limits=c(-17.4,-16.1))+
  #scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(16.5,18.0,0.2), limits=c(16.5,18.0))+
  xlab("C:N")+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(legend.position='none')+
  scale_fill_manual(values=Palette1)+
  scale_shape_manual(values=Shape1)

#To Recap: d13C was sig dif btwn Pm layers, as was C:N, but XC wasn't sig diff, just XN. 
#d15N was sig dif btwn Mn layers, slight sig in C:N, neither XC nor XN sig diff. 
#Nothing sig diff btwn Bp layers. 

#Combine Datasets for full biplot:
require(plyr)
Mn.Pm.Bp <- rbind.fill(Pm1, Mn1, Bp1)
View(Mn.Pm.Bp)

#use doBy to set up data sheet with averages and std.error of each species
library(doBy)
myfun1<-function(x) {c(m=mean(x) , v=var(x))}  
myfun2<-function(x) {c(m=mean(x), sd=sd(x))}

Mn.Pm.Bp2 <-summaryBy(d15N+d13C~Species, data=Mn.Pm.Bp, FUN=myfun2)
View(Mn.Pm.Bp2)

Mn.Pm.Bp2$Species2 <- c("Sperm Whale", "Humpback Whale", "Fin Whale")
View(Mn.Pm.Bp2)

setwd('/Users/laurenwild/Desktop')
tiff(filename="Mn.Pm.Bp.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggplot(Mn.Pm.Bp2,aes(d13C.m,d15N.m, label=Species2, shape=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=Mn.Pm.Bp2$d13C.m+Mn.Pm.Bp2$d13C.sd,xmin=Mn.Pm.Bp2$d13C.m-Mn.Pm.Bp2$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=Mn.Pm.Bp2$d15N.m+Mn.Pm.Bp2$d15N.sd,ymin=Mn.Pm.Bp2$d15N.m-Mn.Pm.Bp2$d15N.sd, width = 0.01))+
  geom_text(hjust=-0.06,vjust = -0.7, size = 5)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.position="none")
dev.off()

#Plotting ellipses #install.packages("ellipse")#library(ellipse)
setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig2_Revision.tiff", height = 12, width = 16, units = 'cm', 
     compression = "lzw", res = 300)
ggplot(Mn.Pm.Bp,aes(d13C,d15N, shape=Species)) + geom_point(size=4) +
  #geom_text(hjust=-0.8,vjust = -0.8, size = 4)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)",sep="")), breaks=seq(-22.0,-15,1), limits=c(-22.0,-15))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)",sep="")), breaks=seq(10,19,1), limits=c(10,19))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  #stat_ellipse(type="norm", level=0.95)+
  theme(legend.position=c(0.24,0.8))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  scale_shape_manual(values=c(16,17,15), labels=c("Sperm Whale", "Humpback Whale", "Fin Whale"))
dev.off()

##############################################################################
###  Layer Data - Sperm whales - All 3 Layers Only ###
###  Only layer samples that have full 3 layers + full sample (n=28) ###
###  THIS DATA USED FOR 10-SAMPLES PAPER, CH 1;
##############################################################################
Pm3<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmIsotopes_Layers1forR.csv',sep=",",header=TRUE)
View(Pm3)

#New Dataframe exclude the "full" sample, just compare layers#
Pm3.2<-Pm3[Pm3$Layer!='Full',]
View(Pm3.2)

#Make doy term
library(date)
class(Pm3.2$Date) #It's a factor variable
Pm3.2$Date2<- as.POSIXlt(Pm3.2$Date, format='%m/%d/%y')
class(Pm3.2$Date2) #Now class is POSIXlt POSIXt
View(Pm3.2) #Should have dates in Date2 Column! 
doy <- strftime(Pm3.2$Date2, format = "%j") # new vector of dates in julian format
doy<-as.numeric(doy) #make julian date numeric
Pm3.2$doy<-doy #add julian date to data frame
View(Pm3.2) # Check everything looks good! 

#EXPLORATORY PLOTS:
##############
#plot showing each layer's variation from the d15N mean of that sample. 
ggplot(Pm3.2,aes(d15N.Avg,d15N)) + 
  geom_point(aes(shape=Layer), size=3)+
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(method='lm', formula=y~x, se=FALSE, color="black", size=0.5, linetype="dashed")+
  xlab(expression(paste(delta^15, "N (\u2030) ", "Mean")))+
  ylab(expression(paste(delta^15, "N (\u2030) ", "per Layer"))) +
  labs(Shape="Layer")+
  theme_bw()
ggplot(Pm3.2,aes(d13C.Avg,d13C)) + 
  geom_point(aes(shape=Layer), size=3)+
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(method='lm', formula=y~x, se=FALSE, color="black", size=0.5, linetype="dashed")+
  xlab(expression(paste(delta^13, "C (\u2030) ", "Mean")))+
  ylab(expression(paste(delta^13, "C (\u2030) ", "per Layer"))) +
  labs(Shape="Layer")+
  theme_bw()

#plot showing each layer's variation from the d13C mean of that sample
#This one has the smoothed line as the predicted d13C values... don't use! 
ggplot() + 
  geom_point(aes(d13C.Avg, d13C, shape=Layer, size=3),data=Pm3.2)+
  geom_smooth(aes(d13C.Avg, Predict), method='lm', se=FALSE, data=Pm3.2[Pm3.2$Layer=="Inner",])+
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(aes(d13C.Avg, d13C), method='lm', formula=y~x, se=FALSE, color="black", size=0.5, linetype="dashed", data=Pm3.2)+
  geom_smooth(aes(d13C.Avg, d13C), method='lm', se=FALSE, color='red', data=Pm3.2[Pm3.2$Layer=="Inner",])+
  xlab(expression(paste(delta^13, "C (\u2030) ", "Mean")))+
  ylab(expression(paste(delta^13, "C (\u2030) ", "per Layer")))+
  labs(Shape="Layer")+
  theme_bw()

ggplot() + 
  geom_point(aes(d13C.Avg, d13C, shape=Layer, size=3),data=Pm3.2)+
  geom_smooth(aes(d13C.Avg, d13C), method='lm', se=FALSE, color='blue', data=Pm3.2[Pm3.2$Layer=="Inner",])+
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(aes(d13C.Avg, d13C), method='lm', formula=y~x, se=FALSE, color="black", size=0.5, linetype="dashed", data=Pm3.2)+
  geom_smooth(aes(d13C.Avg, d13C), method='lm', se=FALSE, color='red', data=Pm3.2[Pm3.2$Layer!="Inner",])+
  xlab(expression(paste(delta^13, "C (\u2030) ", "Mean")))+
  ylab(expression(paste(delta^13, "C (\u2030) ", "per Layer")))+
  labs(Shape="Layer")+
  theme(legend.position='none')
#theme_bw()

#Plot with d13C for each date sampled by month ... 
ggplot(aes(Date,d13C), data=Pm3.2)+
  geom_point() + 
  xlab("Sampling Date")+ 
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45))

### Plots to check out how C:N could be related to isotope... Geraldine suggested these as well
## Show that C:N isn't correlated to delta values. Completely random. So no relationship. 
## Therefore don't need C:N in model, it's not supposed to predict the delta values. 
ggplot(Pm3.2, aes(C.N, d13C)) +geom_point() +
  xlab("C:N")+
  ylab(expression(paste(delta^13, "C \u2030",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

ggplot(Pm3.2, aes(C.N, d15N)) +geom_point() +
  xlab("C:N")+
  ylab(expression(paste(delta^15, "N \u2030",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

ggplot(Pm3.2, aes(d13C, d15N)) +geom_point() +
  xlab(expression(paste(delta^13, "C \u2030",sep="")))+
  ylab(expression(paste(delta^15, "N \u2030",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

###-----------------------------------------------------------------------####
#####  Mike sugested plots showing deviation from mean on a horizontal plane: #####

Pm3.2["d15N.dev"] <- NA
Pm3.2$d15N.dev<- Pm3.2$d15N-Pm3.2$d15N.Avg

ggplot(Pm3.2,aes(as.factor(Sample.Number),d15N.dev)) + 
  geom_point(aes(shape=Layer), size=3)+
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(method='lm', formula=y~x, se=FALSE, color="black", size=0.5, linetype="dashed")+
  xlab("Sample")+
  ylab(expression(paste(delta^15, "N (\u2030) ", "Deviation from Mean")))+
  labs(Shape="Layer")+
  theme(axis.text.x = element_text(angle=45))

Pm3.2["d13C.dev"] <- NA
Pm3.2$d13C.dev<- Pm3.2$d13C-Pm3.2$d13C.Avg

ggplot(Pm3.2,aes(as.factor(Sample.Number),d13C.dev)) + 
  geom_point(aes(shape=Layer), size=3)+
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(method='lm', formula=y~x, se=FALSE, color="black", size=0.5, linetype="dashed")+
  xlab("Sample")+
  ylab(expression(paste(delta^15, "N (\u2030) ", "Deviation from Mean")))+
  labs(Shape="Layer")+
  theme(axis.text.x = element_text(angle=45))

### Look at range of inner to outer layer: How variable is an individual?
#########
#Range of inner to outer layer d13C values, for histograms
#Check variability of each whale's d13C values... max-min
dfI<-Pm3.2[Pm3.2$Layer=="Inner",]
dfM<-Pm3.2[Pm3.2$Layer=="Middle",]
dfO<-Pm3.2[Pm3.2$Layer=="Outer",]
dfI$d13C-dfO$d13C
hist(dfI$d13C-dfO$d13C, xlab="Inner - Outer Layer", ylab = "Frequency", main="", breaks=8)
hist(dfI$d13C-dfM$d13C, xlab="Inner - Middle Layer", ylab = "Frequency", main="")
hist(dfM$d13C-dfO$d13C, xlab="Middle - Outer Layer", ylab = "Frequency", main="")
dfI$IminusO<-dfI$d13C-dfO$d13C
ggplot() +
  geom_point(aes(doy,IminusO), data=dfI)+
  xlab("Day of Year")+
  ylab("Range of Inner Minus Outer Layer")+
  theme(legend.position='none')
#no real trend here with range of inner minus outer layers not getting smaller over the summer
dfI$IminusM<-dfI$d13C-dfM$d13C
ggplot() +
  geom_point(aes(doy,IminusM), data=dfI)+
  xlab("Day of Year")+
  ylab("Range of Inner Minus Middle Layer")+
  theme(legend.position='none')
#no real trend here with range of inner minus middle layers not getting smaller over the summer
dfI$MminusO<-dfM$d13C-dfO$d13C
ggplot() +
  geom_point(aes(doy,MminusO), data=dfI)+
  xlab("Day of Year")+
  ylab("Range of Middle Minus Outer Layer")+
  theme(legend.position='none')
#no real trend here with range of middle minus outer layers not getting smaller over the summer

# Try to subtract lowest from highest for each Sample.Number
dfN<-aggregate(Pm3.2$d15N, list(Pm3.2$Sample.Number), max)
dfN$M<- aggregate(Pm3.2$d15N, list(Pm3.2$Sample.Number), min)
dfN$D<- dfN$x - dfN$M
View(dfN) #D column gives d15N differences for each sample ... 
#min = 0.08 permil max = 2 permil 
dfN2<- c(0.37,0.79,0.08,0.46,0.52,0.37,0.48,0.54,0.40,0.25,0.63,0.23,0.27,0.14,0.76,0.20,0.39,0.41,0.62,0.38,0.33,0.14,0.12,0.18,1.05,0.32,0.41,0.25)
mean(dfN2) #0.4
std.error(dfN2) #0.04
sd(dfN2) #0.22

dfC<-aggregate(Pm3.2$d13C, list(Pm3.2$Sample.Number), max)
dfC$M<- aggregate(Pm3.2$d13C, list(Pm3.2$Sample.Number), min)
dfC$D<- dfC$x - dfC$M
View(dfC) #D column gives d13C differences for each sample ... 
#min = 0.1 permil #max = 1.44 permil 
dfC2<- c(0.5,0.29,0.28,0.51,0.57,0.92,0.58,0.83,1.01,0.05,0.79,1.21,0.57,0.64,1.44,0.95,0.4,0.57,0.30,0.74,0.86,0.27,0.51,0.74,1.03,0.76,0.5,0.37)
mean(dfC2) #0.65
std.error(dfC2) #0.06
sd(dfC2) #0.31

#These might not actually be a good test ... not really saying what I want here...
t.test(dfC$D, mu=0) #p=2.09e-08; t= -6.55; df=55; saying the range of d13C values within an ind.'s layers is sig diff from 0
t.test(dfN$D, mu=0) #p=2.09e-08; t= -6.55; df=55; saying range of d15N values of layers w/in an ind. is sig diff from 0

######## -----------------------------------------

#Full model: Update from Franz Sep 2017: 
#Decided not to use model with C:N because it's not a relevant variable. 
#####
fit2.1<- lme(d13C ~ Layer + C.N, data=Pm3.2, random =~1|Sample.Number, method = 'ML')
summary(fit2.1) #C:N is significantly different from intercept, which is 'inner' layer
summary(aov(fit2.1)) #here we have significance with "layer" only, not "C:N"

### ------------------------------------------
### These models are tested versions comparing Layer to isotope ratios:
#####
fit2.2<- lme(d13C ~ Layer, data=Pm3.2, random =~1|Sample.Number, method = 'ML')
summary(fit2.2)
anova(fit2.2)
dredge(fit2.2)
plot(fit2.2)
plot(fit2.2$fitted, fit2.2$res)
qqnorm(fit2.2$res) #looks pretty good! 
library(RLRsim)
exactRLRT(fit2.2) #Restricted Likelihood Ratio Test - need to change "method" in model to be REML
test<-anova(fit2.2) #run an anova of the lme, in order to do a Tukey post hoc comparison
test2<-aov(fit2.2)
library(multcomp)
summary(glht(fit2.2, linfct=mcp(Layer="Tukey"))) #shows difference btwn specific layers
summary(glht(test2, linfct=mcp(Layer="Tukey"))) #shows difference btwn specific layers

#Try LME with layer & month as predictors...
fit2.22<- lme(d13C ~ Layer*doy + doy^2 + Year, data=Pm3.2, random =~1|Sample.Number, method = 'ML')
summary(fit2.22)
anova(fit2.22)
dredge(fit2.22) # Best model has Layer & Month, AICc=75.5
plot(fit2.22)
plot(fit2.22$fitted, fit2.22$res)
qqnorm(fit2.22$res) #looks pretty good! 
test<-anova(fit2.22) #run an anova of the lme, in order to do a Tukey post hoc comparison
test2<-aov(fit2.22)
library(multcomp)
summary(glht(fit2.22, linfct=mcp(Layer="Tukey"))) #shows difference btwn specific layers
summary(glht(test2, linfct=mcp(Layer="Tukey")))
summary(glht(test2, linfct=mcp(Month="Tukey")))

# LME with just layer as  predictor, but with 2 random effects: sample & date
fit2.24<- lme(d13C ~ Layer, data=Pm3.2, random = list(~1|Sample.Number, ~1|Date), method = 'ML')
summary(fit2.24)
anova(fit2.24)
dredge(fit2.24) #Layer sig, AICc 89.9


#LME with d15N and Layer only 
fit2.3<- lme(d15N ~ Layer, data=Pm3.2, random =~1|Sample.Number, method = 'ML')
summary(fit2.3) #Middle and outer layers sig diff from intercept, which is 'inner' layer
anova(fit2.3) #Layer is significant @ 5% level (p=0.03)
dredge(fit2.3)
plot(fit2.3$fitted,fit2.3$res)
qqnorm(fit2.3$res) #not great ... normality might be an issue... 
library(RLRsim)
exactRLRT(fit2.3) #Restricted Likelihood Ratio Test - need to change "method" in model to be REML
test3<-anova(fit2.3) #run an anova of the lme, in order to do a Tukey post hoc comparison
test4<-aov(fit2.3)
library(multcomp)
summary(glht(test4, linfct=mcp(Layer="Tukey"))) #shows difference btwn specific layers

####################

#Working on LMEs Franz suggested at Oct. 18th meeting. 

###### FINAL MODELS FOR PUBLICATION: 

#FULL MODEL, WITH QUADRATIC, INTERACTION, & ALL TERMS: 
fitOct1<-lme(d13C ~ poly(doy,2)*Layer + poly(doy,2) + Layer, data=Pm3.2, random = list( ~ 1 |Sample.Number, ~1|Year), method="ML")
summary(fitOct1)
anova(fitOct1) #Layer most sig; but also doy^2, and almost year 
dredge(fitOct1) #Best model is Layer, doy^2, and Year interestingly.

#Reduced model from DREDGE:
fitOct2<-lme(d13C ~  Layer + poly(doy,2), data=Pm3.2, random = list( ~ 1 |Sample.Number, ~1|Year), method="ML")
summary(fitOct2)
anova(fitOct2) 
dredge(fitOct2)

#Without Year as a random effect:
fitOct2.5<-lme(d13C ~  Layer + poly(doy,2), data=Pm3.2, random =  ~ 1 |Sample.Number, method="ML")
summary(fitOct2.5)

# Test if random effect of year improves the model: 
anova(fitOct2.5,fitOct2) #shows difference between model w/ and w/out year as rand.eff. not significant!
# Doesn't improve model, leave year out as a random effect: 
# Final best model is fitOct2.5
dredge(fitOct2.5)
r.squaredGLMM(fitOct2.5) #marginal R2 = 0.30; conditional R2 = 0.86
test6<-aov(fitOct2.5)
library(multcomp)
summary(glht(test6, linfct=mcp(Layer="Tukey"))) #shows difference btwn specific layers
plot(fitOct2.5$fitted,fitOct2.5$res)
hist(fitOct2.5$res)
hist(residuals(fitOct2.5))
qqnorm(residuals(fitOct2.5)) #Pretty good
qqline(residuals(fitOct2.5))

#PREDICTIONS: Feed a new dataframe one individual, with one year, and give it 20 or so DOYs, to predict d13C values
LMEpredict1<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/LMEpredict1_Oct.csv',sep=",",header=TRUE)
View(LMEpredict1)
LMEpredict1$Predict<- predict(fitOct2.5, LMEpredict1)
View(LMEpredict1)
plot(LMEpredict1$doy, LMEpredict1$Predict) 
#shows as doy increases, predicted d13C values get less negative... across all layers
#Use this one, with sample number 158134
LMEpredict2<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/LMEpredict2_Oct.csv',sep=",",header=TRUE)
LMEpredict2$Predict<- predict(fitOct2.5, LMEpredict2)
View(LMEpredict2)
plot(LMEpredict2$doy, LMEpredict2$Predict) 
#Try different individual, different year: 
LMEpredict3<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/LMEpredict3_Oct.csv',sep=",",header=TRUE)
LMEpredict3$Predict<- predict(fitOct2.5, LMEpredict3)
View(LMEpredict3)
plot(LMEpredict3$doy, LMEpredict3$Predict) #Get the same pattern, different d13C scale depending on idiv. 


setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig10_d13C_doy_predicted_Final.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggplot(aes(doy,d13C, shape=Layer), data=Pm3.2)+
  geom_point() + 
  geom_line(aes(doy,Predict, linetype=Layer), data=LMEpredict2)+ 
  xlab("Julian Day")+ 
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw()
dev.off()

##### ---------- FULL MODEL d15N ---------------- ######## 
fitOct3<-lme(d15N ~ poly(doy,2)*Layer + poly(doy,2), data=Pm3.2, random = list( ~ 1 |Sample.Number, ~1|Year), method="ML")
summary(fitOct3)
anova(fitOct3) # layer and doy^2 signifcant
dredge(fitOct3) #best model is with doy and doy^2, and Layer;
fitOct3.5<-lme(d15N ~ poly(doy,2)*Layer + poly(doy,2), data=Pm3.2, random = ~1 | Sample.Number, method="ML")
anova(fitOct3.5,fitOct3) #Shows year not needed as a random effect. 

#Reduced model: 
fitOct4<-lme(d15N ~ Layer + poly(doy,2), data=Pm3.2, random = ~ 1 |Sample.Number, method="ML")
summary(fitOct4)
anova(fitOct4) #Layer sig.; doy^2 slightly signifcant (p=0.01)
r.squaredGLMM(fitOct4) #marginal R-sqrd = 026; conditional R-sqrd = 0.92
dredge(fitOct4)
plot(fitOct4$fitted,fitOct4$res)
plot(fitOct4$res)
qqnorm(fitOct4$res) #not great ... normality might be an issue... 
qqline(fitOct4$res) #not great ... normality might be an issue...
test5<-aov(fitOct4)
library(multcomp)
summary(glht(test5, linfct=mcp(Layer="Tukey"))) #shows difference btwn specific layers

ggplot() +
  geom_point(aes(doy, d15N, shape=Layer, size=3), data=Pm3.2) +
  scale_shape_manual(values=c(1,17,4))+
  geom_smooth(aes(doy, d15N), method='lm', se=FALSE, color="blue", data=Pm3.2[Pm3.2$Layer=="Inner",])+
  geom_smooth(aes(doy, d15N), method='lm', se=FALSE, color="red", data=Pm3.2[Pm3.2$Layer=="Middle",])+
  geom_smooth(aes(doy, d15N), method='lm', se=FALSE, color="black", data=Pm3.2[Pm3.2$Layer=="Outer",])+
  xlab("Day of Year")+
  ylab(expression(paste(delta^15, "N (\u2030) ")))+
  labs(Shape="Layer")+
  theme(legend.position='none')

#PREDICTIONS: 
LMEpredict3$Predict2<- predict(fitOct4, LMEpredict3)
View(LMEpredict2)
LMEpredict1$Predict2<- predict(fitOct4, LMEpredict1)
LMEpredict2$Predict2<- predict(fitOct4, LMEpredict2)

#Figure 9, for Publication: 
setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig9_d15N_doy_predicted_Final.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggplot(aes(doy,d15N, shape=Layer), data=Pm3.2)+
  geom_point() + 
  geom_line(aes(doy,Predict2, linetype=Layer), data=LMEpredict3)+ 
  xlab("Julian Day")+ 
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw()
dev.off()

###COMBINE FIGS 9 & 10 FOR PUBLICATION INTO PANEL: 
MonthId <- c("May", "Jun", "Jul", "Aug", "Sep")
MiddleDay <- c(136, 166, 197, 228, 258)
p9<-ggplot(aes(doy,d15N, shape=Layer), data=Pm3.2)+
  geom_point(size=3) + 
  geom_line(aes(doy,Predict2, linetype=Layer), data=LMEpredict3)+
  xlab("")+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_x_continuous(sec.axis = sec_axis(~ . ,  breaks = MiddleDay, labels =MonthId, name = "Month"))+
  theme(axis.text.x=element_blank())+
  theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())

p10<-ggplot(aes(doy,d13C, shape=Layer), data=Pm3.2)+
  geom_point(size=3) + 
  geom_line(aes(doy,Predict, linetype=Layer), data=LMEpredict2)+ 
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  xlab("Julian Day")+ 
  #scale_x_continuous(sec.axis = sec_axis(~ . , labels =MonthId, breaks = MiddleDay))+
  theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="none")

library(egg)
library(gridExtra)
library(reshape)

setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig6_Revised.tiff", height = 18, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggarrange(p9,p10, labels = c("A", "B"), ncol=1, nrow=2)
dev.off()

multiplot(p9,p10, cols=1)

#Iliana requests a few plots: 
#Single biplot with C on x axis and N on y axis, label the layers with different colors or symbols
#Then connect values from 3 layers with a line for each whale
Pm3.2$Group<-c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,
               5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7)
setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig4_Revised.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 100)
ggplot(aes(d13C, d15N), data=Pm3.2)+
  geom_point(aes(shape=Layer), size=3, data=Pm3.2)+ 
  scale_shape_manual(values=c(19,24,0))+
  #geom_line(aes(group=as.factor(Sample.Number)), data=Pm3.2)+
  geom_path(aes(group=as.factor(Sample.Number)), data=Pm3.2)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_discrete(guide="none")+
  theme_bw(base_size = 18, base_family = "Helvetica")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(legend.position=c(0.12,0.82))+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=10))
#facet_wrap(~ Group, ncol=2 )
dev.off()

#Now normalize all C data to inner layer values (like Mike's suggeted plots, but by Inner layer not mean)
Pm3.2$d13C.devInner<- Pm3.2$d13C.inner-Pm3.2$d13C
ggplot()+
  geom_point(aes(as.factor(Sample.Number),d13C.devInner, shape=Layer), size=3, data=Pm3.2) + 
  scale_shape_manual(values=c(1,17,4))+
  #geom_hline(yintercept=0.45)+
  xlab("Sample")+
  ylab(expression(paste(delta^13, "C (\u2030) ", "Deviation from Inner Layer")))+
  labs(Shape="Layer")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45))

Pm3.2$d15N.devInner<- Pm3.2$d15N.inner-Pm3.2$d15N
ggplot()+
  geom_point(aes(as.factor(Sample.Number),d15N.devInner, shape=Layer), size=3, data=Pm3.2) + 
  scale_shape_manual(values=c(1,17,4))+
  #geom_hline(yintercept=0.45)+
  xlab("Sample")+
  ylab(expression(paste(delta^15, "N (\u2030) ", "Deviation from Inner Layer")))+
  labs(Shape="Layer")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45))


## - NEW FINAL PLOTS WITH DEVIATION - FIG 5 (was 7 & 8) ### 
## Inner - Middle, and Middle - Outer, with 0 as the "inner" horizontal line

#First subset and reshape the data so it is "wide" by sample number rather than "long"
Pm3.3 <- Pm3.2[,c(3,5,8,9)]
View(Pm3.3)
Pm3.4<-reshape(Pm3.3, idvar="Sample.Number", timevar="Layer", direction="wide")
View(Pm3.4)
#Then add columns for Inner - Middle and Middle - Outer
Pm3.4$d13C.IM<-Pm3.4$d13C.Middle-Pm3.4$d13C.Inner
Pm3.4$d13C.MO<-Pm3.4$d13C.Outer-Pm3.4$d13C.Middle
Pm3.4$d13C.IO<-Pm3.4$d13C.Outer-Pm3.4$d13C.Inner
mean(Pm3.4$d13C.Middle-Pm3.4$d13C.Inner)
mean(Pm3.4$d13C.Outer-Pm3.4$d13C.Middle)
mean(Pm3.4$d13C.IO)
mean(Pm3.4$d13C.Middle)
mean(Pm3.4$d13C.Outer)
mean(Pm3.4$d13C.Inner)
shapes<-c("triangle"=17, "X"=4)
p8<- ggplot()+
  geom_point(aes(as.factor(Sample.Number),d13C.IM),shape=17,size=3, data=Pm3.4) +
  geom_point(aes(as.factor(Sample.Number),d13C.MO),shape=4,size=3, data=Pm3.4) +
  geom_hline(yintercept=0.0)+
  geom_hline(yintercept=-0.3789286, linetype="longdash")+
  geom_hline(yintercept=-0.1607143, linetype="dotted")+
  xlab("Sample Number")+
  ylab(expression(paste(delta^13, "C (\u2030) ", "Sequential Deviation")))+
  scale_shape_manual(name="legend",breaks=c(17,4),values=shapes,labels=c("Inner-Middle", "Middle-Outer"))+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(axis.text.x = element_text(angle=45))

Pm3.4$d15N.IM<-Pm3.4$d15N.Middle-Pm3.4$d15N.Inner
Pm3.4$d15N.MO<-Pm3.4$d15N.Outer-Pm3.4$d15N.Middle
Pm3.4$d15N.IO<-Pm3.4$d15N.Outer-Pm3.4$d15N.Inner
mean(Pm3.4$d15N.Middle-Pm3.4$d15N.Inner)
mean(Pm3.4$d15N.Outer-Pm3.4$d15N.Middle)
mean(Pm3.4$d15N.IO)
mean(Pm3.4$d15N.Middle)
mean(Pm3.4$d15N.Outer)
mean(Pm3.4$d15N.Inner)
shapes<-c("triangle"=17, "X"=4)
p7<- ggplot()+
  geom_point(aes(as.factor(Sample.Number),d15N.IM),shape=17,size=3, data=Pm3.4) +
  geom_point(aes(as.factor(Sample.Number),d15N.MO),shape=4,size=3, data=Pm3.4) +
  geom_hline(yintercept=0.0)+
  geom_hline(yintercept=-0.1971429, linetype="longdash")+
  geom_hline(yintercept=0.07965286, linetype="dotted")+
  xlab("Sample Number")+
  ylab(expression(paste(delta^15, "N (\u2030) ", "Sequential Deviation")))+
  scale_shape_manual(name="legend",breaks=c(17,4),values=shapes, labels=c("Inner-Middle", "Middle-Outer"))+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig5_Revision.tiff", height = 18, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggarrange(p7,p8, labels = c("A", "B"), ncol=1, nrow=2)
dev.off()

multiplot(p7,p8, cols=1)
ggarrange(p7,p8, labels = c("A", "B"), ncol=1, nrow=2)

citation("lme4")