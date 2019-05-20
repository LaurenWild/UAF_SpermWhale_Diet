# This code was developed to analyze variability in stable isotope ratios of sperm whales, and their groundfish/squid prey. 
# This data is part of Chapter 2 of my Ph.D. dissertation.  
# Species include sperm whales, sablefish, grenadier, shortraker rockfish, spiny dogfish, skates, robust clubhook squid, magister armhook squid, glass squid, and neocalanus copepods.
# Author: Lauren Wild, lawild@alaska.edu
# March 2019

############################################################
# Load necessary libraries:
library(ggplot2)
library(viridis) # color scheme for plots that is easy to read (from simmr)
library(ggplot2)
library(Hmisc)
library(psych)
library(devtools)
library(plotrix)  #library for the std.error function; can use mean_se() as well, get same answer
library(MASS)
library(car)
library(stats)
library(ggpubr)
library(plyr)
library(here)

#### Begin with all layers, Line 32; includes all isotope data for sperm whales with innner layer available
#### Isolate inner layer, begins line 87
#### Prey starts at line 663
#### TL calculations at 1500; will change..

#################################################################################
#### -----------------------------------------------------------------###########
###### --------------------  SPERM WHALE LAYER DATA: -------------------- #######
###### ----------- Contains Inner Layer for isotope analysis ----------- ######## 
#################################################################################
Pm2<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmIsotopes2_forR.csv',sep=",",header=TRUE)
View(Pm2) 

Pm2<-read.table(here::here("PmIsotopes2_forR.csv"),sep=",",header=TRUE)

str(Pm2) #Get the structure of variables (factor, numeric, etc.)
#Need to make Layer a numeric variable: 
Pm2["Layer2"] <- NA #Creates new column
Pm2$Layer2<- 5-as.numeric(Pm2$Layer) # Creates a new column making layer numeric
levels(Pm2$Layer)
str(Pm2$Layer2)
levels(Pm2$Month)

#Explore how all the data looks ... 
ggplot(Pm2,aes(d13C,d15N,label=Layer, color=as.factor(Whale.ID))) + 
  geom_point(size=4)  +  #scale_color_manual(breaks = Pm2$Whale,values=ccodes) +
  geom_text(hjust=0.5,vjust = -2, size = 3)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.position="none")

#Test for normality
hist(Pm2$d13C) 
shapiro.test(Pm2$d13C) #p=0.09494, not sig so it IS normal
hist(Pm2$d15N)
shapiro.test(Pm2$d15N) # p=0.1935, not sig so it IS normal

# Test layers with respect to d13C, and d15N:
# Exclude the "full" sample, just compare layers
fit1.C<- lme(d13C~Layer, data=Pm2[Pm2$Layer!='Full',], random=~1|Sample.Number)
summary(fit1.C) #All layers significantly different
summary(aov(fit1.C))
plot(fit1.C)
fit2.C<-lmer(d13C~Layer + (1 | Sample.Number), data=Pm2[Pm2$Layer!='Full',])
summary(fit2.C)

fit1.N<- lme(d15N~Layer, data=Pm2[Pm2$Layer!='Full',], random=~1|Sample.Number)
summary(fit1.N) #All layers significantly different
summary(aov(fit1.N)) #Now "layer" is not significant?????? 
plot(fit1.N)

ggplot(Pm2,aes(as.factor(Sample.Number),d15N, color=as.factor(Layer))) + 
  geom_point()+
  xlab("Sample Number")+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  labs(color="Layer")
ggplot(Pm2,aes(as.factor(Sample.Number), d13C, color=as.factor(Layer))) + 
  geom_point()+
  xlab("Sample Number")+
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  labs(color="Layer")
######################################################################

######################################################################
### Isolate Inner Layer ### SAMPLE SIZE = 33
######################################################################
Pm2Inner<- Pm2[Pm2$Layer=="Inner",]
dim(Pm2Inner) #length of data set is 33 samples

# Set up month as a numeric value
match(Pm2Inner$Month, month.abb)
sapply(Pm2Inner$Month,function(x) grep(paste("(?i)",x,sep=""),month.abb))
Pm2Inner$Month2<- match(Pm2Inner$Month, month.abb)

#Convert date to julian date and add column to data frame. 
library(date)
class(Pm2Inner$Date) #It's a factor variable
Pm2Inner$Date2<- as.POSIXlt(Pm2Inner$Date, format='%m/%d/%y')
class(Pm2Inner$Date2) #Now class is POSIXlt POSIXt
View(Pm2Inner) #Should have dates in Date2 Column! 
doy <- strftime(Pm2Inner$Date2, format = "%j") # new vector of dates in julian format
doy<-as.numeric(doy) #make julian date numeric
Pm2Inner$doy<-doy #add julian date to data frame
View(Pm2Inner) # Check everything looks good! 

#Plot out all sperm whale inner layers (n=33): 
ggplot(Pm2Inner,aes(d13C,d15N)) + 
  geom_point(size=4)  +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 20, base_family = "Helvetica")

# Explore how d13C and d15N values of inner layers change by month: 
setwd('/Users/laurenwild/Desktop')
tiff(filename="InnerLayerByMonth.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggplot(Pm2Inner,aes(d13C,d15N,label=Month, color=as.factor(Month))) + 
  geom_point(size=4)  +  #scale_color_manual(breaks = Pm2$Whale,values=ccodes) +
  geom_text(hjust=0.5,vjust = -2, size = 2)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(color="Month") #takes off as.factor from legend!
dev.off()

MonthLabs<-c("May", "Jun", "Jul", "Aug", "Sep")
I.fit1.1 <- aov(d13C~Month, data=Pm2Inner)
summary(I.fit1.1) #Significant @ 5% level p=0.018
ggplot(Pm2Inner, aes(as.factor(Month2), d13C))+
  geom_boxplot()+
  scale_x_discrete(labels=c('5'='May', '6'='Jun', '7'='Jul', '8'='Aug', '9'='Sep'))+
  xlab("Month")+
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw()
ggplot(Pm2Inner, aes(as.factor(Month2), d13C)) + geom_boxplot() +
  xlab("Month")+
  ylab(expression(paste(delta^13, "C (\u2030)",sep=""))) +
  scale_x_discrete(labels=MonthLabs)+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

I.fit1.2 <- aov(d15N~Month, data=Pm2Inner)
summary(I.fit1.2) #Not significant  p=0.051
ggplot(Pm2Inner, aes(as.factor(Month2), d15N))+
  geom_boxplot()+
  scale_x_discrete(labels=c('5'='May', '6'='Jun', '7'='Jul', '8'='Aug', '9'='Sep'))+
  xlab("Month")+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw()
ggplot(Pm2Inner, aes(as.factor(Month2), d15N)) + geom_boxplot() +
  xlab("Month")+
  ylab(expression(paste(delta^15, "N (\u2030)",sep=""))) +
  scale_x_discrete(labels=MonthLabs)+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


##Pairwise t-test to test if different months have sig diff in d13C & d15N values;
## NEED TO LOOK UP p.adj TREATMENTS
tapply(Pm2Inner$d13C,INDEX=list(Pm2Inner$Month), FUN=mean)
pairwise.t.test(Pm2Inner$d13C, Pm2Inner$Month, p.adj="none")
tapply(Pm2Inner$d15N,INDEX=list(Pm2Inner$Month), FUN=mean)
pairwise.t.test(Pm2Inner$d15N, Pm2Inner$Month, p.adj="none")

###Plot ISOTOPE RATIOS by YEAR to see if any potential baseline shifts over the sampling period
ggplot(Pm2Inner, aes(as.factor(Year), d13C)) + geom_boxplot() +
  xlab("Year")+
  ylab(expression(paste(delta^13, "C (\u2030)",sep=""))) +
  #scale_x_discrete(labels=MonthLabs)+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
I.fit1.3 <- aov(d13C~Year, data=Pm2Inner)
summary(I.fit1.3) #Significant at 5% level p=0.0396
p<-plot(as.factor(Pm2Inner$Year), Pm2Inner$d13C, main="Inner Layer d13C by Year", xlab="Year", ylab=expression(paste(delta^13, "C (\u2030)",sep="")))
plot(as.factor(Pm2Inner$Year), Pm2Inner$d13C, main="Inner Layer d13C by Year", xlab="Year", ylab=expression(paste(delta^13, "C (\u2030)",sep="")))   
text(seq_along(Pm2Inner), p$stats[3,], p$n) 
#Probably not enough data points per year ....

ggplot(Pm2Inner, aes(as.factor(Year), d15N)) + geom_boxplot() +
  xlab("Year")+
  ylab(expression(paste(delta^15, "N (\u2030)",sep=""))) +
  #scale_x_discrete(labels=MonthLabs)+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
I.fit1.4 <- aov(d15N~Year, data=Pm2Inner)
summary(I.fit1.4) #Not significant
q<-plot(as.factor(Pm2Inner$Year), Pm2Inner$d15N,main="Inner Layer d15N by Year", xlab="Year", ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
plot(as.factor(Pm2Inner$Year), Pm2Inner$d15N,main="Inner Layer d15N by Year", xlab="Year", ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
text(seq_along(Pm2Inner), q$stats[3,], q$n) 

#SET UP AVG YEARLY ISOTOPE RATIOS FOR INNER LAYER SAMPLES BIPLOT
PmInnerAvg <- data.frame(cbind(Pm2Inner$d15N, Pm2Inner$d13C))
d15N.In.Mean <- aggregate(d15N~Year, data=Pm2Inner, mean)
d13C.In.Mean <- aggregate(d13C~Year, data=Pm2Inner, mean)
d15N.In.SD <- aggregate(d15N~Year, data=Pm2Inner, sd)
d13C.In.SD <- aggregate(d13C~Year, data=Pm2Inner, sd)
PmInnerAvgYr <- data.frame(cbind(d15N.In.Mean,d13C.In.Mean,d15N.In.SD, d13C.In.SD))
View(PmInnerAvgYr)
PmInnerAvgYr$Year.1=NULL
PmInnerAvgYr$Year.2=NULL
PmInnerAvgYr$Year.3=NULL
colnames(PmInnerAvgYr) <- c("Year", "d15N", "d13C", "SD.N", "SD.C")
View(PmInnerAvgYr)

setwd('/Users/laurenwild/Desktop')
tiff(filename="PmInner_By_Year.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
ggplot(PmInnerAvgYr,aes(d13C,d15N, label=Year)) + geom_point(size=4) +
  geom_errorbarh(aes(xmax=PmInnerAvgYr$d13C+PmInnerAvgYr$SD.C,xmin=PmInnerAvgYr$d13C-PmInnerAvgYr$SD.C, height = 0.01)) +
  geom_errorbar(aes(ymax=PmInnerAvgYr$d15N+PmInnerAvgYr$SD.N,ymin=PmInnerAvgYr$d15N-PmInnerAvgYr$SD.N, width = 0.01))+
  geom_text(hjust=0.5,vjust = -2, size = 3)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Year") #takes off as.factor from legend! 
dev.off()

ggplot(Pm2Inner,aes(d13C,d15N, color=as.factor(Year), label=Year)) + geom_point(size=4) +
  geom_text(hjust=0.5,vjust = -2, size = 3)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Year") #takes off as.factor from legend!

ggplot(Pm2Inner,aes(d13C,d15N,color=Season)) + 
  geom_point(size=4)  +  #scale_color_manual(breaks = Pm2$Whale,values=ccodes) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d(breaks=c("Early","Mid","Late"))+
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(color="Season") #takes off as.factor from legend!

#Plot Frequent versus non-frequent depredators
ggplot(Pm2Inner,aes(d13C,d15N, color=as.factor(Frequent))) + geom_point(size=4) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_manual(values=c("Frequent"="purple", "Non-Frequent"="blue", "Unk"="grey"))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Frequent") + theme(legend.position = c(0.15, 0.8))+
  theme(legend.title=element_blank(),legend.text=element_text(size=8)) #legend.direction = "horizontal"

Pm2InnerSerial<-subset(Pm2Inner, !Frequent=="Unk")
ggplot(Pm2InnerSerial,aes(d13C,d15N, color=as.factor(Frequent))) + geom_point(size=4) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_manual(values=c("Frequent"="purple", "Non-Frequent"="blue"))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Frequent") + theme(legend.position = c(0.15, 0.8))+
  theme(legend.title=element_blank(),legend.text=element_text(size=8)) #legend.direction = "horizontal"

ggplot(Pm2Inner,aes(d13C,d15N, color=as.factor(Recent))) + geom_point(size=4) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Recent") + theme(legend.position = c(0.15, 0.8))+
  theme(legend.title=element_blank(),legend.text=element_text(size=8))+ #legend.direction = "horizontal"
  scale_color_discrete(breaks=c("Old", "Recent"),
                       labels=c("Old (<2010)", "Recent (>2010)"))

Pm2InnerRecent<-subset(Pm2Inner, Recent=="Recent")
ggplot(Pm2InnerRecent,aes(d13C,d15N, color=as.factor(Frequent))) + geom_point(size=4) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_manual(values=c("Frequent"="purple", "Non-Frequent"="blue"))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Frequent") + theme(legend.position = c(0.85, 0.2))+
  theme(legend.title=element_blank(),legend.text=element_text(size=8)) #legend.direction = "horizontal"

InnerOldRecentC<-aov(d13C~Recent, data=Pm2Inner) #sig @ 5% level recent vs. old for d13C
InnerOldRecentN<-aov(d15N~Recent, data=Pm2Inner) #not sig recent vs. old for d15N

fitdoyC<-lm(d13C~doy, data=Pm2Inner) #doy not significant
fitdoyN<-lm(d15N~doy, data=Pm2Inner) #doy not significant

write.table(Pm2Inner, file="PmSIMMData.csv", sep=",")

#################################################################
### Mixed Effects models for d15N:  DAY OF YEAR, REGION, YEAR ###
#################################################################
mod1<-lme(d15N ~ doy + I(doy^2) + as.factor(Region), data=Pm2Inner, random = ~ 1 |as.factor(Year))
mod2<-lme(d15N ~ doy + as.factor(Region), data=Pm2Inner, random = ~ 1|as.factor(Year))
mod3<-lme(d15N ~ doy + I(doy^2) + as.factor(Region) + Year, data=Pm2Inner, random = ~ 1 |as.factor(Sample.Number))
mod4<-lme(d15N ~ doy + as.factor(Region) + Year, data=Pm2Inner, random= ~ 1|as.factor(Sample.Number))
mod5<-lm(d15N ~ doy + as.factor(Region) + Year, data=Pm2Inner)
anova(mod1) #doy^2 term is significant (p=0.0297), which means there is slightly significant seasonal variation. 
anova(mod2) #neither variable significant
anova(mod3) #doy^2 significant p=0.0085; Year almost p=0.0659
anova(mod4) #no variables significant
anova(mod5) #no variables significant
dredge(mod1) #best model is intercept-only model; AICc=86.3, next lowest 92.5(w/ region)
dredge(mod2) #best model is intercept-only model; AICc=86.3, next lowest 92.5(w/ region)
dredge(mod3) #best model is intercept-only model; AICc=87.5, next lowest 94.5(w/ region)
dredge(mod4) #best model is intercept-only model; AICc=87.5, next lowest 94.5(w/ region)

setwd('/Users/laurenwild/Desktop')
tiff(filename="DOY_d15N_Plot.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 300)
ggplot(aes(doy,d15N), data=Pm2Inner)+
  geom_point() + 
  #geom_line(aes(x=Pm2Inner$doy, y), data=Pm2Inner) +
  #geom_line(aes(doy,fitted), data=Pm2Inner)
  stat_smooth(method="lm", formula = y~ x + I(x^2), size = 0.5)+
  xlab("Julian Day")+ 
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw()
dev.off()

plot(mod1)
plot(mod2)
plot(Pm2Inner$doy, residuals(mod1)) #look @ residuals over time to make sure they look random
plot(Pm2Inner$doy, residuals(mod2))

#################################################################
### Mixed Effects models for d13C:  DAY OF YEAR, REGION, YEAR ###
#################################################################
mod6<-lme(d13C ~ doy + I(doy^2) + as.factor(Region), data=Pm2Inner, random = ~ 1 |as.factor(Year))
mod7<-lme(d13C ~ doy + as.factor(Region), data=Pm2Inner, random = ~ 1|as.factor(Year))
mod8<-lme(d13C ~ doy + I(doy^2) + as.factor(Region) + Year, data=Pm2Inner, random = ~ 1 |as.factor(Sample.Number))
mod9<-lme(d13C ~ doy + as.factor(Region) + Year, data=Pm2Inner, random= ~ 1|as.factor(Sample.Number))
anova(mod6) #no variables significant (doy is close, p=0.086). 
anova(mod7) #neither variable significant (doy is close, p=0.0998)
anova(mod8) #doy is significant (p=0.0239) and region (p=0.0463) @ 5% level
anova(mod9) #no variables significant
dredge(mod6) #best model is intercept-only model; AICc = 65.0
dredge(mod7) #best model is intercept-only model; AICc=65.0
dredge(mod8) #best model is intercept-only model; AICc=70.0
dredge(mod9) #best model is intercept-only model; AICc=70.0

ggplot(aes(doy,d13C), data=Pm2Inner)+
  geom_point() + 
  #geom_line(aes(x=Pm2Inner$doy, y), data=Pm2Inner) +
  #geom_line(aes(doy,fitted), data=Pm2Inner)
  stat_smooth(method="lm", formula = y~ x + I(x^2), size = 0.5)+
  xlab("Julian Day")+ 
  ylab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  theme_bw()

#########################################################################
##### BUILD SEPARATE Pm DATAFRAMES FOR SERIAL, NONSERIAL, RECENT, & OLD
##### To input later for use with Mixing Models #######
#########################################################################
#keep<- c("GOA.091", "GOA-091", "GOA-064", "GOA-010", "GOA-085", "GOA-026")
#Pm2InnerSerial<-Pm2Inner[Pm2Inner$Whale.ID %in% keep,]
Pm2InnerSerial<-subset(Pm2Inner, Frequent=="Frequent")
View(Pm2InnerSerial)
write.table(Pm2InnerSerial, file="PmSISerial.csv", sep=",")

#Pm2InnerNonSerial<- Pm2Inner[!Pm2Inner$Whale.ID %in% keep,]
Pm2InnerNonSerial<- subset(Pm2Inner, Frequent=="Non-Frequent")
View(Pm2InnerNonSerial)
write.table(Pm2InnerNonSerial, file="PmSINonSerial.csv", sep=",")

Pm2InnerRecent<-Pm2Inner[Pm2Inner$Year>2009,]
View(Pm2InnerRecent)
write.table(Pm2InnerRecent, file="PmSIRecent.csv", sep=",")

Pm2InnerOld<-Pm2Inner[Pm2Inner$Year<2010,]
View(Pm2InnerOld)
write.table(Pm2InnerOld, file="PmSIOld.csv", sep=",")

###############################################################
## -------------------- PREY ------------------------------- ##
##### --------------------------------------------------- #####

Prey2<-read.csv('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/Prey.LE.final.Isotopes.forR4.csv',header=TRUE, sep=",")
View(Prey2)

#Subest each species
Or<-subset(Prey2, Species=='Clubhook Squid')
nrow(Or)
Ia<-subset(Prey2, Species=='Ragfish')
nrow(Ia)
Af<-subset(Prey2, Species=="Sablefish")
nrow(Af)
Cy<-subset(Prey2, Species=='Grenadier')
nrow(Cy)
Ap<-subset(Prey2, Sub.Species=='Giant Grenadier')
nrow(Ap)
Cy2<-subset(Prey2, Sub.Species=='Grenadier')
nrow(Cy2)
Sb<-subset(Prey2, Species=="Shortraker Rockfish")
nrow(Sb)
Sa<-subset(Prey2, Species=="Spiny Dogfish")
nrow(Sa)
Rb<-subset(Prey2, Species=="Skate")
nrow(Rb)
Rr<-subset(Prey2, Sub.Species=="Longnose Skate")
nrow(Rr)
Bm<-subset(Prey2, Species=='Magister Squid')
nrow(Bm)
Gp<-subset(Prey2, Species=='Glass Squid')
nrow(Gp)

#Calculate differences between bulk & LE, to determine if there is difference in LE; 
Or$D.d15N<-Or$d15N.bulk - Or$d15N.LE
Or$D.d13C<-Or$d13C.bulk - Or$d13C.LE
Cy$D.d15N<-Cy$d15N.bulk - Cy$d15N.LE
Cy$D.d13C<-Cy$d13C.bulk - Cy$d13C.LE
Af$D.d15N<-Af$d15N.bulk - Af$d15N.LE
Af$D.d13C<-Af$d13C.bulk - Af$d13C.LE
Sb$D.d15N<-Sb$d15N.bulk - Sb$d15N.LE
Sb$D.d13C<-Sb$d13C.bulk - Sb$d13C.LE
Bm$D.d15N<-Bm$d15N.bulk - Bm$d15N.LE
Bm$D.d13C<-Bm$d13C.bulk - Bm$d13C.LE
Sa$D.d15N<-Sa$d15N.bulk - Sa$d15N.LE
Sa$D.d13C<-Sa$d13C.bulk - Sa$d13C.LE
Rb$D.d15N<-Rb$d15N.bulk - Rb$d15N.LE
Rb$D.d13C<-Rb$d13C.bulk - Rb$d13C.LE
Gp$D.d15N<-Gp$d15N.bulk - Gp$d15N.LE
Gp$D.d13C<-Gp$d13C.bulk - Gp$d13C.LE
Ia$D.d15N<-Ia$d15N.bulk - Ia$d15N.LE
Ia$D.d13C<-Ia$d13C.bulk - Ia$d13C.LE

####------------------------------------------
#### PLOTS to see effect of LE on d15N values: 
####------------------------------------------
#Onykia plots to look at LE vs NLE effect on d15N, d13C, and C:N values: 
ggplot(Or, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Or.N <- lm(d15N.bulk~d15N.LE, data=Or)
summary(Or.N)
t.test(Or$D.d15N) #p=0.049, not significant @ 5% level, so bulk N is no different from LE. 

ggplot(Or, aes(d13C.bulk, d13C.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^13, "C (\u2030)", " bulk")))+
  ylab(expression(paste(delta^13, "C (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Or.C <- lm(d13C.bulk~d13C.LE, data=Or)
summary(Or.C)
t.test(Or$D.d13C) #p=0.004 Confirms big difference in LE vs NLE for carbon, justifies lipid-extracting

#Grenadier Plots to look at LE vs NLE of d15N, d13C, and C:N values
ggplot(Cy, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Cy.N<-lm(d15N.bulk ~ d15N.LE, data=Cy)
summary(Cy.N) #p=<0.0001, so d15N is sig different between LE & NLE
t.test(Cy$D.d15N) #p=0.0004 - so should probably LE 1/2 and NLE other half

#Sablefish Plots to look at LE vs NLE of d15N values
ggplot(Af, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  #scale_x_continuous(name=expression(paste(delta^15, "N (\u2030)", " bulk")), breaks=seq(12.5,18,1), limits=c(12.5,18))+
  #scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)", " extracted")), breaks=seq(12.5,18.5,1), limits=c(12.5,18.5))+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Af.N <- lm(d15N.bulk~d15N.LE, data=Af)
summary(Af.N)
t.test(Af$D.d15N)  #p<0.0001 - so should LE for C and NLE for N

ggplot(Af, aes(d13C.bulk, d13C.LE)) + geom_point()+
  geom_smooth(method='lm')+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)", " bulk")), breaks=seq(-22,-16,1), limits=c(-22,-16))+
  scale_y_continuous(name=expression(paste(delta^13, "C (\u2030)", " extracted")), breaks=seq(-22,-16,1), limits=c(-22,-16))+
  xlab(expression(paste(delta^13, "C (\u2030)", " bulk")))+
  ylab(expression(paste(delta^13, "C (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Af.C <- lm(d13C.bulk~d13C.LE, data=Af)
summary(Af.C)
t.test(Af$D.d13C) #p=7.98e-09, sig diff btwn LE & NLE, so need to LE

#Rockfish Plots to look at LE vs NLE of d15N values
ggplot(Sb, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Sb.N<-lm(d15N.bulk ~ d15N.LE, data=Sb)
summary(Sb.N)
t.test(Sb$D.d15N) #p<0.0001 - so should LE for C and NLE for N

## Berryteuthis Plots to look at LE vs NLE of d15N 
ggplot(Bm, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Bm.N <- lm(d15N.bulk~d15N.LE, data=Bm)
summary(Bm.N)
t.test(Bm$D.d15N) #p<0.0001, so need to use bulk d15N 

ggplot(Bm, aes(d13C.LE, d15N.bulk)) + geom_point() + 
  xlab(expression(paste(delta^13, "C (\u2030)"))) +
  ylab(expression(paste(delta^15, "N (\u2030)"))) +
  theme_bw(base_size = 24, base_family = "Helvetica")

## Spiny Dogfish Plots to look at LE vs NLE of d15N
ggplot(Sa, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)  #Plot really looks off - like we need to use bulk d15N
Sa.N <- lm(d15N.bulk~d15N.LE, data=Sa)
summary(Sa.N)
t.test(Sa$D.d15N) #p=0.226, so don't need to use bulk d15N

## Skate Plots to look at LE vs NLE of d15N
ggplot(Rb, aes(d15N.bulk, d15N.LE)) + geom_point()+
  geom_smooth(method='lm')+
  xlab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " extracted")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  geom_abline(intercept=0)
Rb.N <- lm(d15N.bulk~d15N.LE, data=Rb)
summary(Rb.N)
t.test(Rb$D.d15N) #p=0.988, so don't need to use bulk d15N

### Ragfish
ggplot(Ia, aes(d13C.LE, d15N.bulk)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
#PWS shallower, but larger, ragfish are 1permil higher in both d13C & d15N...

###########################################################
### Explore each species plot for outliers:
###########################################################
ggplot(Af, aes(d13C.LE, d15N.bulk)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none') #one potential outlier more negative d13C,row 162

ggplot(Cy, aes(d13C.LE, d15N.bulk, color=Sub.Species)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none') #no outliers
#theme(legend.position=c(0.2,0.82))
## plot shows Giant and Pacific Grenadier have similar isotope ratios, and it's a huge spread

ggplot(Sb, aes(d13C.LE, d15N.bulk)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none') #one outlier in low d15N,  low d13C, row 256

ggplot(Sa, aes(d13C.LE, d15N.bulk)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica") #no outliers for bulk values

ggplot(Rb, aes(d13C.LE, d15N.LE, color=Sub.Species)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
## plot shows longnose and skate have similar isotope ratios, and it's a huge spread in Nitrogen; 
#two outliers in d13C, lines 310 & 311
######--------------------------------
######--------------------------------
### DEAL WITH OUTLIERS:
Prey3<- Prey2[-c(55, 162, 254, 318, 324),] #all outliers removed. 
#Prey2<- Prey2[-c(55),] #one outlier with grenadier length 45
#Prey2<- Prey2[-c(310,311),] #two outliers in skates with very negative d13C values
#Prey2<- Prey2[-c(256),] #one outlier in Shortraker rockfish less than 12permil, row 256
#Prey2<- Prey2[-c(162,170),] #two sablefish outliers
#Prey2<- Prey2[-c(382,383,384),] #three outliers for spiny dogfish, low d15N LE
View(Prey3)
write.table(Prey3, file="Prey-outliers-removed-Sep2018.csv", sep=",")

Prey3<-read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/Prey-outliers-removed-Sep2018.csv',sep=",",header=TRUE)
View(Prey3)
Prey3$Depth.Strata<-as.factor(Prey3$Depth.Strata) #Needs to be factor for significance testing

Or<-subset(Prey3, Species=='Clubhook Squid')
Ia<-subset(Prey3, Species=='Ragfish')
Af<-subset(Prey3, Species=="Sablefish")
Cy<-subset(Prey3, Species=='Grenadier')
Ap<-subset(Prey3, Sub.Species=='Giant Grenadier')
Cy2<-subset(Prey3, Sub.Species=='Grenadier')
Sb<-subset(Prey3, Species=="Shortraker Rockfish")
Sa<-subset(Prey3, Species=="Spiny Dogfish")
Rb<-subset(Prey3, Species=="Skate")
Rr<-subset(Prey3, Sub.Species=="Longnose Skate")
Bm<-subset(Prey3, Species=='Magister Squid')
Gp<-subset(Prey3, Species=='Glass Squid')

###--------------------------------------------
#Species that need separate bulk N and LE for C: Sablefish, Grenadier, Shortraker, Magister Squid, Glass Squid, Ragfish, Spiny Dogfish; Species that don't: Skate, clubhook squid
#Need to make columns that have d13C LE values and d15N bulk values where necessary
#Also d13C LE values and d15N LE values where allowed as well...
#For species that have bulk values, reduce spreadsheet to just the rows where bulk was done:
Af2<-Af[!is.na(Af$d15N.bulk),]
Af2$d15N<-Af2$d15N.bulk
nrow(Af2) #45
ncol(Af2) #30
Cy2<-Cy[!is.na(Cy$d15N.bulk),]
Cy2$d15N<-Cy2$d15N.bulk
nrow(Cy2) #44
ncol(Cy2) #30
Sb2<-Sb[!is.na(Sb$d15N.bulk),]
Sb2$d15N<-Sb2$d15N.bulk
nrow(Sb2) #44
ncol(Sb2) #30
Sa2<-Sa[!is.na(Sa$d15N.bulk),]
Sa2$d15N<-Sa2$d15N.bulk
nrow(Sa2) #34
ncol(Sa2) #30
Bm2<-Bm[!is.na(Bm$d15N.bulk),]
Bm2$d15N<-Bm2$d15N.bulk
nrow(Bm2) #44
ncol(Bm2) #30
Ia$d15N<-Ia$d15N.bulk
nrow(Ia) #3
ncol(Ia) #30
Gp$d15N<-Gp$d15N.bulk
nrow(Gp) #2
ncol(Gp) #30
Or$d15N<-Or$d15N.bulk
nrow(Or) #10
ncol(Or) #30
Rb$d15N<-Rb$d15N.LE
nrow(Rb) #66
ncol(Rb) #30
mean(Ia$d15N.bulk)
sd(Ia$d15N.bulk)
mean(Ia$d13C.LE)
sd(Ia$d13C.LE)

#Combine dogfish & ragfish as a source (b/c simmr only allows 2 sp to be combined, and need Af, Sa, & Ia)
SaIa<-rbind(Sa2, Ia)
nrow(SaIa) #37, 3 ragfish & 34 dogfish
SaIa$Species<-'Dogfish.Ragfish'

###---------------------------------------------------------------------
### Explore each species (length/isotope relationship, subspecies/isotope relationship, etc.)
### 
library(jmv) #Has the mancova function
## GRENADIER #### 
Cy2$d13C.abs<-Cy2$d13C.LE*-1
Grennie<-mancova(data=Cy2, deps=vars(d15N.bulk, d13C.LE), factors=vars(Depth.Strata, Sub.Species), covs=Length)
Grennie<-manova(cbind(d15N.bulk,d13C.LE)~Length+Depth.Strata+Sub.Species, data=Cy2)
summary(Grennie) #only length, 
summary.aov(Grennie) #d13C significant at 1% and d15N almost at 5% level for Length
summary(aov(d15N.bulk~Length, data=Cy2)) #confirmed not over %
summary(aov(d13C.LE~Length, data=Cy2)) #confirmed significant
ggplot(Cy, aes(Length, d15N.bulk, color=Sub.Species)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')  #Length 45cm is an outlier. Remove row 55?
ggplot(Cy, aes(Length, d13C.LE, color=Sub.Species)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Cy, aes(Depth.Strata, d13C.LE, color=Sub.Species)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Cy, aes(Depth.Strata, d15N.bulk, color=Sub.Species)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')

### SKATE #### 
Skate<-manova(cbind(d15N.bulk,d13C.LE)~Length+Depth.Strata+Sub.Species, data=Rb)
summary(Skate) #Length & maybe depth strata significant
summary.aov(Skate) # Depth & d15N at 5% level; p=0.45; Length at d13C, p<0.002
summary(aov(d15N.bulk~Depth.Strata, data=Rb)) #p=0.44
summary(aov(d13C.LE~Length, data=Rb)) #p=0.004
ggplot(Rb, aes(Length, d15N.LE, color=Sub.Species)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Rb, aes(Length, d13C.LE, color=Sub.Species)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Rb, aes(Depth.Strata, d15N.LE, color=Sub.Species)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Rb, aes(Depth.Strata, d13C.LE, color=Sub.Species)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')

## SHORTRAKER ROCKFISH ###
SrRock<-manova(cbind(d15N.bulk,d13C.LE)~Length+Depth.Strata+Year, data=Sb2)
summary(SrRock) # Nothing significant
summary.aov(SrRock) #Nothing Significant
ggplot(Sb, aes(Length, d15N.bulk)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")
ggplot(Sb, aes(Length, d13C.LE)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")
ggplot(Sb, aes(Depth.Strata, d13C.LE)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Sb, aes(Depth.Strata, d15N.bulk)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')

#### SABLEFISH ###
Sable<-manova(cbind(d15N.bulk,d13C.LE)~Length+Depth.Strata, data=Af2)
summary(Sable) # Nothing significant
summary.aov(Sable) #Nothing significant
ggplot(Af, aes(Length, d15N.bulk)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")
ggplot(Af, aes(Length, d13C.LE)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

### SPINY DOGFISH  ##### 
Dogfish<-manova(cbind(d15N.bulk,d13C.LE)~Length+Depth.Strata, data=Sa2)
summary(Dogfish) # Depth.Strata significant
summary.aov(Dogfish) #Length for d15N and Depth for d13C both at 5% level
summary(aov(d15N.bulk~Length, data=Af2)) # at 5% level
summary(aov(d13C.LE~Depth.Strata, data=Af2)) #Not significant
ggplot(Sa, aes(Length, d15N.bulk)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")
ggplot(Sa, aes(Length, d13C.LE)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")
ggplot(Sa, aes(Depth.Strata, d15N.bulk)) + geom_point()+
  xlab("Depth")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")
ggplot(Sa, aes(Depth.Strata, d13C.LE)) + geom_point()+
  xlab("Depth")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")

### BERRYTEUTHIS MAGISTER ####
Berry<-manova(cbind(d15N.bulk,d13C.LE)~Length*Depth.Strata, data=Bm2)
summary(Berry) #Length is significant
summary.aov(Berry) #Length and d15N signficant
summary(aov(d15N.bulk~Length, data=Bm2)) #Significant
summary(aov(d13C.LE~Depth.Strata, data=Bm2)) #significant at 5% level
ggplot(Bm, aes(Length, d15N.bulk)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')  
ggplot(Bm, aes(Length, d13C.LE)) + geom_point()+
  xlab("Length")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Bm, aes(Depth.Strata, d13C.LE)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^13, "C (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')
ggplot(Bm, aes(Depth.Strata, d15N.bulk)) + geom_point()+
  xlab("Depth Strata")+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position='none')


Berry<-manova(cbind(d15N.bulk,d13C.LE)~Length+Depth.Strata, data=Bm2)
summary(Berry) #Length is significant
summary.aov(Berry) #Length and d15N signficant
summary(aov(d15N.bulk~Length, data=Bm2)) #Significant
summary(aov(d13C.LE~Depth.Strata, data=Bm2)) #significant at 5% level


#Then combine everything back to new reduced dataset
Prey3.2<-rbind(Or,Cy2,Af2,Sb2,Sa2,Rb,Bm2,Ia,Gp) 
Prey3.3<-rbind(Or,Cy2,Af2,Sb2,Sa2,Rb,Bm2)
Prey3.4<-rbind(Or,Cy2,Af2,Sb2,Sa2,Rb,Bm2,Gp)
Prey3.5<-rbind(Or,Cy2,Af2,SaIa,Sb2,Rb,Bm2,Gp)

#Now want the LE values of d13C ... 
Prey3.2$d13C<- Prey3.2$d13C.LE
Prey3.3$d13C<- Prey3.3$d13C.LE
Prey3.4$d13C<- Prey3.4$d13C.LE
Prey3.5$d13C<- Prey3.5$d13C.LE

color2<-c("red","forestgreen","cyan1","black","blue","purple","darkslategray3", "darkgoldenrod")
ggplot(Prey3.2, aes(d13C, d15N, color=Species)) + geom_point(size=3)+
  xlab(expression(paste(delta^13, "C (\u2030)", sep='')))+
  ylab(expression(paste(delta^15, "N (\u2030)", sep='')))+
  scale_colour_viridis_d()
color3<-c("red","forestgreen","black","blue","purple","darkslategray3", "darkgoldenrod")
ggplot(Prey3.3, aes(d13C, d15N, color=Species)) + geom_point(size=3)+
  xlab(expression(paste(delta^13, "C (\u2030)", sep='')))+
  ylab(expression(paste(delta^15, "N (\u2030)", sep='')))+
  scale_colour_viridis_d()
#Add elipses
ggplot(Prey3.2, aes(d13C, d15N, color=Species)) + geom_point()+
  xlab(expression(paste(delta^13, "C (\u2030)", sep='')))+
  ylab(expression(paste(delta^15, "N (\u2030)", sep='')))+
  #theme_bw(base_size = 24, base_family = "Helvetica")+
  scale_colour_viridis_d()+
  stat_ellipse(type="norm")

#Put all d15N values in the same column (bulk for those that need bulk, and LE for those without)
#Now the d15N.final column has the accurate number for that species. Same for d13C
#Prey3$d15N<- Prey3$d15N.bulk
#View(Prey3)
#for (i in 1:nrow(Prey3)) {
# Prey3$d15N[is.na(Prey3$d15N)] <- Prey3$d15N.LE[is.na(Prey3$d15N)]
#}
#View(Prey3)
#Prey2.2<-Prey2[!is.na(Prey2$d13C.LE),]
#Prey2.2<-Prey2.2[!is.na(Prey2.2$d15N.LE),]
########################################################
########################################################
library(doBy)
myfun1<-function(x) {c(m=mean(x) , v=var(x))}
myfun2<-function(x) {c(m=mean(x), sd=sd(x))}
All.Prey.Sum<- summaryBy(d15N+d13C~Species, data=Prey3.3, FUN=myfun2) #Just 7 sp
All.Prey.Sum2<- summaryBy(d15N+d13C~Species, data=Prey3.2, FUN=myfun2) #GpIa
All.Prey.Sum3<- summaryBy(d15N+d13C~Species, data=Prey3.4, FUN=myfun2) #Just Gp
All.Prey.Sum4<- summaryBy(d15N+d13C~Species, data=Prey3.5, FUN=myfun2) #Just Gp, SaIa compbined
write.table(All.Prey.Sum, file="PmSources.csv", sep=",")
write.table(All.Prey.Sum2, file="PmSources-GpIa.csv", sep=",")
write.table(All.Prey.Sum3, file="PmSources-Gp.csv", sep=",")
write.table(All.Prey.Sum4, file="PmSources_Gp_IaSaCombined.csv", sep=",")

ggplot(All.Prey.Sum, aes(d13C.m, d15N.m, color=Species)) + geom_point(size=3)+
  #geom_text(hjust=0.3,vjust = -0.7, size = 5)+
  geom_errorbarh(aes(xmax=All.Prey.Sum$d13C.m+All.Prey.Sum$d13C.sd,xmin=All.Prey.Sum$d13C.m-All.Prey.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Prey.Sum$d15N.m+All.Prey.Sum$d15N.sd,ymin=All.Prey.Sum$d15N.m-All.Prey.Sum$d15N.sd, width = 0.01))+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.title=element_text(size=10))+
  theme(legend.text=element_text(size=8))
#theme(legend.position = c(0.9, 0.75))

### This plot has each prey as a different color, and black text within plot
ggplot(All.Prey.Sum, aes(d13C.m, d15N.m, color=Species, label=Species)) + geom_point(size=3)+
  geom_errorbarh(aes(xmax=All.Prey.Sum$d13C.m+All.Prey.Sum$d13C.sd,xmin=All.Prey.Sum$d13C.m-All.Prey.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Prey.Sum$d15N.m+All.Prey.Sum$d15N.sd,ymin=All.Prey.Sum$d15N.m-All.Prey.Sum$d15N.sd, width = 0.01))+
  geom_text(color="black",hjust=-0.05,vjust = -0.7, size = 5)+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")

Shape1<- c(3,8,17,18,7,19,15,13)
ggplot(All.Prey.Sum,aes(d13C.m,d15N.m, label=Species, shape=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Prey.Sum$d13C.m+All.Prey.Sum$d13C.sd,xmin=All.Prey.Sum$d13C.m-All.Prey.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Prey.Sum$d15N.m+All.Prey.Sum$d15N.sd,ymin=All.Prey.Sum$d15N.m-All.Prey.Sum$d15N.sd, width = 0.01))+
  geom_text(hjust=-0.06,vjust = -0.7, size = 5)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_shape_manual(values=Shape1)+
  theme_bw(base_size = 24, base_family = "Helvetica")

### Add glass squid & ragfish prey species; 
### Prey as a different color, and black text within plot
ggplot(All.Prey.Sum2, aes(d13C.m, d15N.m, color=Species, label=Species)) + geom_point(size=3)+
  geom_errorbarh(aes(xmax=All.Prey.Sum2$d13C.m+All.Prey.Sum2$d13C.sd,xmin=All.Prey.Sum2$d13C.m-All.Prey.Sum2$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Prey.Sum2$d15N.m+All.Prey.Sum2$d15N.sd,ymin=All.Prey.Sum2$d15N.m-All.Prey.Sum2$d15N.sd, width = 0.01))+
  geom_text(color="black",hjust=-0.05,vjust = -0.7, size = 5)+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")

#Add humboldt squid for biplot:
All.Prey.Sum5<-read.csv('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources-GpIaDg.csv',sep=",",header=TRUE)
ggplot(All.Prey.Sum5, aes(Mean.d13C, Mean.d15N, color=Species, label=Species)) + geom_point(size=3)+
  geom_errorbarh(aes(xmax=All.Prey.Sum5$Mean.d13C+All.Prey.Sum5$SD.d13C,xmin=All.Prey.Sum5$Mean.d13C-All.Prey.Sum5$SD.d13C, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Prey.Sum4$Mean.d15N+All.Prey.Sum5$SD.d15N,ymin=All.Prey.Sum5$Mean.d15N-All.Prey.Sum5$SD.d15N, width = 0.01))+
  geom_text(color="black",hjust=-0.05,vjust = -0.7, size = 5)+
  xlab(expression(paste(delta^13, "C (\u2030)")))+
  ylab(expression(paste(delta^15, "N (\u2030)")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")
#----------------------------------------------------------------
### Prey data - descriptive stats: 
#Describe range, mean, median, etc. d15N, d13C, and C:N ratios for each species

###---------------------------------------------------------
# Add sperm whales in to the top of it
Pm.Prey <- cbind(Prey3.3)
View(Pm.Prey)
library(gtools)
library(doBy)
#use sperm whale inner layer data frame
Pm2Inner$Species<-"Sperm Whale" 
Pm2Inner2<-Pm2Inner[,c(1,3,4)]
PmInnerAvg<-summaryBy(d15N + d13C ~ Species, data=Pm2Inner2, FUN=myfun2)

Prey3.6<-Prey3.3[,c(1,30,31)]
Prey3.7<-Prey3.2[,c(1,30,31)]
Prey3.8<-Prey3.4[,c(1,30,31)]
Pm.Prey2<-rbind(Prey3.8, Pm2Inner2) #All 7 plus just Gp
Pm.Prey3<-rbind(Prey3.7, Pm2Inner2) #GpIa
Pm.Prey4<-rbind(Prey3.6, Pm2Inner2) #Just 7 main species
str(Pm.Prey2)
str(Pm.Prey3)
str(Pm.Prey4)

#Reduce data frame to just the columns I want:
#keeps <- c("SampleName", "d15N", "d13C", "Species")
#Pm.Prey2<-Pm.Prey2[keeps]

#use doBy to set up data sheet with averages and std.error of each species
myfun2<-function(x) {c(m=mean(x), sd=sd(x))}
All.Sp.Sum <-summaryBy(d15N+d13C~Species, data=Pm.Prey2, FUN=myfun2) #All 7 +Gp 
View(All.Sp.Sum)
All.Sp.Sum2 <-summaryBy(d15N+d13C~Species, data=Pm.Prey3, FUN=myfun2) #All7+GpIa
View(All.Sp.Sum2)
All.Sp.Sum3 <-summaryBy(d15N+d13C~Species, data=Pm.Prey4, FUN=myfun2) #All 7
View(All.Sp.Sum3)
All.Sp.Sum4<- All.Sp.Sum2[-c(5), ] #Another way to take out ragfish... ? 
View(All.Sp.Sum4)
All.Sp.Sum4.5<-rbind(All.Prey.Sum3, PmInnerAvg) #Another way to do it? all 7 sp + Gp

#Another way to add sperm whales in .... 
#PmInnerSums<-summaryBy(d15N+d13C~Layer, data=Pm2Inner, FUN=myfun2)
#colnames(PmInnerSums)<- c("Species", "d15N.m", "d15N.sd", "d13C.m", "d13C.sd")
#PmInnerSums$Species <- as.character(PmInnerSums$Species)
#PmInnerSums$Species[PmInnerSums$Species == "Inner"] <- "Sperm Whale"
#All.Sp.Sum<-rbind(All.Prey.Sum, PmInnerSums)
#View(All.Sp.Sum)

Shape1<- c(3,8,17,18,7,19,15,13)
color3<- c('lightskyblue', 'darkorchid2','blue1','darkturquoise', 'yellow2', 'deeppink',  "coral1", 'darkgoldenrod1')
#Just Main 7 species:
setwd('/Users/laurenwild/Desktop')
tiff(filename="AllSp.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
ggplot(All.Sp.Sum,aes(d13C.m,d15N.m, label=Species, color=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum$d13C.m+All.Sp.Sum$d13C.sd,xmin=All.Sp.Sum$d13C.m-All.Sp.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum$d15N.m+All.Sp.Sum$d15N.sd,ymin=All.Sp.Sum$d15N.m-All.Sp.Sum$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")
dev.off()

# Glass Squid, Ragfish included:
color4<- c('lightskyblue', 'lavenderblush3', 'darkorchid2','blue1','green4', 'darkturquoise', 'yellow2', 'deeppink',  "coral1", 'darkgoldenrod1')
ggplot(All.Sp.Sum2,aes(d13C.m,d15N.m, label=Species, color=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum2$d13C.m+All.Sp.Sum2$d13C.sd,xmin=All.Sp.Sum2$d13C.m-All.Sp.Sum2$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum2$d15N.m+All.Sp.Sum2$d15N.sd,ymin=All.Sp.Sum2$d15N.m-All.Sp.Sum2$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")
setwd('/Users/laurenwild/Desktop')
tiff(filename="AllSp.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
ggplot(All.Sp.Sum2,aes(d13C.m,d15N.m, label=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum2$d13C.m+All.Sp.Sum2$d13C.sd,xmin=All.Sp.Sum2$d13C.m-All.Sp.Sum2$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum2$d15N.m+All.Sp.Sum2$d15N.sd,ymin=All.Sp.Sum2$d15N.m-All.Sp.Sum2$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")
dev.off()

setwd('/Users/laurenwild/Desktop')
tiff(filename="AllSpGp.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
ggplot(All.Sp.Sum,aes(d13C.m,d15N.m, label=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum$d13C.m+All.Sp.Sum$d13C.sd,xmin=All.Sp.Sum$d13C.m-All.Sp.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum$d15N.m+All.Sp.Sum$d15N.sd,ymin=All.Sp.Sum$d15N.m-All.Sp.Sum$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")
dev.off()

color5<- c('lightskyblue', 'lavenderblush3', 'darkorchid2','blue1', 'darkturquoise', 'yellow2', 'deeppink',  "coral1", 'darkgoldenrod1')
ggplot(All.Sp.Sum3,aes(d13C.m,d15N.m, label=Species, color=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum3$d13C.m+All.Sp.Sum3$d13C.sd,xmin=All.Sp.Sum3$d13C.m-All.Sp.Sum3$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum3$d15N.m+All.Sp.Sum3$d15N.sd,ymin=All.Sp.Sum3$d15N.m-All.Sp.Sum3$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")

ggplot(All.Sp.Sum4,aes(d13C.m,d15N.m, label=Species, color=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum4$d13C.m+All.Sp.Sum4$d13C.sd,xmin=All.Sp.Sum4$d13C.m-All.Sp.Sum4$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum4$d15N.m+All.Sp.Sum4$d15N.sd,ymin=All.Sp.Sum4$d15N.m-All.Sp.Sum4$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")

##############################################################################
##### ------- Add Baseline Data, First Import and Reorganize it: ------- #####
##############################################################################

### COPEPODS ### BASELINE ### 
Base<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/Baseline.Isotopes_forR.csv',sep=",",header=TRUE)
View(Base)
#Lipid Normalization to mathematically correct d13C values:
Base$d13C.LN<-(-3.32+0.99*Base$C.N.bulk)+Base$d13C.bulk #Uses Post 2007
Base$d13C.protein<-Base$d13C.bulk + (-6.39*(3.76-Base$C.N.bulk))/Base$C.N.bulk  #Uses Hoffman et al. 2010
Base<-Base[,c(1:5,25,26,6:24)]

range(Base$d15N.bulk)
median(Base$d15N.bulk)
range(Base$d13C.LN)
median(Base$d13C.LN)

ggplot(Base, aes(d13C.LN, d15N.bulk, color=Station)) + geom_point(size=5)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)", " bulk")))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)", " bulk")))+
  xlab(expression(paste(delta^13, "C (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  scale_color_viridis_d()+
  theme(legend.justification= c(1,1), legend.position=c(1,1))

ggplot(Base, aes(d13C.LN, d15N.bulk, shape=Station)) + geom_point(size=3)+
  #geom_text(hjust=-0.06,vjust = -0.7, size = 5)+
  xlab(expression(paste(delta^13, "C (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.justification= c(1,1), legend.position=c(1,1))

#Average each station
library(doBy)
myfun1<-function(x) {c(m=mean(x) , v=var(x))}  
myfun2<-function(x) {c(m=mean(x), sd=sd(x))}
Base.Avg <-summaryBy(d15N.bulk+d13C.LN~Station, data=Base, FUN=myfun2)
View(Base.Avg)

ggplot(Base.Avg, aes(d13C.LN.m, d15N.bulk.m, label=Station)) + geom_point(size=3)+
  geom_text(hjust=-0.06,vjust = -0.7, size = 5)+
  geom_errorbarh(aes(xmax=Base.Avg$d13C.LN.m+Base.Avg$d13C.LN.sd,xmin=Base.Avg$d13C.LN.m-Base.Avg$d13C.LN.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=Base.Avg$d15N.bulk.m+Base.Avg$d15N.bulk.sd,ymin=Base.Avg$d15N.bulk.m-Base.Avg$d15N.bulk.sd, width = 0.01))+
  xlab(expression(paste(delta^13, "C (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.justification= c(1,1), legend.position=c(1,1))

Base.Avg$d15N<-Base.Avg$d15N.bulk.m
Base.Avg$d13C<-Base.Avg$d13C.LN.m
Base.Avg$d15N.sd<-Base.Avg$d15N.bulk.sd
Base.Avg$d13C.sd<-Base.Avg$d13C.LN.sd

Base.Final<-Base.Avg[,c(1,6,8,7,9)]
View(Base.Final)

ggplot(Base.Final, aes(d13C, d15N, label=Station)) + geom_point(size=3)+
  geom_text(hjust=-0.06,vjust = -0.7, size = 5)+
  geom_errorbarh(aes(xmax=Base.Final$d13C+Base.Final$d13C.sd,xmin=Base.Final$d13C-Base.Final$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=Base.Final$d15N+Base.Final$d15N.sd,ymin=Base.Final$d15N-Base.Final$d15N.sd, width = 0.01))+
  xlab(expression(paste(delta^13, "C (\u2030)", " bulk")))+
  ylab(expression(paste(delta^15, "N (\u2030)", " bulk")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.justification= c(1,1), legend.position=c(1,1))

Base.Final3<-summaryBy(d15N+d13C~Species, data=Base.Final, FUN=myfun2)
View(Base.Final3)
#colnames(Base.Final3)<- c("Species", "d15N.m", "d15N.sd", "d13C.m", "d13C.sd")
Base.Final3$Species <- 'Neocalanus sp.'
Base.Final3<-Base.Final3[,c(5,1,2,3,4)]
Base.Final3$d13C.sd<-0.959673
View(Base.Final3)

All.Sp.Base.Sum<-rbind(All.Sp.Sum,Base.Final3)
View(All.Sp.Base.Sum)

color5<- c('lightskyblue', 'darkorchid2','blue1','darkturquoise', 'yellow2', 'deeppink',  "coral1", 'darkgoldenrod1', 'grey2')
ggplot(All.Sp.Base.Sum,aes(d13C.m,d15N.m, label=Species, color=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Base.Sum$d13C.m+All.Sp.Base.Sum$d13C.sd,xmin=All.Sp.Base.Sum$d13C.m-All.Sp.Base.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Base.Sum$d15N.m+All.Sp.Base.Sum$d15N.sd,ymin=All.Sp.Base.Sum$d15N.m-All.Sp.Base.Sum$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_manual(values=color4)+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.position="none")

### ------------------------------------------------------ ###
### -------------------------------------------------------- ###
##### ANNIE'S DATA #### 
myData <- myData[-c(2, 4, 6), ]
All.Sp.Sum.Annie<-All.Sp.Sum[-c(1,2,4,5,6,7),]
Shape3<- c(3,9,17)
setwd('/Users/laurenwild/Desktop')
tiff(filename="AnnieData.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
ggplot(All.Sp.Sum.Annie,aes(d13C.m,d15N.m, label=Species, shape=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=All.Sp.Sum.Annie$d13C.m+All.Sp.Sum.Annie$d13C.sd,xmin=All.Sp.Sum.Annie$d13C.m-All.Sp.Sum.Annie$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum.Annie$d15N.m+All.Sp.Sum.Annie$d15N.sd,ymin=All.Sp.Sum.Annie$d15N.m-All.Sp.Sum.Annie$d15N.sd, width = 0.01))+
  geom_text(hjust=-0.06,vjust = -0.7, size = 3)+
  scale_x_continuous(name=expression(paste(delta^13, "C (\u2030)")), breaks=seq(-26,-16,2), limits=c(-26,-16))+
  scale_y_continuous(name=expression(paste(delta^15, "N (\u2030)")), breaks=seq(7.5,18,2.5), limits=c(7.5,18))+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_shape_manual(values=Shape3)+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.position="none")
dev.off()

###########################################################################
#---------------------------------------------------------------------------
###############################################################################

# Trophic level calculations:

TL<-NULL
for (i in 1:nrow(All.Sp.Base.Sum)) {
  TL[i]<- 2 + ((All.Sp.Base.Sum$d15N.m[i] - All.Sp.Base.Sum$d15N.m[9]) / 2.12) 
}
View(TL)
All.Sp.Base.Sum$TrophicLevel<-TL
write.table(All.Sp.Base.Sum, file="AllSpeciesBaseAvgsTL.csv", sep=",")

#Using Bree's equation, and my 2.12 enrichment factor, and Pm inner layer skin:
Pm2Inner$TL<-2 + (Pm2Inner$d15N - Base.Final3$d15N.m) / 2.12
mean(Pm2Inner$TL)
sd(Pm2Inner$sd)

mean(2 + (Ia$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 4.25 for ragfish
sd(2 + (Ia$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.24
mean(2 + (Af2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 4.59 sablefish
sd(2 + (Af2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.40
mean(2 + (Cy2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 4.45 grenadier
sd(2 + (Cy2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.44
mean(2 + (Sb2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 4.91 shortraker
sd(2 + (Sb2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.0.44
mean(2 + (Sa2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 4.36 dogfish
sd(2 + (Sa2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.44
mean(2 + (Bm2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 3.50 magister squid
sd(2 + (Bm2$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.51
mean(2 + (Rb$d15N - Base.Final3$d15N.m) / 2.12) #TL = 5.36 skate
sd(2 + (Rb$d15N - Base.Final3$d15N.m) / 2.12) #0.32
mean(2 + (Gp$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 3.98 glass squid
sd(2 + (Gp$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.24
mean(2 + (Or$d15N.bulk - Base.Final3$d15N.m) / 2.12) #TL = 5.66 clubhook squid
sd(2 + (Or$d15N.bulk - Base.Final3$d15N.m) / 2.12) #0.38

#Using Geraldine's 1.6 enrichment factor, and Pm Inner Layer skin: 
2 + (PmInnerSums$d15N.m - Base.Final3$d15N.m) / 1.6   #TL = 6.74
#Using Geraldine's 1.6 enrichment factor, and Pm Full: 
2 + (PmFull.mean - Base.Final3$d15N.m) / 1.6  #TL = 6.77


