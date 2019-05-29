# This code was developed to analyze variability in stable isotope ratios of sperm whales, and their groundfish/squid prey. 
# This data is part of Chapter 2 of my Ph.D. dissertation.  
# Species include sperm whales, sablefish, grenadier, shortraker rockfish, spiny dogfish, skates, robust clubhook squid, magister armhook squid, glass squid, and neocalanus copepods.
# Author: Lauren Wild, lawild@alaska.edu; March 2019

library(ggplot2)
library(Hmisc)
library(psych)
library(devtools)
library(plotrix)  #library for the std.error function; can use mean_se() too, same answer
library(MASS)
library(MuMIn) #for dredge
library(lme4)
library(car)
library(nlme)
library(multcomp) #for glht to compare lme parameters
library(AICcmodavg)
library(stats)
library(ggpubr)
library(dplyr)
library(viridis) # color scheme for plots that is easy to read (from simmr)
library(jmv) #Has the mancova function
library(doBy)
library(gtools)
library(here)

#### Begin with all layers, Line 33; includes all isotope data with innner layer available
#### Isolate inner layer, begins line 90
#### Prey starts at line 317
#### TL calculations at 1500; will change..

#################################################################################
#### -----------------------------------------------------------------###########
###### ------ SPERM WHALE SKIN; all samples that include nnner Layer ----- ######
##### --------------------------------------------------------------- ########### 
#################################################################################
Pm2<-read.table(here::here("PmIsotopes2_forR.csv"),sep=",",header=TRUE)
View(Pm2) 

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
#--------------------------------------------------------------------

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

# Explore how repeat samples look: (GOA-010, GOA-024,GOA-027, GOA-045, GOA-091, GOA-136)
n_occur <- data.frame(table(Pm2Inner$Whale.ID))
n_occur[n_occur$Freq > 1,]
Repeat<-Pm2Inner[Pm2Inner$Whale.ID %in% n_occur$Var1[n_occur$Freq > 1],]
Repeat<-Repeat[!(Repeat$Whale.ID=="Unk"),]
ggplot(Repeat,aes(d13C,d15N,label=as.factor(Year), color=as.factor(Whale.ID))) + 
  geom_point(size=4)  +
  geom_text(hjust=0.4,vjust = -0.6, size=4) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(color="Month") #takes off as.factor from legend!

# Explore how d13C and d15N values of inner layers change by month: 
InnerByMonth<-ggplot(Pm2Inner,aes(d13C,d15N,label=Month, color=as.factor(Month))) + 
  geom_point(size=5)  + 
  geom_text(hjust=0.5,vjust = -1, size = 8)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(color="Month") #takes off as.factor from legend!
#setwd('/Users/laurenwild/Desktop')
ggsave(filename="InnerLayerByMonth.png", plot=InnerByMonth, dpi=500, width=17, height=12, units="in")

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
PmInnerAvgYr$Year.1=NULL
PmInnerAvgYr$Year.2=NULL
PmInnerAvgYr$Year.3=NULL
colnames(PmInnerAvgYr) <- c("Year", "d15N", "d13C", "SD.N", "SD.C")
View(PmInnerAvgYr)

InLayerByYear<-ggplot(PmInnerAvgYr,aes(d13C,d15N, label=Year)) + geom_point(size=4) +
  geom_errorbarh(aes(xmax=PmInnerAvgYr$d13C+PmInnerAvgYr$SD.C,xmin=PmInnerAvgYr$d13C-PmInnerAvgYr$SD.C, height = 0.01)) +
  geom_errorbar(aes(ymax=PmInnerAvgYr$d15N+PmInnerAvgYr$SD.N,ymin=PmInnerAvgYr$d15N-PmInnerAvgYr$SD.N, width = 0.01))+
  geom_text(hjust=0.5,vjust = -2, size = 3)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Year") #takes off as.factor from legend! 

#Isotope ratios colored by year
ggplot(Pm2Inner,aes(d13C,d15N, color=as.factor(Year), label=Year)) + geom_point(size=4) +
  geom_text(hjust=0.5,vjust = -2, size = 3)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Year") #takes off as.factor from legend!

#Isotope ratios colored by season
ggplot(Pm2Inner,aes(d13C,d15N,color=Season)) + 
  geom_point(size=4)  +  #scale_color_manual(breaks = Pm2$Whale,values=ccodes) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d(breaks=c("Early","Mid","Late"))+
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(color="Season") #takes off as.factor from legend!

#Isotope ratios for Frequent versus non-frequent depredators
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

#Isotope ratios for old vs. recent samples
ggplot(Pm2Inner,aes(d13C,d15N, color=as.factor(Recent))) + geom_point(size=4) +
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  labs(color="Recent") + theme(legend.position = c(0.15, 0.8))+
  theme(legend.title=element_blank(),legend.text=element_text(size=8))+ #legend.direction = "horizontal"
  scale_color_discrete(breaks=c("Old", "Recent"),
                       labels=c("Old (<2010)", "Recent (>2010)"))

InnerOldRecentC<-aov(d13C~Recent, data=Pm2Inner) #sig @ 5% level recent vs. old for d13C
InnerOldRecentN<-aov(d15N~Recent, data=Pm2Inner) #not sig recent vs. old for d15N

fitdoyC<-lm(d13C~doy, data=Pm2Inner) #doy not significant
fitdoyN<-lm(d15N~doy, data=Pm2Inner) #doy not significant

write.table(Pm2Inner, file="PmSIMMData.csv", sep=",")

#########################################################################
##### BUILD SEPARATE Pm DATAFRAMES FOR SERIAL, NONSERIAL, RECENT, & OLD
##### To input later for use with Mixing Models #######
#########################################################################
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

################################################################
## ---------------------------------------------------------- ##
## -----------------------  PREY  --------------------------- ##
################################################################
# Start with prey spreadsheet that has already had outliers omitted, and testing for how bulk d15N vales may be changed by lipid extraction: 

Prey3<-read.table(here::here('data/Prey-outliers-removed-Sep2018.csv'),sep=",",header=TRUE)
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

# Make sure d15N bulk value is used for each species, and create a final d15N and d13C column
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

###---------------------------------------------------------------------
### Explore each species (length/isotope relationship, subspecies/isotope relationship, etc.)

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

############################################################

#Combine all prey back to new reduced dataset
Prey3.1<-rbind(Or,Cy2,Af2,Sb2,Sa2,Rb,Bm2)

#Now want the LE values of d13C ... 
Prey3.1$d13C<- Prey3.1$d13C.LE

myfun1<-function(x) {c(m=mean(x), sd=sd(x))}
All.Prey.Sum<- summaryBy(d15N+d13C~Species, data=Prey3.1, FUN=myfun1) #Main 7 species

write.table(All.Prey.Sum, file="PmSources.csv", sep=",")

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

#----------------------------------------------------------------

# Add sperm whales in to the top of it
Pm.Prey <- cbind(Prey3.1)

#use sperm whale inner layer data frame
Pm2Inner$Species<-"Sperm Whale" 
Pm2Inner2<-Pm2Inner[,c(1,3,4)]
Prey3.1<-Prey3.1[,c(1,30,31)]
Pm.Prey<-rbind(Prey3.1, Pm2Inner2) #Just 7 main species
str(Pm.Prey)

#set up data sheet with averages and std.error of each species
All.Sp.Sum <-summaryBy(d15N+d13C~Species, data=Pm.Prey, FUN=myfun1) #All 7  
View(All.Sp.Sum)


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

Fig2<-ggplot(All.Sp.Sum,aes(d13C.m,d15N.m, label=Species)) + geom_point(size=8) +
  geom_errorbarh(aes(xmax=All.Sp.Sum$d13C.m+All.Sp.Sum$d13C.sd,xmin=All.Sp.Sum$d13C.m-All.Sp.Sum$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum$d15N.m+All.Sp.Sum$d15N.sd,ymin=All.Sp.Sum$d15N.m-All.Sp.Sum$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 10)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+ 
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=28), axis.title=element_text(size=32,face="bold"))+
  theme(legend.position="none")
ggsave(filename="Fig2_AllSp_Isotopes.png", plot=Fig2, dpi=500, width=17, height=12, units="in")


##############################################################################
##### ------- Add Baseline Data, For Trophic Level Calculations: ------- #####
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

