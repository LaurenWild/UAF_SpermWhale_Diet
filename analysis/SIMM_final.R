# This code uses SIMMR and Mix-SIAR models to assess proportional contribution of groundfish and squid to sperm whale diet in the Gulf of Alaska
# Created by Lauren Wild, March 2019
# lawild@alaska.edu
# Packages used require "jags" (just another gibs sampler), downloaded from the internet. Search "jags for r" on google

#### NOTES: ####
#load('') is a better way (than read.table/.csv) to load your data in without needing Users/laurenwild/desktop...etc. 
#if you save the script within an R-project, that project automatically looks within that directory.
#load("../") this code allows you to go up a folder in where r looks for the file... 
### Line 78, SIMMR package analysis begins
### Line 328, MixSIAR package analysis begins

#Load libraries: 
library(ggplot2)
library(Hmisc)
library(psych)
library(devtools)
library(MASS)
library(car)
library(stats)
library(ggpubr)
library(dplyr)
library(viridis) # color scheme for plots that is easy to read (from simmr)
#install.packages('rjags')
library(rjags)
library(siar)
library(simmr)
library(MixSIAR)  ### Code starts line 853
library(patternplot) #to put patterns in the fill of boxplots, etc. 

### First, load data into R:
#Sperm whale (i.e. "mixtures") data:
PmSimmData<-read.table(here::here("data/PmSIMMData.csv"),sep=",",header=TRUE)
View(PmSimmData) 

PmSimmData<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSIMMData.csv',sep=",",header=TRUE) #All whales, full spreadsheet with all columns
View(PmSimmData)
PmSimmData2<-PmSimmData[,c(2,3,9,16,17,18)] #Pare down the data frame to just what we're interested in
PmSimmData3<-PmSimmData2[,1:2]
PmSimmData3<-as.matrix(PmSimmData3[,c(2,1)])
#Break up data-frame into subsets of whale groupings, then reduce so it is just the isotope ratios: 
PmSerialSimm<-subset(PmSimmData2, Frequent=="Frequent")
PmSerialSimm<-as.matrix(PmSerialSimm[,1:2])
PmNonSerialSimm<-subset(PmSimmData2, Frequent == "Non-Frequent")
PmNonSerialSimm<-as.matrix(PmNonSerialSimm[,1:2])
PmFreqNFreqSimm<-rbind(PmSerialSimm, PmNonSerialSimm) #data frame with just the freq & non-freq depr
PmRecentSimm<-subset(PmSimmData2, Recent=="Recent" )
PmRecentSimm<-as.matrix(PmRecentSimm[,1:2])
PmOldSimm<-subset(PmSimmData2, Recent=="Old")
PmOldSimm<-as.matrix(PmOldSimm[,1:2])
PmSeasonSimm<-PmSimmData2[,c(1,2,6)]
PmSeasonSimm$Grp[PmSeasonSimm$Season=='Early'] <- 1   
PmSeasonSimm$Grp[PmSeasonSimm$Season=='Mid'] <- 2
PmSeasonSimm$Grp[PmSeasonSimm$Season=='Late'] <- 3
PmSeasonSimm<-PmSeasonSimm[,c(1,2,4)]
PmSeasonSimm<-PmSeasonSimm[order(PmSeasonSimm$Grp),] #order by groups (to then separate into 2 matrices of just the isotope values and just the group number)
PmSeasonSimmData<-as.matrix(PmSeasonSimm[,1:2]) #Break into just consumer matrix (both in same order!)
PmSeasonSimmGrpData<-as.matrix(PmSeasonSimm[,3]) #and just group matrix (both in same order!)
names(PmSeasonSimmGrpData)[1]<-"Group"
grp<-as.integer(c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3))

#Separate the serial depredators into early, mid, and late summer samples:
PmSerialLate<-PmSerialSimm[c(1,2),]
PmSerialLate<-as.matrix(PmSerialFall[,1:2])
PmSerialEarly<-PmSerialSimm[c(3,4),]
PmSerialEarly<-as.matrix(PmSerialEarlySumm[,1:2])
PmSerialMid<-PmSerialSimm[c(5:10),]
PmSerialMid<-as.matrix(PmSerialSummer[,1:2])

### Load prey isotope data into R: 
Source <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources.csv', sep=",", header=TRUE)
SourceGp_AfSaIaCombined <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources-GpIa3.csv', sep=",", header=TRUE)
SourceGp_AfSaIaCombined <- SourceGp_AfSaIaCombined[,c(1,4,5,2,3,6)]
colnames(SourceGp_AfSaIaCombined)<-c("Species","d15N.m", "d15N.sd", "d13C.m", "d13C.sd")
SourceGp_AfSaIaCombined<-SourceGp_AfSaIaCombined[,c(1:5)]
SourceGp_AfSaIaCombined$Species <- c("Clubhook Squid","Glass Squid","Grenadier","Magister Squid","Shortraker Rockfish","Skate", "Sablefish.Dogfish")

### Load Source TEF data into R: 
SourceTEFsd<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd.csv', sep=",", header=TRUE)
SourceTEFsd_GpIa<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd-GpIa.csv', sep=",", header=TRUE)
SourceTEFsd_Gp_AfSaIa_Comb<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd-GpIa-AfIaSaComb.csv', sep=",", header=TRUE)
colnames(SourceTEFsd_Gp_AfSaIa_Comb)<-c("Species","d15N.m", "d15N.sd", "d13C.m", "d13C.sd")
SourceTEFsd_Gp_AfSaIa_Comb$Species <- c("Clubhook Squid","Glass Squid","Grenadier","Magister Squid","Shortraker Rockfish","Skate", "Sablefish.Dogfish")

####################################################################################
##### SIMMR package analysis - update from SIAR package #####
####################################################################################
# Load data into simmr, plot isospace, Run simmr_mcmc() model, run summary diagnostics, plot output (boxplot, matrix, histogram, density, convergence)

### SIMMR 1 - All whales non-grouped, all 9 prey, AfSaIa combined, TEF of 2.12. with s.d.
simmr_1_comb = simmr_load(mixtures=PmSimmData3,  
                          source_names=as.character(Source[,1]),
                          source_means=Source[,c(2,4)],
                          source_sds=Source[,c(3,5)],
                          correction_means=SourceTEFsd[,c(2,4)],
                          correction_sds=SourceTEFsd[,c(3,5)])
#Plot isospace of the data
setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig4_Isospace.tiff", height = 12, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
plot(simmr_1_comb, title=NULL, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
dev.off()

#Run the SIMMR model using mcmc
simmr_1_comb_out = simmr_mcmc(simmr_1_comb)
summary(simmr_1_comb_out, type='statistics')
plot(simmr_1_comb_out, type='boxplot', title=NULL)
plot(simmr_1_comb_out, type = 'matrix', title=NULL)
plot(simmr_1_comb_out, type = 'histogram', title=NULL)
plot(simmr_1_comb_out, type="density", title=NULL)
plot(simmr_1_comb_out, type="convergence", title=NULL)

### SIMMR 2 - Frequent depredating whales, all 9 prey, AfSaIa combined, TEF of 2.12. with s.d.
simmr_2_comb = simmr_load(mixtures=PmSerialSimm,
                          source_names=as.character(SourceGp_AfSaIaCombined[,1]),
                          source_means=SourceGp_AfSaIaCombined[,c(2,4)],
                          source_sds=SourceGp_AfSaIaCombined[,c(3,5)],
                          correction_means=SourceTEFsd_Gp_AfSaIa_Comb[,c(2,4)],
                          correction_sds=SourceTEFsd_Gp_AfSaIa_Comb[,c(3,5)])
#Plot isospace of the data
plot(simmr_2_comb, title='Isospace of frequent depredators') + coord_flip()

#Run the SIMMR model:
simmr_2_comb_out = simmr_mcmc(simmr_2_comb) #Serial Depredators (n=10)
summary(simmr_2_comb_out, type='diagnostics')
summary(simmr_1_comb_out, type='statistics')
summary(simmr_2_comb_out)
plot(simmr_2_comb_out, type = 'boxplot') 
plot(simmr_2_comb_out, type = 'matrix')
plot(simmr_2_comb_out, type = 'histogram') # or type = 'density', or type = 'convergence'

### SIMMR 3 - Non-frequent depredating whales, all 9 prey, AfSaIa combined, TEF of 2.12. with s.d.
simmr_3_comb = simmr_load(mixtures=PmNonSerialSimm,
                          source_names=as.character(SourceGp_AfSaIaCombined[,1]),
                          source_means=SourceGp_AfSaIaCombined[,c(2,4)],
                          source_sds=SourceGp_AfSaIaCombined[,c(3,5)],
                          correction_means=SourceTEFsd_Gp_AfSaIa_Comb[,c(2,4)],
                          correction_sds=SourceTEFsd_Gp_AfSaIa_Comb[,c(3,5)])
#Plot isospace of the data
plot(simmr_3_comb, title='Isospace of non-frequent depredators') + coord_flip()

#Run SIMMR model:
simmr_3_comb_out = simmr_mcmc(simmr_3_comb) #Non-Serial Depredators (n=20)
summary(simmr_3_comb_out, type='diagnostics')
summary(simmr_3_comb_out, type='statistics')
plot(simmr_3_comb_out, type = 'boxplot')  #non-serial depredators with sd
plot(simmr_3_comb_out, type = 'matrix')  # or type = 'histogram' / 'density' / or 'convergence'

### SIMMR 4 - Recent samples, all 9 prey, AfSaIa combined, TEF of 2.12. with s.d.
simmr_4_comb = simmr_load(mixtures=PmRecentSimm,
                          source_names=as.character(SourceGp_AfSaIaCombined[,1]),
                          source_means=SourceGp_AfSaIaCombined[,c(2,4)],
                          source_sds=SourceGp_AfSaIaCombined[,c(3,5)],
                          correction_means=SourceTEFsd_Gp_AfSaIa_Comb[,c(2,4)],
                          correction_sds=SourceTEFsd_Gp_AfSaIa_Comb[,c(3,5)])
#Plot isospace of the data
plot(simmr_4_comb, title='Isospace of recently sampled whales') + coord_flip()

#Run SIMMR model:
simmr_4_comb_out = simmr_mcmc(simmr_4_comb) #Recent samples  (n=23)
summary(simmr_4_comb_out, type='diagnostics')
summary(simmr_4_comb_out, type='statistics')
plot(simmr_4_comb_out, type = 'boxplot')  
plot(simmr_4_comb_out, type = 'matrix')  # or type = 'histogram' / 'density' / or 'convergence'

### SIMMR 5 - Old samples, all 9 prey, AfSaIa combined, TEF of 2.12. with s.d.
simmr_5_comb = simmr_load(mixtures=PmOldSimm,
                          source_names=as.character(SourceGp_AfSaIaCombined[,1]),
                          source_means=SourceGp_AfSaIaCombined[,c(2,4)],
                          source_sds=SourceGp_AfSaIaCombined[,c(3,5)],
                          correction_means=SourceTEFsd_Gp_AfSaIa_Comb[,c(2,4)],
                          correction_sds=SourceTEFsd_Gp_AfSaIa_Comb[,c(3,5)])
#Plot isospace of the data
plot(simmr_5_comb, title='Isospace of older whale samples') + coord_flip()

#Run SIMMR model:
simmr_5_comb_out = simmr_mcmc(simmr_5_comb) #Old samples  (n=10)
summary(simmr_5_comb_out, type='diagnostics')
summary(simmr_5_comb_out, type='statistics')
plot(simmr_5_comb_out, type = 'density')  
plot(simmr_5_comb_out, type = 'matrix')  # or type = 'histogram' / 'density' / or 'convergence'

### SIMMR 6 - Seasonal samples, all 9 prey, AfSaIa combined, TEF of 2.12. with s.d.
simmr_6_comb = simmr_load(mixtures=PmSeasonSimmData,
                          source_names=as.character(SourceGp_AfSaIaCombined[,1]),
                          source_means=SourceGp_AfSaIaCombined[,c(2,4)],
                          source_sds=SourceGp_AfSaIaCombined[,c(3,5)],
                          correction_means=SourceTEFsd_Gp_AfSaIa_Comb[,c(2,4)],
                          correction_sds=SourceTEFsd_Gp_AfSaIa_Comb[,c(3,5)],
                          group=as.integer(PmSeasonSimmGrpData))

#Plot isospace of the data
plot(simmr_6_comb, group=1:3,title='Isospace of seasonal whale samples') + coord_flip()

#Run SIMMR model:
simmr_6_comb_out = simmr_mcmc(simmr_6_comb) #Seasonal samples
summary(simmr_6_comb_out, type='diagnostics')
summary(simmr_6_comb_out, type='statistics', group=1:3)
plot(simmr_6_comb_out, type = 'boxplot',group=1:3)  
plot(simmr_6_comb_out, type = 'matrix')  # or type = 'histogram' / 'density' / or 'convergence'
compare_groups(simmr_6_comb_out, source='Sablefish.Dogfish.Ragfish', groups=1:3) #builds boxplots for each group
compare_groups(simmr_6_comb_out, source='Shortraker Rockfish', groups=1:3)
compare_groups(simmr_6_comb_out, source='Clubhook Squid', groups=1:3)
compare_groups(simmr_6_comb_out, source='Glass Squid', groups=1:3)
compare_groups(simmr_6_comb_out, source='Skate', groups=1:3)
compare_groups(simmr_6_comb_out, source='Grenadier', groups=1:3)
compare_groups(simmr_6_comb_out, source='Magister Squid', groups=1:3)

plot(simmr_6_comb_out,type='boxplot',title='Late season whales',group=3)
plot(simmr_6_comb_out,type='boxplot',title='Mid season whales',group=2)
plot(simmr_6_comb_out,type='boxplot',title='Early season whales',group=1)

###################################################################################
### This code builds a data frame with the mean for each prey, and a column with predator
.AllWhalesComb_out<-as.data.frame(as.matrix(simmr_1_comb_out$output[[1]]))
AllWhalesPropComb<- as.data.frame(as.matrix(colMeans(.AllWhales_out))) #should be same values as if you did summary(simmr_x_out, type='statistics)
write.table(AllWhalesPropComb, file="AllWhalesPropComb.csv", sep=",")

.FreqDepComb_out<-as.data.frame(as.matrix(simmr_2_comb_out$output[[1]]))
.NonFreqDepComb_out<-as.data.frame(as.matrix(simmr_3_comb_out$output[[1]]))
FreqDepPropComb<- as.data.frame(as.matrix(colMeans(.FreqDepComb_out)))
NonFreqDepPropComb<- as.data.frame(as.matrix(colMeans(.NonFreqDepComb_out)))
write.table(FreqDepPropComb, file="FreqDepPropComb.csv", sep=",")
write.table(NonFreqDepPropComb, file="NonFreqDepPropComb.csv", sep=",")

.OldComb_out<-as.data.frame(as.matrix(simmr_5_comb_out$output[[1]]))
.RecentComb_out<-as.data.frame(as.matrix(simmr_4_comb_out$output[[1]]))
OldPropComb<- as.data.frame(as.matrix(colMeans(.OldComb_out)))
RecentPropComb<- as.data.frame(as.matrix(colMeans(.RecentComb_out)))
write.table(OldPropComb, file="OldPropComb.csv", sep=",")
write.table(RecentPropComb, file="RecentPropComb.csv", sep=",")

##### Build proportion of prey to predator diet plots (first used for AMSS talk, Jan 2019)
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

setwd("~/Desktop/UAF/Thesis/StableIsotopes/Data/")
#read in prop diet data:
APMComb = read.csv("AllWhalesPropComb.csv")  #All 9 samples (AfSaIa Combined)
FPMComb = read.csv("FreqDepPropComb.csv")  #All 9 samples (AfSaIa Combined)
NFPMComb = read.csv("NonFreqDepPropComb.csv")  #All 9 samples (AfSaIa Combined)
OPMComb = read.csv("OldPropComb.csv")  #All 9 samples (AfSaIa Combined)
RPMComb = read.csv("RecentPropComb.csv")  #All 9 samples (AfSaIa Combined)
ESPMComb = read.csv("EarlyDepPropComb.csv") #All 9 samples (AfSaIa Combined)
MSPMComb = read.csv("MidDepPropComb.csv") #All 9 samples (AfSaIa Combined)
LSPMComb = read.csv("LateDepPropComb.csv") #All 9 samples (AfSaIa Combined)
#rename predators:
APMComb$Predator = "All Sperm Whales"
FPMComb$Predator = "Frequent Depredators"
NFPMComb$Predator = "Non-Frequent Depredators"
OPMComb$Predator = "Older Samples"
RPMComb$Predator = "Recent Samples"
ESPMComb$Predator = "Early Season"
MSPMComb$Predator = "Mid Season"
LSPMComb$Predator = "Late Season"
#combine all species in a dataframe:
FNFComb = rbind(FPMComb,NFPMComb)
ORPMComb = rbind(OPMComb,RPMComb)
SeasonPMComb= rbind(ESPMComb,MSPMComb,LSPMComb)
require(dplyr)
# Order prey items by phylogeny:
APMComb$Prey.Species = ordered(APMComb$Prey.Species, levels= c("Magister Squid","Glass Squid", "Clubhook Squid", "Skate", "Grenadier", "Shortraker Rockfish", "Sablefish.Dogfish.Ragfish"))
FNFComb$Prey.Species = ordered(FNFComb$Prey.Species, levels= c("Magister Squid","Glass Squid", "Clubhook Squid", "Skate", "Grenadier", "Shortraker Rockfish", "Sablefish.Dogfish.Ragfish"))
ORPMComb$Prey.Species = ordered(ORPMComb$Prey.Species, levels= c("Magister Squid","Glass Squid", "Clubhook Squid", "Skate", "Grenadier", "Shortraker Rockfish", "Sablefish.Dogfish.Ragfish"))
SeasonPMComb$Prey.Species = ordered(SeasonPMComb$Prey.Species, levels= c("Magister Squid","Glass Squid", "Clubhook Squid", "Skate", "Grenadier", "Shortraker Rockfish", "Sablefish.Dogfish.Ragfish"))
SeasonPMComb$Predator = ordered(SeasonPMComb$Predator, levels = c("Early Season", 'Mid Season', "Late Season"))
FNFComb$Predator = ordered(FNFComb$Predator, levels = c('Non-Frequent Depredators', 'Frequent Depredators'))

require(ggplot2)
AllWhalesPlot3 = ggplot(APMComb, aes(x = Predator, y = propPrey, fill = Prey.Species, order = -as.numeric(Prey.Species))) + 
  geom_bar(position = "fill", stat = "identity", col='black', width=0.25) +
  theme_bw() +
  ggtitle("") +
  theme(axis.text = element_text(family="Arial", size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=14))+
  scale_fill_manual(values = c("Glass Squid"="lightpink1", "Clubhook Squid"="firebrick", "Magister Squid"="mistyrose", "Skate"="mediumpurple1", "Grenadier"="lightcyan2", "Sablefish.Dogfish.Ragfish"="deepskyblue", "Shortraker Rockfish"="cadetblue1"), drop=FALSE) +
  labs(x="", y="Proportion by Weight") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) 
AllWhalesPlot3
ggsave(filename="AllWhaleComb.png", plot=AllWhalesPlot3, dpi=500, width=11.5, height=7, units="in")

library(stringr)
FNFComb$Predator2 = str_wrap(FNFComb$Predator, width = 10)
FNFComb$Predator3 = ordered(FNFComb$Predator2, levels = c('Non-Frequent Depredators', 'Frequent Depredators'))###Doesn't work?! Supposed to order NFN first, on the left, for presentation...
DepredatorsPlot3 = ggplot(FNFComb, aes(x = Predator2, y = propPrey, fill = Prey.Species, order = -as.numeric(Prey.Species))) + 
  geom_bar(position = "fill", stat = "identity", col='black', width=0.5) +
  theme_bw() +
  ggtitle("") +
  theme(axis.text = element_text(family="Arial", size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  theme(axis.title.x = element_text(vjust=1.1, size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=14))+
  scale_fill_manual(values = c("Glass Squid"="lightpink1", "Clubhook Squid"="firebrick", "Magister Squid"="mistyrose", "Skate"="mediumpurple1", "Grenadier"="lightcyan2", "Sablefish.Dogfish.Ragfish"="deepskyblue", "Shortraker Rockfish"="cadetblue1"), drop=FALSE) +
  labs(x="", y="Proportion by Weight") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) 
DepredatorsPlot3
ggsave(filename="FreqNonFreqComb.png", plot=DepredatorsPlot3, dpi=500, width=11.5, height=7, units="in")

OldRecentPlot3 = ggplot(ORPMComb, aes(x = Predator, y = propPrey, fill = Prey.Species, order = -as.numeric(Prey.Species), col=black)) + 
  geom_bar(position = "fill", stat = "identity", col="black", width=0.5) +
  theme_bw() +
  ggtitle("") +
  theme(axis.text = element_text(family="Arial", size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  theme(axis.title.x = element_text(vjust=1.1, size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=14))+
  scale_fill_manual(values = c("Glass Squid"="lightpink1", "Clubhook Squid"="firebrick", "Magister Squid"="mistyrose", "Skate"="mediumpurple1", "Grenadier"="lightcyan2", "Sablefish.Dogfish.Ragfish"="deepskyblue", "Shortraker Rockfish"="cadetblue1"), drop=FALSE) +
  labs(x="", y="Proportion by Weight") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) 
OldRecentPlot3
ggsave(filename="OldRecentComb.png", plot=OldRecentPlot3, dpi=500, width=11.5, height=7, units="in")

SeasonalPlot3 = ggplot(SeasonPMComb, aes(x = Predator, y = propPrey, fill = Prey.Species, order = -as.numeric(Prey.Species))) + 
  geom_bar(position = "fill", stat = "identity", col="black", width=0.6) +
  theme_bw() +
  ggtitle("") +
  theme(axis.text = element_text(family="Arial", size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  theme(axis.title.x = element_text(vjust=1.1, size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=14))+
  scale_fill_manual(values = c("Glass Squid"="lightpink1", "Clubhook Squid"="firebrick", "Magister Squid"="mistyrose", "Skate"="mediumpurple1", "Grenadier"="lightcyan2", "Sablefish.Dogfish.Ragfish"="deepskyblue", "Shortraker Rockfish"="cadetblue1"), drop=FALSE) +
  labs(x="", y="Proportion by Weight") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) 
SeasonalPlot3
ggsave(filename="SeasonalComb.png", plot=SeasonalPlot3, dpi=500, width=11.5, height=7, units="in")

##########################################################################################
##########################################################################################

#####################################################
########### COMPARE MODELS WITH MIX-SIAR ############
#####################################################
# This code was created by Lauren Wild to estimate the proportional contribution of prey items to sperm whale diet in the Gulf of Alaska.
# Original Authors: Brian Stock [cre, aut], Brice Semmens [aut], Eric Ward [ctb], Andrew Parnell [ctb], Andrew Jackson [ctb], Donald Phillips [ctb]
# Maintained by Brian Stock, b1stock@ucsd.edu
# Uses rjags and JAGS
# https://github.com/brianstock/MixSIAR
# cite: Stock, B.C. and B.X. Semmens (2016). MixSIAR GUI User Manual. Version 3.1. https://github.com/brianstock/MixSIAR/. doi: 10.5281/zenodo.47719
# https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf
browseVignettes("MixSIAR") #shows examples with scripts

mixsiar.dir <- find.package("MixSIAR")
paste0(mixsiar.dir,"/example_scripts")
#source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_killerwhale.R")) #Needs 'splancs' package to run calc_area #run sexample scripts for killer whale example! 
#setwd('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/')

#Load in data for MixSIAR:
#Full data set, no factors. 
mixPm <- load_mix_data(filename="data/PmSimmData2.csv",
                       iso_names=c("d13C","d15N"),
                       factors=NULL,
                       fac_random=NULL,
                       fac_nested=NULL,
                       cont_effects=NULL) #Needs to be in format of just d13C and d15N columns, and "factor" column if there is one. 
#Temporal data with column of "old" or "recent" sample
mixPmTemporal <- load_mix_data(filename="data/PmMixSIARDataTemporal.csv",
                               iso_names=c("d13C","d15N"),
                               factors='Temporal',
                               fac_random=TRUE,
                               fac_nested=FALSE,
                               cont_effects=NULL)
# "Seasonal" with Early, Mid, and Late season factor variable
mixPmSeasonal <- load_mix_data(filename="data/PmMixSIARDataSeasonal.csv",
                               iso_names=c("d13C","d15N"),
                               factors='Season',
                               fac_random=TRUE,
                               fac_nested=FALSE,
                               cont_effects=NULL)
mixPmFreqNFreq <- load_mix_data(filename="data/PmMixSIARDataSerialNSerial.csv",
                                iso_names=c("d13C","d15N"),
                                factors='Frequent',
                                fac_random=TRUE,
                                fac_nested=FALSE,
                                cont_effects=NULL)

#Load in source data for MixSIAR:
sourcePm.Comb <- load_source_data(filename="data/PmSources_SaAf_Comb.csv",
                                  source_factors=NULL,
                                  conc_dep=FALSE,
                                  data_type="means",
                                  mixPm)   #Source data frame needs to be in format: "Species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", "n" (sample size each source)
sourcePm.TempComb <- load_source_data(filename="data/PmSources_SaAf_Comb.csv",
                                      source_factors=NULL,
                                      conc_dep=FALSE,
                                      data_type="means",
                                      mixPmTemporal)
sourcePm.SeasonComb <- load_source_data(filename="data/PmSources_SaAf_Comb.csv",
                                        source_factors=NULL,
                                        conc_dep=FALSE,
                                        data_type="means",
                                        mixPmSeasonal)
sourcePm.FreqNFreqComb <- load_source_data(filename="data/PmSources_SaAf_Comb.csv",
                                           source_factors=NULL,
                                           conc_dep=FALSE,
                                           data_type="means",
                                           mixPmFreqNFreq)
discrPm.Comb <- load_discr_data(filename="data/PmSourceTEFsd-AfSaComb.csv", mixPm)
discrPm.TempComb <- load_discr_data(filename="data/PmSourceTEFsd-AfSaComb.csv", mixPmTemporal)
discrPm.SeasonComb <- load_discr_data(filename="data/PmSourceTEFsd-AfSaComb.csv", mixPmSeasonal)
discrPm.FreqNFreqComb <- load_discr_data(filename="data/PmSourceTEFsd-AfSaComb.csv", mixPmFreqNFreq)

# Make an isospace plot
isospace<-plot_data(filename="isospace_plot_combined", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPm,sourcePm.Comb,discrPm.Comb) 
coord_flip(isospace)
plot_data(filename="isospace_plot_combined_Temporal", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmTemporal,sourcePm.TempComb,discrPm.TempComb)
plot_data(filename="isospace_plot_combined_Seasonal", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmSeasonal,sourcePm.SeasonComb,discrPm.SeasonComb)
plot_data(filename="isospace_plot_combined_FreqNFreq", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmFreqNFreq,sourcePm.FreqNFreqComb,discrPm.FreqNFreqComb)

# Plot uninformative prior
plot_prior(alpha.prior=1, sourcePm.Comb, filename = "prior_plot_pm_comb_uninf", plot_save_pdf=FALSE)
plot_prior(alpha.prior=1, sourcePm.TempComb, filename = "prior_plot_pmTempComb_uninf", plot_save_pdf=TRUE)
plot_prior(alpha.prior=1, sourcePm.SeasonComb, filename = "prior_plot_pmSeasonComb_uninf", plot_save_pdf=FALSE)
plot_prior(alpha.prior=1, sourcePm.FreqNFreqComb, filename = "prior_plot_pmFnFComb_uninf", plot_save_pdf=FALSE)

# Define model structure and write JAGS model file
model_Pm.Comb <- "MixSIAR_model_Pm_Comb_uninf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_Pm.Comb, resid_err, process_err, mixPm, sourcePm.Comb)
model_Pm_TempComb <- "MixSIAR_model_Pm_TempComb_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_TempComb, resid_err, process_err, mixPmTemporal, sourcePm.TempComb)
model_Pm_SeasonComb <- "MixSIAR_model_Pm_SeasonComb_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_SeasonComb, resid_err, process_err, mixPmSeasonal, sourcePm.SeasonComb)
model_Pm_FreqNFreqComb <- "MixSIAR_model_Pm_FreqNFreqComb_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_FreqNFreqComb, resid_err, process_err, mixPmFreqNFreq, sourcePm.FreqNFreqComb)

# Run the JAGS model ("test" took ~20sec; "normal" should be good... about 35min, with chain length 100,000; burn in 50,000, and thin 50 with 3 chains and it seems to converge nicely) 
# Pg 14 of MixSIAR manual lists how many MCMC burn-ins and chains each one uses ("very long" is 1,000,000 chain length, 500,000 burn-in, and thin by 500, should take about 5hr to run)
jags.uninf.Pm.Comb <- run_model(run="normal",mixPm,sourcePm.Comb,discrPm.Comb,model_Pm.Comb,alpha.prior = 1, resid_err, process_err) #took 35 min to run.
jags.uninf.Pm.TempComb <- run_model(run="normal",mixPmTemporal,sourcePm.TempComb,discrPm.TempComb, model_Pm_TempComb,alpha.prior = 1, resid_err, process_err) #'normal' took 35min 
jags.uninf.Pm.SeasonComb <- run_model(run="normal",mixPmSeasonal,sourcePm.SeasonComb,discrPm.SeasonComb,model_Pm_SeasonComb,alpha.prior = 1, resid_err, process_err) #Took 35min
jags.uninf.Pm.FreqNFreqComb <- run_model(run="normal",mixPmFreqNFreq,sourcePm.FreqNFreqComb, discrPm.FreqNFreqComb,model_Pm_FreqNFreqComb,alpha.prior = 1, resid_err, process_err) #Took 35min to run

# get posterior medians for new source groupings
#apply(combined$post, 2, median)
#summary_stat(jags.uninf.Pm.Comb, meanSD=FALSE, quantiles=c(.025,.5,.975), savetxt=FALSE)

# Process diagnostics, summary stats, and posterior plots
jags_output_fullmodel<- output_JAGS(jags.uninf.Pm.Comb, mixPm, sourcePm.Comb) #DIC=119.4686
jags_output_tempcomb<- output_JAGS(jags.uninf.Pm.TempComb, mixPmTemporal, sourcePm.TempComb) #DIC 118.4287 - model doesn't converge super well, might need to run more iterations! 
jags_output_seasoncomb<- output_JAGS(jags.uninf.Pm.SeasonComb, mixPmSeasonal, sourcePm.SeasonComb) #DIC 120.8179 - model did converge pretty well,
jags_output_fnfcomb<- output_JAGS(jags.uninf.Pm.FreqNFreqComb, mixPmFreqNFreq, sourcePm.FreqNFreqComb) #DIC 114.8703 - model converged pretty well

#Using saveRDS to save jags model output - isn't working, but leaving code here:
saveRDS(output_JAGS(jags.uninf.Pm.TempComb, mixPmTemporal, sourcePm.TempComb), file = "jags_output_tempcomb.rds") #Saves the jags model output! 
saveRDS(jags_output_seasoncomb, file = "jags_output_seasoncomb.rds") #Saves the jags model output!
readRDS(file="jags_output_tempcomb.rds")

### Use this code to access posterior density outputs from the models:
jags.uninf.Pm.TempComb$BUGSoutput #Look at the model output to isolate posterior densities... 
jags.uninf.Pm.TempComb$BUGSoutput$sims.list$p.global[,1] #gives posterior output for first species
jags.uninf.Pm.TempComb$BUGSoutput$mean$p.global #gives means for each species
jags.uninf.Pm.TempComb$BUGSoutput$mean$p.fac1 #gives mean for each species and each whale group
jags.uninf.Pm.TempComb$BUGSoutput$sims.list$p.fac1[,,5] #gives posterior for each whale group of just my 5th prey item (sablefish group)

# Find out which group is "Frequent" depredators:
freq_level = which(levels(mixPmFreqNFreq$data$Frequent) == "Frequent") #first column
# Find out which group is "Non-Frequent" depredators:
nonfreq_level = which(levels(mixPmFreqNFreq$data$Frequent) == "Non-Frequent") #2nd column
# Find out which source is "Sablefish.Dogfish"
source_1 = which(sourcePm.FreqNFreqComb$source_names == "Sablefish.Dogfish")
# Find out which source is "Skates"
source_2 = which(sourcePm.FreqNFreqComb$source_names == "Skate")

which(levels(mixPmTemporal$data$Temporal)=="Old") # [1] 1
which(levels(mixPmTemporal$data$Temporal)=="Recent") # [1] 2
which(sourcePm.TempComb$source_names == "Clubhook Squid") # [1] 1
which(sourcePm.TempComb$source_names == "Grenadier") # [1] 2
which(sourcePm.TempComb$source_names == "Magister Squid") # [1] 3
which(sourcePm.TempComb$source_names == "Sablefish.Dogfish") #[1] 4
which(sourcePm.TempComb$source_names == "Shortraker Rockfish") #[1] 5
which(sourcePm.TempComb$source_names == "Skate") #[1] 6
which(sourcePm.SeasonComb$source_names == "Clubhook Squid") # [1] 1
which(sourcePm.SeasonComb$source_names == "Grenadier") # [1] 2
which(sourcePm.SeasonComb$source_names == "Magister Squid") # [1] 3
which(sourcePm.SeasonComb$source_names == "Sablefish.Dogfish") #[1] 4
which(sourcePm.SeasonComb$source_names == "Shortraker Rockfish") #[1] 5
which(sourcePm.SeasonComb$source_names == "Skate") #[1] 6

#Prey List: 1=Clubhook squid, 2=Grenadier, 3=Magister Squid, 4=Sablefish Group, 5=Shortraker rockfish, 6=Skate
#Predator Factor List: TempComb (1=Old, 2=Recent); SeasonComb (1=Early, 2=Late, 3=Mid); FreqNFreqComb (1=Freq, 2=Non-Freq, 3=Unk)
AllCombSummary<-jags.uninf.Pm.Comb$BUGSoutput$summary
write.table(AllCombSummary, file="AllCombSum.csv",sep=',')
TempCombSummary<-jags.uninf.Pm.TempComb$BUGSoutput$summary
TempCombSummary2<-TempCombSummary[c("p.fac1[1,1]","p.fac1[2,1]","p.fac1[1,2]","p.fac1[2,2]","p.fac1[1,3]","p.fac1[2,3]","p.fac1[1,4]","p.fac1[2,4]","p.fac1[1,5]","p.fac1[2,5]","p.fac1[1,6]","p.fac1[2,6]"),]
write.table(TempCombSummary2, file="TempCombSum.csv", sep=',')
SeasonCombSummary<-jags.uninf.Pm.SeasonComb$BUGSoutput$summary
SeasonCombSummary2<-SeasonCombSummary[c("p.fac1[1,1]","p.fac1[2,1]","p.fac1[3,1]","p.fac1[1,2]","p.fac1[2,2]","p.fac1[3,2]","p.fac1[1,3]","p.fac1[2,3]","p.fac1[3,3]","p.fac1[1,4]","p.fac1[2,4]","p.fac1[3,4]","p.fac1[1,5]","p.fac1[2,5]","p.fac1[3,5]","p.fac1[1,6]","p.fac1[2,6]","p.fac1[3,6]"),]
write.table(SeasonCombSummary2, file="SeasonCombSum2.csv", sep=',')
FNFCombSummary<-jags.uninf.Pm.FreqNFreqComb$BUGSoutput$summary
FNFCombSummary2<-FNFCombSummary[c("p.fac1[1,1]","p.fac1[2,1]","p.fac1[3,1]","p.fac1[1,2]","p.fac1[2,2]","p.fac1[3,2]","p.fac1[1,3]","p.fac1[2,3]","p.fac1[3,3]","p.fac1[1,4]","p.fac1[2,4]","p.fac1[3,4]","p.fac1[1,5]","p.fac1[2,5]","p.fac1[3,5]","p.fac1[1,6]","p.fac1[2,6]","p.fac1[3,6]"),]
write.table(FNFCombSummary2, file="FNFCombSum2.csv", sep=',')

jags.uninf.Pm.TempComb$BUGSoutput$median$p.fac1 #gives just median, but shouldn't need, as it is also 50% in the full output summary above... 
jags.uninf.Pm.FreqNFreqComb$BUGSoutput$median$p.fac1
jags.uninf.Pm.SeasonComb$BUGSoutput$median$p.fac1

AllCombSumBox<-read.table('data/AllCombSumBox.csv', sep=",",header=TRUE)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')
AllComb<-ggplot(AllCombSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity", fill="grey", color="black") +
  xlab("Species") +
  ylab("Proportion of Diet")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_bw()+
  theme(axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), axis.title.x=element_text(size=28), axis.title.y=element_text(size=28))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
AllComb
ggsave(filename="All_Comb_Boxplot.png", plot=AllComb, dpi=500, width=17, height=12, units="in")

TempCombSumBox<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/TempCombSumBox3.csv',sep=",",header=TRUE)
TempCombSumBox$Middle<-as.numeric(TempCombSumBox$Middle)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')
A<-ggplot(TempCombSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax, color=Group, fill=Group), stat="identity", color="black") +
  scale_fill_manual(values=c("grey",'white'))+
  ylab("Proportion of Diet")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=14))+
  theme(legend.justification=c(0.8,0), legend.position=c(0.18,0.6))+
  theme(legend.box.background=element_rect(color="black"))+
  theme(legend.title=element_blank(), legend.text=element_text(size=11))+
  theme(legend.key.size=unit(0.8,'cm'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

A

SeasonCombSumBox<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/SeasonCombSumBox2.csv',sep=",",header=TRUE)
SeasonCombSumBox$Group<-factor(SeasonCombSumBox$Group, levels=c("Early","Mid","Late"), ordered=TRUE)
B<-ggplot(SeasonCombSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax, color=Group, fill=Group), stat="identity", color="black") +
  scale_fill_manual(values=c("grey",'white',"gray47"))+
  ylab("Proportion of Diet")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=14))+
  theme(legend.justification=c(0.8,0), legend.position=c(0.15,0.5))+
  theme(legend.box.background=element_rect(color="black"))+
  theme(legend.title=element_blank(), legend.text=element_text(size=11))+
  theme(legend.key.size=unit(0.8,'cm'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
B

FNFCombSumBox<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/FNFCombSumBox2.csv',sep=",",header=TRUE)
C<- ggplot(FNFCombSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax, color=Group,fill=Group), stat="identity", color="black") +
  scale_fill_manual(values=c("grey",'white'))+
  xlab("Species") +
  ylab("Proportion of Diet")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=12),axis.text.y=element_text(size=12), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))+
  theme(legend.justification=c(0.8,0), legend.position=c(0.24,0.6))+
  theme(legend.box.background=element_rect(color="black"))+
  theme(legend.title=element_blank(), legend.text=element_text(size=11))+
  theme(legend.key.size=unit(0.8,'cm'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
C

#library(egg)
#library(gridExtra)
#library(reshape)

Fig3<-ggarrange(A,B,C, labels = c("A", "B", "C"), ncol=1, nrow=3)
#setwd('/Users/laurenwild/Desktop')
tiff(filename="Fig3_boxplots.tiff", height = 24, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
Fig3
dev.off()

#ggsave(filename="Fig3_MixingModels.png", plot=Fig3, dpi=500, width=16, height=32, units="in")



###################################################################
# Compare Proportion of "Sablefish.Dogfish.Ragfish" vs "Skate"  
# This computes the probability that the proportion of "Skates" of "frequent deperedator" diets is bigger than the dietary proportion of "Sablefish.group" for "freq depr", given the data

prop_Sablefish.Dogfish.Ragfish_freq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,freq_level,source_1]
prop_Skate_freq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,freq_level,source_2]
prop_Sablefish.Dogfish.Ragfish_nfreq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,nonfreq_level,source_1]
prop_Skate_nfreq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,nonfreq_level,source_2]
prop_Clubhook_freq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,freq_level,1]
jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,nonfreq_level,1]

# Create a plot of the differences
diff = prop_Skate_freq - prop_Sablefish.Dogfish.Ragfish_freq
qplot(diff, geom = 'histogram',
      xlab = 'Posterior difference in dietary proportions between Skate and Sablefish/Dogfish for Frequent depredators')
diff2 = prop_Sablefish.Dogfish.Ragfish_nfreq - prop_Skate_nfreq
qplot(diff2, geom = 'histogram',
      xlab = 'Posterior difference in dietary proportions between Skate and Sablefish/Dogfish for non-frequent depredators')

diff3 = prop_Skate_freq - prop_Skate_nfreq

# Find the probability one is bigger than the other

sum(diff2>0)/length(diff2) # 0.822
#In the non-frequent depredator group, there is an 82% chance that whales are eating a higher proportion of sablefish than skates

sum(diff>0)/length(diff)
#In the frequent depredator group, there is an 80% chance that whales are eating a higher proportion of skates than sablefish

sum(diff3>0)/length(diff3) #0.87
#87% chance that the proportion of skates in freq depr diets is greater than than non-freq depredator diets, given the data. 

#Another way to compute it, this gives the proportion of time that frequent depredators eat more skates than non-freq
skatefnfcomb<-as.data.frame(as.matrix(jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,,7]))
sablefishfnfcomb<-as.data.frame(as.matrix(jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,,5]))

mean(skatefnfcomb$V1)
mean(skatefnfcomb$V2)
summary(skatefnfcomb$V1>skatefnfcomb$V2) #V1=freq, > V2(non-freq) 87% of the time... 
# So 2611 times out of 3000, the V1 column is greater than V2 column, or 87% of the time the V1 (freq depr) group of whales eats more skates than the V2 (non-freq) group of whales.

mean(sablefishfnfcomb$V1)
mean(sablefishfnfcomb$V2)
summary(sablefishfnfcomb$V1>sablefishfnfcomb$V2) #84% of the time... the V2 (non-freq) group of whales eats more sablefish than the V1 (freq) group of whales

#Look at old vs. recent model: 
sablefishtempcomb<-as.data.frame(as.matrix(jags.uninf.Pm.TempComb$BUGSoutput$sims.list$p.fac1[,,5]))
summary(sablefishtempcomb$V1>sablefishtempcomb$V2)
# So 2573 times out of 3000, the V1(recent) column is greater than V2(old) column, or 86% of the time the recent group of whales eats more sablefish than the old group of whales.

skatetempcomb<-as.data.frame(as.matrix(jags.uninf.Pm.TempComb$BUGSoutput$sims.list$p.fac1[,,7]))
summary(skatetempcomb$V1>skatetempcomb$V2)
sablefishtempcomb<-as.data.frame(as.matrix(jags.uninf.Pm.TempComb$BUGSoutput$sims.list$p.fac1[,,5]))
summary(sablefishtempcomb$V1>sablefishtempcomb$V2)


compare_models(x=list(mod.uninf = jags.uninf.Pm.Comb, mod.uninf.temporalcomb = jags.uninf.Pm.TempComb), loo=TRUE) #Temp model doesn't help explain data better, higher LOOic value
compare_models(x=list(mod.uninf = jags.uninf.Pm.Comb, mod.uninf.seasonalcomb = jags.uninf.Pm.SeasonComb), loo=TRUE) #Season model doesn't help explain data better, higher LOOic value
compare_models(x=list(mod.uninf = jags.uninf.Pm.Comb, mod.unif.FreqNFreqcomb = jags.uninf.Pm.FreqNFreqComb), loo=TRUE) #Freq dep model DOES help explain data better, lower LOOic value
