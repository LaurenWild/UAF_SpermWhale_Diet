#### Chapter 2 SIAR models to assess proportional contribution to Diet
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
library(dplyr)
library(doBy)
#library(plyr) #plyr probably gives errors, has been replaced by dplyr, but it's different
library(viridis) # color scheme for plots that is easy to read (from simmr)

### Package MixSIAR begins line 853 ####

### Package SIMMR; Need rjags, but first have to download jags from the internet. search "jags" for r on google
#install.packages('rjags')
library(rjags)
library(siar)
library(simmr)
library(MixSIAR)  ### Code starts line 853

#Load data into r
PmSimmData<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSIMMData.csv',sep=",",header=TRUE)
View(PmSimmData)
#load('') is a better way to load your data in without needing Users/laurenwild/desktop...etc. 
#if you save the script within an R-project, that project automatically looks within that directory. 
#load("../") this code allows you to go up a folder in where r looks for the file... 
Source <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources.csv', sep=",", header=TRUE)
View(Source)
SourceGp <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources-Gp.csv', sep=",", header=TRUE)
SourceGpIa <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources-GpIa.csv', sep=",", header=TRUE)
SourceGpIaDg <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources-GpIaDg.csv', sep=",", header=TRUE)
SourceGpDg_IaSaCombined <- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSources-GpDg_IaSaCombined.csv', sep=",", header=TRUE)
View(SourceGpDg_IaSaCombined)
SourceTEFse<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFse.csv', sep=",", header=TRUE)
View(SourceTEFse)
SourceTEFsd<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd.csv', sep=",", header=TRUE)
View(SourceTEFsd)
SourceTEFsd.GpIaDg<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd-GpIaDg.csv', sep=",", header=TRUE)
SourceTEFsd.Gp<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd-Gp.csv', sep=",", header=TRUE)
SourceTEFsd.GpIa<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd-GpIa.csv', sep=",", header=TRUE)
SourceTEF2.4.GpIa<- read.csv ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEF2.4-GpIa.csv', sep=",", header=TRUE)
SourceTEFsd.GpDg.IaSaCombined<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceTEFsd-GpDg-IaSaCombined.csv', sep=",", header=TRUE)
View(SourceTEFsd.GpDg.IaSaCombined)
SourceConcDepen<- read.table ('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSourceConcDepen.csv', sep=",", header=TRUE)
View(SourceConcDepen)

#Other Groupings of Sperm Whales
PmSerialSimm<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSISerialData.csv',sep=",",header=TRUE)
View(PmSerialSimm)
PmSerialSimm$Code<-1
PmSerialSimm<-PmSerialSimm[,c(2,3,23)]
PmSerialSimm<-PmSerialSimm[,c(3,1,2)]
View(PmSerialSimm)
PmNonSerialSimm<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSINonSerial.csv',sep=",",header=TRUE)
PmNonSerialSimm$Code<-2
PmNonSerialSimm<-PmNonSerialSimm[,c(2,3,23)]
PmNonSerialSimm<-PmNonSerialSimm[,c(3,1,2)]
#PmNonSerialSimm<-PmNonSerialSimm[,c(1,2,4)]
View(PmNonSerialSimm)
PmFreqSimm<-rbind(PmSerialSimm, PmNonSerialSimm)
View(PmFreqSimm)
PmRecentSimm<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSIRecent.csv',sep=",",header=TRUE)
PmRecentSimm$Code<-1
PmRecentSimm<-PmRecentSimm[,c(3,4,27)]
PmRecentSimm<-PmRecentSimm[,c(3,1,2)]
View(PmRecentSimm)
PmOldSimm<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSIOld.csv',sep=",",header=TRUE)
PmOldSimm$Code<-2
PmOldSimm<-PmOldSimm[,c(3,4,27)]
PmOldSimm<-PmOldSimm[,c(3,1,2)]
View(PmOldSimm)
PmROSimm<-rbind(PmRecentSimm, PmOldSimm)
View(PmROSimm)
PmSeasonSimm<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/PmSISeason.csv',sep=",",header=TRUE)
PmSeasonSimm$Grp[PmSeasonSimm$Season=='Early'] <- 1   
PmSeasonSimm$Grp[PmSeasonSimm$Season=='Mid'] <- 2
PmSeasonSimm$Grp[PmSeasonSimm$Season=='Late'] <- 3
PmSeasonSimm<-PmSeasonSimm[,c(3,4,24)]
PmSeasonSimm<-PmSeasonSimm[,c(3,1,2)]
PmSeasonSimm<-PmSeasonSimm[order(PmSeasonSimm$Grp),] #order by groups to separate into 2 matrices
View(PmSeasonSimm)
PmSeasonSimmData<-PmSeasonSimm[,2:3] #Break into just consumer matrix (both in same order!)
PmSeasonSimmGrpData<-PmSeasonSimm[,1] #and just group matrix (both in same order!)
grp<-as.integer(c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,
                  2,2,2,2,2,2,2,2,3,3,3,3,3,3,
                  3,3,3,3,3))

View(PmSerialSimm)
PmSerialFall<-PmSerialSimm[c(1,2),]
View(PmSerialFall)
PmSerialEarlySumm<-PmSerialSimm[c(3,4),]
PmSerialSummer<-PmSerialSimm[c(5:10),]

# Make each whale dataframe a matrix, and select only the d13C and d15N columns:
PmSimmData2<-as.matrix(PmSimmData[,2:3])
PmSerialSimmData<-as.matrix(PmSerialSimm[,2:3])
PmSerialFallData<-as.matrix(PmSerialFall[,2:3])
PmSerialEarlySummData<-as.matrix(PmSerialEarlySumm[,2:3])
PmSerialSummerData<-as.matrix(PmSerialSummer[,2:3])
PmNonSerialSimmData<-as.matrix(PmNonSerialSimm[,2:3])
PmRecentSimmData<-as.matrix(PmRecentSimm[,2:3])
PmOldSimmData<-as.matrix(PmOldSimm[,1:2])
PmSeasonSimmData<-PmSeasonSimm[,2:3] #Break into just consumer matrix (both in same order!)
PmSeasonSimmGrp<-PmSeasonSimm[,1] #and just group matrix (both in same order!)
names(PmSeasonSimmGrp)[1]<-"Group"
View(PmSeasonSimmGrp)
PmFreqSimm<-as.matrix(PmFreqSimm)
PmROSimm<-as.matrix(PmROSimm)
#SourceGpIaDg<-as.matrix(SourceGpIaDg)
#PmSourceTEFsd.GpIaDg<-as.matrix(PmSourceTEFsd.GpIaDg)

#CHECK HOW ISOSPACE ENCOMPASSES PREY:
NewSources<-SourceGpIa
NewSources$d15N.m<-NewSources$d15N.m+2.12
NewSources$d13C.m<-NewSources$d13C.m+0.95
NewSources.Pm<-rbind(NewSources, PmInnerAvg) #Need PmInnerAvg from SI_Pm_CH2.R, line 1347...
ggplot(NewSources.Pm,aes(d13C.m,d15N.m, label=Species, color=Species)) + geom_point(size=5) +
  geom_errorbarh(aes(xmax=NewSources.Pm$d13C.m+NewSources.Pm$d13C.sd,xmin=NewSources.Pm$d13C.m-NewSources.Pm$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=NewSources.Pm$d15N.m+NewSources.Pm$d15N.sd,ymin=NewSources.Pm$d15N.m-NewSources.Pm$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 4)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  scale_color_viridis_d()+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")
NewSources2<-SourceGpIa
NewSources2$d15N.m<-NewSources2$d15N.m+2.4
NewSources2$d13C.m<-NewSources2$d13C.m+1.01
NewSources2.Pm<-rbind(NewSources2, PmInnerAvg)

### From here go to SIMMR if want that... 

#FROM ANDREW PARNELL:
# SIAR query for Lauren Wild 17/5/18
## Used with SIAR, Now using SIMMR

# Clear the workspace
rm(list = ls())

# Load in packages
library(siar)

#In SIAR, Build the models:
?siarmcmcdirichletv4
model1<-siarmcmcdirichletv4(PmSimmData, Source, SourceTEF)
modelSerial<-siarmcmcdirichletv4(PmSerialSimmData, Source, SourceTEF)
modelNonSerial<-siarmcmcdirichletv4(PmNonSerialSimmData, Source, SourceTEF)
modelRecent<-siarmcmcdirichletv4(PmRecentSimmData, Source, SourceTEF)
modelOld<-siarmcmcdirichletv4(PmOldSimmData, Source, SourceTEF)

siarplotdata(model1)
siarplotdata(modelSerial)
siarplotdata(modelNonSerial)
siarplotdata(modelRecent)
siarplotdata(modelOld)

siarmatrixplot(model1)

siarhistograms(model1)
####################################################################################
# SIMMR version - much simpler (and more beautiful) ---------------------
####################################################################################
# Load data into simmr, plot isospace, Run simmr_mcmc() model, run summary diagnostics, plot output (boxplot, matrix, histogram, density, convergence)

### SIMMR 1 - All Whales ###
# All whales non-grouped, Standard 7 prey, TEF of 2.12 with more variable sd
simmr_1 = simmr_load(mixtures=PmSimmData2,
                     source_names=as.character(Source[,1]),
                     source_means=Source[,c(2,4)],
                     source_sds=Source[,c(3,5)],
                     correction_means=SourceTEFsd[,c(2,4)],
                     correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_1, title='Isospace of all whales (inner layer)') + coord_flip()
simmr_1_out = simmr_mcmc(simmr_1)  
summary(simmr_1_out, type='diagnostics')
simmr_combine1_out<-combine_sources(simmr_1_out, 
                                    to_combine=c('Sablefish', 'Spiny Dogfish'),
                                    new_source_name='Sablefish+SpinyDogfish')
plot(simmr_1_out, type = 'boxplot')
plot(simmr_1_out, type = 'matrix')
plot(simmr_1_out, type = 'histogram')
plot(simmr_1_out, type="density")
plot(simmr_1_out, type="convergence")
plot(simmr_combine1_out,type='boxplot',title='All Whales, combined sablefish & spiny dogfish')
plot(simmr_combine1_out, type='matrix')

#ALl whales non-grouped, standard 7 prey, TEF of 2.12 with smaller less variabl s.e.
simmr_1.5 = simmr_load(mixtures=PmSimmData2,
                       source_names=as.character(Source[,1]),
                       source_means=Source[,c(2,4)],
                       source_sds=Source[,c(3,5)],
                       correction_means=SourceTEFse[,c(2,4)],
                       correction_sds=SourceTEFse[,c(3,5)])
plot(simmr_1.5) + coord_flip()
simmr_1.5_out = simmr_mcmc(simmr_1.5)
summary(simmr_1.5_out, type='diagnostics')
plot(simmr_1.5_out, type = 'boxplot')

#All whales non-grouped, All 7 prey with added Glass Squid only, TEF of 2.12 with s.d.
simmr_1.Gp = simmr_load(mixtures=PmSimmData2,  
                        source_names=as.character(SourceGp[,1]),
                        source_means=SourceGp[,c(2,4)],
                        source_sds=SourceGp[,c(3,5)],
                        correction_means=SourceTEFsd.Gp[,c(2,4)],
                        correction_sds=SourceTEFsd.Gp[,c(3,5)])
simmr_1_Gp_out = simmr_mcmc(simmr_1.Gp)
simmr_combine1_Gp_out<-combine_sources(simmr_1_Gp_out, 
                                       to_combine=c('Sablefish', 'Spiny Dogfish'),
                                       new_source_name='Sablefish+SpinyDogfish')
plot(simmr_combine1_Gp_out$input, title="Isospace of all whales, combine Sblfsh & Dgfsh") + coord_flip()

#All whales non-grouped, all 9 prey, TEF of 2.12. with s.d.
simmr_1.GpIa = simmr_load(mixtures=PmSimmData2,  
                          source_names=as.character(SourceGpIa[,1]),
                          source_means=SourceGpIa[,c(2,4)],
                          source_sds=SourceGpIa[,c(3,5)],
                          correction_means=SourceTEFsd.GpIa[,c(2,4)],
                          correction_sds=SourceTEFsd.GpIa[,c(3,5)])
simmr_1_GpIa_out = simmr_mcmc(simmr_1.GpIa)
summary(simmr_1_GpIa_out, type='statistics')
plot(simmr_1_GpIa_out, type='boxplot', title=NULL)

#All whales non-grouped, All 9 main prey; TEF of 2.4 with s.d.
simmr_1.GpIa.TEF2.4 = simmr_load(mixtures=PmSimmData2,  
                                 source_names=as.character(SourceGpIa[,1]),
                                 source_means=SourceGpIa[,c(2,4)],
                                 source_sds=SourceGpIa[,c(3,5)],
                                 correction_means=SourceTEF2.4.GpIa[,c(2,4)],
                                 correction_sds=SourceTEF2.4.GpIa[,c(3,5)])
simmr_1_GpIa_TEF2.4_out = simmr_mcmc(simmr_1.GpIa.TEF2.4)

#All whales non-grouped, All 9 prey with Humboldt Squid added in too; TEF of 2.12 with s.d.
simmr_1.GpIaDg = simmr_load(mixtures=PmSimmData2,  
                            source_names=as.character(SourceGpIaDg[,1]),
                            source_means=SourceGpIaDg[,c(2,4)],
                            source_sds=SourceGpIaDg[,c(3,5)],
                            correction_means=SourceTEFsd.GpIaDg[,c(2,4)],
                            correction_sds=SourceTEFsd.GpIaDg[,c(3,5)])
plot(simmr_1.GpIaDg, title='Isospace with all whales & all prey') +coord_flip()
simmr_1_GpIaDg_out = simmr_mcmc(simmr_1.GpIaDg)
summary(simmr_1_GpIaDg_out, type='diagnostics')
simmr_combine1GpIaDg_out<-combine_sources(simmr_1_GpIaDg_out, 
                                          to_combine=c('Sablefish', 'Spiny Dogfish'),
                                          new_source_name='Sablefish+SpinyDogfish')
### SIMRR 2 - Serial Depredators ###
#Serial depredating whales, Main 7 Prey, TEF of 2.12 with s.d.
simmr_2 = simmr_load(mixtures=PmSerialSimmData,
                     source_names=as.character(Source[,1]),
                     source_means=Source[,c(2,4)],
                     source_sds=Source[,c(3,5)],
                     correction_means=SourceTEFsd[,c(2,4)],
                     correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_2, title='Isospace of serial depredators - sd') +coord_flip()
simmr_2_out = simmr_mcmc(simmr_2) #Serial Depredators (n=10)
summary(simmr_2_out, type='diagnostics')
summary(simmr_2_out)
simmr_combine2_out<-combine_sources(simmr_2_out, 
                                    to_combine=c('Sablefish', 'Spiny Dogfish'),
                                    new_source_name='Sablefish+SpinyDogfish')
plot(simmr_combine2_out$input, title="Isospace of Serial Depredators, combine Sblfsh & Dgfsh") + coord_flip()
plot(simmr_2_out, type = 'boxplot') 
plot(simmr_2_out, type = 'matrix')
plot(simmr_2_out, type = 'histogram') # or type = 'density', or type = 'convergence'
plot(simmr_combine2_out,type='boxplot',title='Serial Depredators, combined sablefish & spiny dogfish')

#Serial Depr from Early Summer (n=2), Main 7 prey, TEF of 2.12 and s.d.
simmr_2.2 = simmr_load(mixtures=PmSerialEarlySummData,
                       source_names=as.character(Source[,1]),
                       source_means=Source[,c(2,4)],
                       source_sds=Source[,c(3,5)],
                       correction_means=SourceTEFsd[,c(2,4)],
                       correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_2.2, title='Isospace of Early Summmer serial depredators')+coord_flip()
simmr_2.2_out = simmr_mcmc(simmr_2.2) #early season serial depredators (n=2)
summary(simmr_2.5_out, type='diagnostics')
plot(simmr_2.2_out, type = 'boxplot', title='Early summer serial depredators') #Early summer serial depredators 

#Serial Depr from Mid-Summer (n=6), Main 7 prey, TEF of 2.12 and s.d.
simmr_2.3 = simmr_load(mixtures=PmSerialSummerData,
                       source_names=as.character(Source[,1]),
                       source_means=Source[,c(2,4)],
                       source_sds=Source[,c(3,5)],
                       correction_means=SourceTEFsd[,c(2,4)],
                       correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_2.3, title='Isospace of Mid Summer serial depredators')+coord_flip()
simmr_2.3_out = simmr_mcmc(simmr_2.3) #summer serial depredators (n=6)
plot(simmr_2.3_out, type = 'boxplot', title='Summer serial depredators') #Summer serial depredators - looks like these guys really drive the skate preference

#Serial Depr from Late-Summer (n=2), Main 7 prey, TEF of 2.12 and s.d.
simmr_2.4 = simmr_load(mixtures=PmSerialFallData,
                       source_names=as.character(Source[,1]),
                       source_means=Source[,c(2,4)],
                       source_sds=Source[,c(3,5)],
                       correction_means=SourceTEFsd[,c(2,4)],
                       correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_2.4, title='Isospace of Late Summer serial depredators')+coord_flip()
simmr_2.4_out = simmr_mcmc(simmr_2.4) # Late Summer serial depredators (n=2)
plot(simmr_2.4_out, type = 'boxplot', title='Late Summer serial depredators') 

#Serial Depredators, Main 7 prey with Glass Squid, TEF of 2.12 and s.d.
simmr_2_Gp = simmr_load(mixtures=PmSerialSimmData,
                        source_names=as.character(SourceGp[,1]),
                        source_means=SourceGp[,c(2,4)],
                        source_sds=SourceGp[,c(3,5)],
                        correction_means=SourceTEFsd.Gp[,c(2,4)],
                        correction_sds=SourceTEFsd.Gp[,c(3,5)])
simmr_2_Gp_out=simmr_mcmc(simmr_2_Gp)
simmr_combine2_Gp_out<-combine_sources(simmr_2_Gp_out, 
                                       to_combine=c('Sablefish', 'Spiny Dogfish'),
                                       new_source_name='Sablefish+SpinyDogfish')

#Serial Depredators, All 9 prey, TEF of 2.12 and s.d.
simmr_2_GpIa = simmr_load(mixtures=PmSerialSimmData,
                          source_names=as.character(SourceGpIa[,1]),
                          source_means=SourceGpIa[,c(2,4)],
                          source_sds=SourceGpIa[,c(3,5)],
                          correction_means=SourceTEFsd.GpIa[,c(2,4)],
                          correction_sds=SourceTEFsd.GpIa[,c(3,5)])
simmr_2_GpIa_out=simmr_mcmc(simmr_2_GpIa)
summary(simmr_2_GpIa_out, type='diagnostics')
plot(simmr_2_GpIa_out, type='boxplot', title=NULL)

# Serial Depredators, All 9 prey + Humboldt Squid, TEF of 2.12 and s.d.
simmr_2_GpIaDg = simmr_load(mixtures=PmSerialSimmData,
                            source_names=as.character(SourceGpIaDg[,1]),
                            source_means=SourceGpIaDg[,c(2,4)],
                            source_sds=SourceGpIaDg[,c(3,5)],
                            correction_means=SourceTEFsd.GpIaDg[,c(2,4)],
                            correction_sds=SourceTEFsd.GpIaDg[,c(3,5)])
plot(simmr_2_GpIaDg, title='Isospace of serial depredators w/ all 9 prey + D.gigas')+coord_flip()
simmr_2_GpIaDg_out=simmr_mcmc(simmr_2_GpIaDg)
summary(simmr_2_GpIaDg_out, type='diagnostics')
simmr_combine2GpIaDg_out<-combine_sources(simmr_2_GpIaDg_out, 
                                          to_combine=c('Sablefish', 'Spiny Dogfish'),
                                          new_source_name='Sablefish+SpinyDogfish')
plot(simmr_2_GpIaDg_out, type='boxplot')
plot(simmr_combine2GpIaDg_out,type='boxplot',title='Serial Depr, combined sablefish & spiny dogfish')

### SIMMR 3 - Non-Serial Depredators ###
#Non-Serial Depredators, Main 7 Prey, TEF of 2.12 and s.d.
simmr_3 = simmr_load(mixtures=PmNonSerialSimmData,
                     source_names=as.character(Source[,1]),
                     source_means=Source[,c(2,4)],
                     source_sds=Source[,c(3,5)],
                     correction_means=SourceTEFsd[,c(2,4)],
                     correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_3, title='Isospace of non-serial depredators - sd')+coord_flip()
simmr_3_out = simmr_mcmc(simmr_3) #Non-Serial Depredators
summary(simmr_3_out, type='diagnostics')
simmr_combine3_out<-combine_sources(simmr_3_out, 
                                    to_combine=c('Sablefish', 'Spiny Dogfish'),
                                    new_source_name='Sablefish+SpinyDogfish')
plot(simmr_combine3_out$input, title="Isospace of Non-Serial Depredators, combine Sblfsh & Dgfsh") + coord_flip()
plot(simmr_3_out, type = 'boxplot')  #non-serial depredators with sd
plot(simmr_3_out, type = 'matrix')  # or type = 'histogram' / 'density' / or 'convergence'
plot(simmr_combine3_out,type='boxplot',title='Non-Serial Depr, combined sablefish & spiny dogfish')

#Non-Serial Depredators, Main 7 Prey + Glass Squid, TEF of 2.12 and s.d.
simmr_3_Gp = simmr_load(mixtures=PmNonSerialSimmData,
                        source_names=as.character(SourceGp[,1]),
                        source_means=SourceGp[,c(2,4)],
                        source_sds=SourceGp[,c(3,5)],
                        correction_means=SourceTEFsd.Gp[,c(2,4)],
                        correction_sds=SourceTEFsd.Gp[,c(3,5)])
simmr_3_Gp_out=simmr_mcmc(simmr_3_Gp)
summary(simmr_3_Gp_out, type='diagnostics')
simmr_combine3_Gp_out<-combine_sources(simmr_3_Gp_out, 
                                       to_combine=c('Sablefish', 'Spiny Dogfish'),
                                       new_source_name='Sablefish+SpinyDogfish')

simmr_3_GpIa = simmr_load(mixtures=PmNonSerialSimmData,
                          source_names=as.character(SourceGpIa[,1]),
                          source_means=SourceGpIa[,c(2,4)],
                          source_sds=SourceGpIa[,c(3,5)],
                          correction_means=SourceTEFsd.GpIa[,c(2,4)],
                          correction_sds=SourceTEFsd.GpIa[,c(3,5)])
simmr_3_GpIa_out=simmr_mcmc(simmr_3_GpIa)
summary(simmr_3_GpIa_out, type='diagnostics')
plot(simmr_3_GpIa_out, type = 'boxplot', title=NULL)

#Non-Serial Depredators, All 9 prey + Humboldt Squid, TEF of 2.12 and s.d.:
simmr_3_GpIaDg = simmr_load(mixtures=PmNonSerialSimmData,
                            source_names=as.character(SourceGpIaDg[,1]),
                            source_means=SourceGpIaDg[,c(2,4)],
                            source_sds=SourceGpIaDg[,c(3,5)],
                            correction_means=SourceTEFsd.GpIaDg[,c(2,4)],
                            correction_sds=SourceTEFsd.GpIaDg[,c(3,5)])
simmr_3_GpIaDg_out=simmr_mcmc(simmr_3_GpIaDg)
summary(simmr_3_GpIaDg_out, type='diagnostics')
plot(simmr_3_GpIaDg_out,type='boxplot',title='Non-Serial Depredators, All Prey')

### SIMMR 4 - Recent Samples ###
# Recent Samples, Main 7 Prey, TEF of 2.12 and s.d.
simmr_4 = simmr_load(mixtures=PmRecentSimmData,
                     source_names=as.character(Source[,1]),
                     source_means=Source[,c(2,4)],
                     source_sds=Source[,c(3,5)],
                     correction_means=SourceTEFsd[,c(2,4)],
                     correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_4, title='Isospace of recent samples - sd')+coord_flip()
simmr_4_out = simmr_mcmc(simmr_4) #Recent (2010-2018) samples
summary(simmr_4_out, type='diagnostics')
simmr_combine4_out<-combine_sources(simmr_4_out, 
                                    to_combine=c('Sablefish', 'Spiny Dogfish'),
                                    new_source_name='Sablefish+SpinyDogfish')
plot(simmr_4_out, type = 'boxplot') 
plot(simmr_4_out, type = 'matrix')
plot(simmr_4_out, type = 'histogram')
plot(simmr_4_out, type="density")
plot(simmr_4_out, type="convergence")
plot(simmr_combine4_out,type='boxplot',title='Recent Whales (2010-2018), combined sablefish & spiny dogfish')

#Recent samples, Main 7 prey + Glass Squid, TEF of 2.12 and s.d.
simmr_4.Gp = simmr_load(mixtures=PmRecentSimmData,
                        source_names=as.character(SourceGp[,1]),
                        source_means=SourceGp[,c(2,4)],
                        source_sds=SourceGp[,c(3,5)],
                        correction_means=SourceTEFsd.Gp[,c(2,4)],
                        correction_sds=SourceTEFsd.Gp[,c(3,5)])
simmr_4_Gp_out = simmr_mcmc(simmr_4.Gp)
simmr_combine4_Gp_out<-combine_sources(simmr_4_Gp_out, 
                                       to_combine=c('Sablefish', 'Spiny Dogfish'),
                                       new_source_name='Sablefish+SpinyDogfish')

#Recent samples, All 9 prey, TEF of 2.12 and s.d.
simmr_4.GpIa = simmr_load(mixtures=PmRecentSimmData,
                          source_names=as.character(SourceGpIa[,1]),
                          source_means=SourceGpIa[,c(2,4)],
                          source_sds=SourceGpIa[,c(3,5)],
                          correction_means=SourceTEFsd.GpIa[,c(2,4)],
                          correction_sds=SourceTEFsd.GpIa[,c(3,5)])
plot(simmr_4.GpIa, title='Isospace of recent samples - full 9 prey')+coord_flip()
simmr_4_GpIa_out = simmr_mcmc(simmr_4.GpIa)
summary(simmr_4_GpIa_out, type='diagnostics')
plot(simmr_4_GpIa_out, type = 'boxplot', title=NULL)

### SIMMR 5 - Old Samples ###
#Old samples, main 7 prey, TEF of 2.12 and s.d.
simmr_5 = simmr_load(mixtures=PmOldSimmData,
                     source_names=as.character(Source[,1]),
                     source_means=Source[,c(2,4)],
                     source_sds=Source[,c(3,5)],
                     correction_means=SourceTEFsd[,c(2,4)],
                     correction_sds=SourceTEFsd[,c(3,5)])
plot(simmr_5, title='Isospace of old samples')+coord_flip()
simmr_5_out = simmr_mcmc(simmr_5) #Old (2003-2009) samples
simmr_combine5_out<-combine_sources(simmr_5_out, 
                                    to_combine=c('Sablefish', 'Spiny Dogfish'),
                                    new_source_name='Sablefish+SpinyDogfish')
plot(simmr_5_out, type = 'boxplot')  
plot(simmr_5_out, type = 'matrix')
plot(simmr_5_out, type = 'histogram')
plot(simmr_5_out, type="density")
plot(simmr_5_out, type="convergence")
plot(simmr_combine5_out,type='boxplot',title='Older samples (2003-2009), combined sablefish & spiny dogfish')

#Old samples, main 7 prey + Glass Squid, TEF of 2.12 and s.d.
simmr_5.Gp = simmr_load(mixtures=PmOldSimmData,
                        source_names=as.character(SourceGp[,1]),
                        source_means=SourceGp[,c(2,4)],
                        source_sds=SourceGp[,c(3,5)],
                        correction_means=SourceTEFsd.Gp[,c(2,4)],
                        correction_sds=SourceTEFsd.Gp[,c(3,5)])
simmr_5_Gp_out = simmr_mcmc(simmr_5.Gp)
simmr_combine5_Gp_out<-combine_sources(simmr_5_Gp_out, 
                                       to_combine=c('Sablefish', 'Spiny Dogfish'),
                                       new_source_name='Sablefish+SpinyDogfish')

#Old samples, All 9 prey, TEF of 2.12 and s.d. 
simmr_5.GpIa = simmr_load(mixtures=PmOldSimmData,
                          source_names=as.character(SourceGpIa[,1]),
                          source_means=SourceGpIa[,c(2,4)],
                          source_sds=SourceGpIa[,c(3,5)],
                          correction_means=SourceTEFsd.GpIa[,c(2,4)],
                          correction_sds=SourceTEFsd.GpIa[,c(3,5)])
plot(simmr_5.GpIa, title='Isospace of old samples - full 9 prey')+coord_flip()
simmr_5_GpIa_out = simmr_mcmc(simmr_5.GpIa) #Old (2003-2009) samples
plot(simmr_5_GpIa_out, type='boxplot', title=NULL)

### SIMMR 6 - Seasonal Whales ### 
#Seasonal Whales, Main 7 Prey, TEF of 2.12 and s.d.
PmSeasonSimmData2 <- as.matrix(PmSeasonSimmData)
simmr_6 = simmr_load(mixtures=PmSeasonSimmData2,
                     source_names=as.character(Source[,1]),
                     source_means=Source[,c(2,4)],
                     source_sds=Source[,c(3,5)],
                     correction_means=SourceTEFsd[,c(2,4)],
                     correction_sds=SourceTEFsd[,c(3,5)],
                     #concentration_means=SourceConcDepen[,c(3,5)],
                     group=as.integer(PmSeasonSimmGrp))
plot(simmr_6,group=1:3,ylab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     xlab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of Sperm Whale seasonal group data', mix_name='Pm')+coord_flip()
simmr_6_out = simmr_mcmc(simmr_6)
summary(simmr_6_out, type='diagnostics')
summary(simmr_6_out, type='quantiles',group=1:3)
simmr_6_groups_out_combineAfSa = combine_sources(simmr_6_out, to_combine=c('Sablefish','Spiny Dogfish'),new_source_name='Sablefish+SpinyDogfish')
plot(simmr_6_groups_out_combineAfSa$input,group=1:3, title="Isospace of all whales grouped + combined sblfsh & dgfsh") + coord_flip()
summary(simmr_6_groups_out_combineAfSa, type='quantiles', group=1:3)
compare_groups(simmr_6_out, source='Sablefish', groups=1:3) #builds boxplots for each group
compare_groups(simmr_6_out, source='Shortraker Rockfish', groups=1:3)
compare_groups(simmr_6_out, source='Clubhook Squid', groups=1:3)
compare_groups(simmr_6_out, source='Spiny Dogfish', groups=1:3)
compare_groups(simmr_6_out, source='Skate', groups=1:3)
compare_groups(simmr_6_out, source='Grenadier', groups=1:3)
compare_groups(simmr_6_out, source='Magister Squid', groups=1:3)

plot(simmr_6_groups_out_combineAfSa,type='boxplot',title='Late season whales (combined Sablefish+SpinyDogfish)',group=3)
plot(simmr_6_groups_out_combineAfSa,type='boxplot',title='Mid season whales (combined Sablefish+SpinyDogfish)',group=2)
plot(simmr_6_groups_out_combineAfSa,type='boxplot',title='Early season whales (combined Sablefish+SpinyDogfish)',group=1)
compare_groups(simmr_6_groups_out_combineAfSa, source="Sablefish+SpinyDogfish", groups=1:3)

#Seasonal Whales, Main 7 Prey + Glass Squid, TEF of 2.12 and s.d.
simmr_6.Gp = simmr_load(mixtures=PmSeasonSimmData2,
                        source_names=as.character(SourceGp[,1]),
                        source_means=SourceGp[,c(2,4)],
                        source_sds=SourceGp[,c(3,5)],
                        correction_means=SourceTEFsd.Gp[,c(2,4)],
                        correction_sds=SourceTEFsd.Gp[,c(3,5)],
                        #concentration_means=SourceConcDepen[,c(3,5)],
                        group=as.integer(PmSeasonSimmGrp))
simmr_6_Gp_out = simmr_mcmc(simmr_6.Gp)
summary(simmr_6_Gp_out, type='diagnostics')
simmr_6_Gp_out_combineAfSa = combine_sources(simmr_6_Gp_out, to_combine=c('Sablefish','Spiny Dogfish'),new_source_name='Sablefish+SpinyDogfish')
summary(simmr_6_Gp_out, group=1, type='statistics')
summary(simmr_6_Gp_out, group=2, type='statistics')
summary(simmr_6_Gp_out, group=3, type='statistics')
summary(simmr_6_Gp_out_combineAfSa, group=1, type='statistics')
summary(simmr_6_Gp_out_combineAfSa, group=2, type='statistics')
summary(simmr_6_Gp_out_combineAfSa, group=3, type='statistics')

#Seasonal whales, All 9 Prey, TEF of 2.12 and s.d.
simmr_6.GpIa = simmr_load(mixtures=PmSeasonSimmData2,
                          source_names=as.character(SourceGpIa[,1]),
                          source_means=SourceGpIa[,c(2,4)],
                          source_sds=SourceGpIa[,c(3,5)],
                          correction_means=SourceTEFsd.GpIa[,c(2,4)],
                          correction_sds=SourceTEFsd.GpIa[,c(3,5)],
                          #concentration_means=SourceConcDepen[,c(3,5)],
                          group=as.integer(PmSeasonSimmGrp))
simmr_6_GpIa_out = simmr_mcmc(simmr_6.GpIa)
summary(simmr_6_GpIa_out, type='diagnostics')
plot(simmr_6_GpIa_out,type='boxplot',title='Late season whales',group=3)
plot(simmr_6_GpIa_out,type='boxplot',title='Mid season whales',group=2)
plot(simmr_6_GpIa_out,type='boxplot',title='Early season whales',group=1)
.Early_out<-as.data.frame(as.matrix(simmr_6_GpIa_out$output[[1]]))
.Mid_out<-as.data.frame(as.matrix(simmr_6_GpIa_out$output[[1]]))
.Late_out<-as.data.frame(as.matrix(simmr_6_GpIa_out$output[[1]]))
summary(simmr_6_GpIa_out, group=1, type='statistics')
summary(simmr_6_GpIa_out, group=2, type='statistics')
summary(simmr_6_GpIa_out, group=3, type='statistics')

simmr_6_GpIaDg = simmr_load(mixtures=PmSeasonSimmData2,
                            source_names=as.character(SourceGpIaDg[,1]),
                            source_means=SourceGpIaDg[,c(2,4)],
                            source_sds=SourceGpIaDg[,c(3,5)],
                            correction_means=SourceTEFsd.GpIaDg[,c(2,4)],
                            correction_sds=SourceTEFsd.GpIaDg[,c(3,5)],
                            group=as.integer(PmSeasonSimmGrp))
plot(simmr_6_GpIaDg,group=1:3,ylab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     xlab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of Sperm Whale seasonal group data with All Potential Prey',mix_name='Pm') +coord_flip()
simmr_6_GpIaDg_out = simmr_mcmc(simmr_6_GpIaDg)
simmr_6_GpIaDg_groups_out_combineSaIa = combine_sources(simmr_6_GpIaDg_out, to_combine=c('Spiny Dogfish','Ragfish'),new_source_name='SpinyDogfish+Ragfish')
plot(simmr_6_GpIaDg_groups_out_combineSaIa$input,group=1:3, title="Isospace of all whales grouped + combined dgfsh &rgfsh") + coord_flip()
plot(simmr_1_GpIaDg_out, type='boxplot')
plot(simmr_combine1GpIaDg_out,type='boxplot',title='All Whales, combined sablefish & spiny dogfish')
plot(simmr_1_GpIaDg_out, type='matrix')

#Seasonal Whales, All Prey including Humboldt Squid with Ragfish & Dogfish Combined, TEF of 2.12 and s.d.
simmr_6_GpDg_IaSa = simmr_load(mixtures=PmSeasonSimmData2,
                               source_names=as.character(SourceGpDg_IaSaCombined[,1]),
                               source_means=SourceGpDg_IaSaCombined[,c(2,4)],
                               source_sds=SourceGpDg_IaSaCombined[,c(3,5)],
                               correction_means=SourceTEFsd.GpDg.IaSaCombined[,c(2,4)],
                               correction_sds=SourceTEFsd.GpDg.IaSaCombined[,c(3,5)],
                               group=as.integer(PmSeasonSimmGrp))
plot(simmr_6_GpDg_IaSa,group=1:3,ylab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     xlab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of Sperm Whale seasonal group data with All Potential Prey, Sa&Ia Combined',mix_name='Pm')+coord_flip()
simmr_6_GpDg_IaSa_out = simmr_mcmc(simmr_6_GpDg_IaSa)
simmr_6_GpDg_IaSa_groups_out_combineAfSaIa = combine_sources(simmr_6_GpDg_IaSa_out,to_combine=c('Sablefish','Dogfish & Ragfish'),new_source_name='Sablefish & Dogfish & Ragfish')
plot(simmr_6_GpDg_IaSa_groups_out_combineAfSaIa$input,group=1:3, title="Isospace of all whales grouped + combined dgfsh &rgfsh") + coord_flip()

###############################################
# Write tables with output of proportion prey #
### --------------------------------------- ###

### SIMMR 1 Data - ALL WHALES ####
.AllWhales_out<-as.data.frame(as.matrix(simmr_1_GpIa_out$output[[1]]))
AllWhalesProp.GpIa<- as.data.frame(as.matrix(colMeans(.AllWhales_out))) #should be same values as if you did summary(simmr_x_out, type='statistics)
write.table(AllWhalesProp.GpIa, file="AllWhalesProp.csv", sep=",")
AllWhalesProp.GpIa2 <- subset(AllWhalesProp.GpIa, select=c('Clubhook Squid', "Grenadier", "Magister Squid", "Sablefish", "Shortraker","Skate", "Spiny Dogfish", "Ragfish", "Glass Squid"))
View(AllWhalesProp.GpIa)
.AllWhales2_out<-as.data.frame(as.matrix(simmr_1_Gp_out$output[[1]]))
AllWhalesProp.Gp<- as.data.frame(as.matrix(colMeans(.AllWhales2_out)))
write.table(AllWhalesProp.Gp, file="AllWhalesProp2.csv", sep=",")
View(AllWhalesProp.Gp)
.AllWhalesComb_out<-as.data.frame(as.matrix(simmr_combine1_out$output[[1]]))
AllWhalesProp.Comb<- as.data.frame(as.matrix(colMeans(.AllWhalesComb_out)))
write.table(AllWhalesProp.Comb, file="AllWhalesPropComb.csv", sep=",")
View(AllWhalesProp.Comb)

### SIMMR 2 - Data - FREQUENT DEPREDATORS ###
.FreqDep_out<-as.data.frame(as.matrix(simmr_2_GpIa_out$output[[1]]))
.NonFreqDep_out<-as.data.frame(as.matrix(simmr_3_GpIa_out$output[[1]]))
FreqDepProp.GpIa<- as.data.frame(as.matrix(colMeans(.FreqDep_out)))
NonFreqDepProp.GpIa<- as.data.frame(as.matrix(colMeans(.NonFreqDep_out)))
write.table(FreqDepProp.GpIa, file="FreqDepProp.csv", sep=",")
write.table(NonFreqDepProp.GpIa, file="NonFreqDepProp.csv", sep=",")
.FreqDep2_out<-as.data.frame(as.matrix(simmr_2_Gp_out$output[[1]]))
.NonFreqDep2_out<-as.data.frame(as.matrix(simmr_3_Gp_out$output[[1]]))
FreqDepProp.Gp<- as.data.frame(as.matrix(colMeans(.FreqDep2_out)))
NonFreqDepProp.Gp<- as.data.frame(as.matrix(colMeans(.NonFreqDep2_out)))
write.table(FreqDepProp.Gp, file="FreqDepProp2.csv", sep=",")
write.table(NonFreqDepProp.Gp, file="NonFreqDepProp2.csv", sep=",")

.FreqDepComb_out<-as.data.frame(as.matrix(simmr_combine2_Gp_out$output[[1]]))
.NonFreqDepComb_out<-as.data.frame(as.matrix(simmr_combine3_Gp_out$output[[1]]))
FreqDepProp.Comb<- as.data.frame(as.matrix(colMeans(.FreqDepComb_out)))
NonFreqDepProp.Comb<- as.data.frame(as.matrix(colMeans(.NonFreqDepComb_out)))
write.table(FreqDepProp.Comb, file="FreqDepPropComb.csv", sep=",")
write.table(NonFreqDepProp.Comb, file="NonFreqDepPropComb.csv", sep=",")

### SIMRR 5 Data -  OLD SAMPLES ### 
.Old_out<-as.data.frame(as.matrix(simmr_5_GpIa_out$output[[1]]))
.Recent_out<-as.data.frame(as.matrix(simmr_4_GpIa_out$output[[1]]))
OldProp.GpIa<- as.data.frame(as.matrix(colMeans(.Old_out)))
RecentProp.GpIa<- as.data.frame(as.matrix(colMeans(.Recent_out)))
write.table(OldProp.GpIa, file="OldProp.csv", sep=",")
write.table(RecentProp.GpIa, file="RecentProp.csv", sep=",")
.Old2_out<-as.data.frame(as.matrix(simmr_5_Gp_out$output[[1]]))
.Recent2_out<-as.data.frame(as.matrix(simmr_4_Gp_out$output[[1]]))
OldProp.Gp<- as.data.frame(as.matrix(colMeans(.Old2_out)))
RecentProp.Gp<- as.data.frame(as.matrix(colMeans(.Recent2_out)))
write.table(OldProp.Gp, file="OldProp2.csv", sep=",")
write.table(RecentProp.Gp, file="RecentProp2.csv", sep=",")
.OldComb_out<-as.data.frame(as.matrix(simmr_combine5_Gp_out$output[[1]]))
.RecentComb_out<-as.data.frame(as.matrix(simmr_combine4_Gp_out$output[[1]]))
OldProp.Comb<- as.data.frame(as.matrix(colMeans(.OldComb_out)))
RecentProp.Comb<- as.data.frame(as.matrix(colMeans(.RecentComb_out)))
write.table(OldProp.Comb, file="OldPropComb.csv", sep=",")
write.table(RecentProp.Comb, file="RecentPropComb.csv", sep=",")


######################################
### PLOT POINT ESTIMATES OVER TIME ### From Mike Sigler instructions for a plot... 
### ############################## ###

allwhale.prop.output <- data.frame(matrix(ncol=6, nrow=7))
.x<- c("Species", "AllWhales", "Serial", "Non-Serial", "Old", "Recent")
colnames(allwhale.prop.output)<-.x
allwhale.prop.output$Species <- c("Skate", "Sablefish", "Dogfish", "Grenadier", "Clubhook", "Shortraker", "Magister")
.skate_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$Skate
.sablefish_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$Sablefish
.dogfish_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$"Spiny Dogfish" #because there is a space in the species, need to put it in quotations
.grenadier_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$Grenadier
.clubhook_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$"Clubhook Squid"
.shortraker_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$"Shortraker Rockfish"
.magister_out<-as.data.frame(as.matrix(simmr_1_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$AllWhales<-c(mean(.skate_out),mean(.sablefish_out),mean(.dogfish_out),mean(.grenadier_out),mean(.clubhook_out),mean(.shortraker_out),mean(.magister_out))
.skate_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$Skate
.sablefish_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$Sablefish
.dogfish_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$"Spiny Dogfish"
.grenadier_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$Grenadier
.clubhook_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$"Clubhook Squid"
.shortraker_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$"Shortraker Rockfish"
.magister_2out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$Serial<-c(mean(.skate_2out),mean(.sablefish_2out),mean(.dogfish_2out),mean(.grenadier_2out),mean(.clubhook_2out),mean(.shortraker_2out),mean(.magister_2out))
.skate_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$Skate
.sablefish_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$Sablefish
.dogfish_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$"Spiny Dogfish"
.grenadier_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$Grenadier
.clubhook_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$"Clubhook Squid"
.shortraker_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$"Shortraker Rockfish"
.magister_3out<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$`Non-Serial`<-c(mean(.skate_3out),mean(.sablefish_3out),mean(.dogfish_3out),mean(.grenadier_3out),mean(.clubhook_3out),mean(.shortraker_3out),mean(.magister_3out))
.skate_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$Skate
.sablefish_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$Sablefish
.dogfish_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$"Spiny Dogfish"
.grenadier_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$Grenadier
.clubhook_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$"Clubhook Squid"
.shortraker_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$"Shortraker Rockfish"
.magister_4out<-as.data.frame(as.matrix(simmr_4_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$Recent<-c(mean(.skate_4out),mean(.sablefish_4out),mean(.dogfish_4out),mean(.grenadier_4out),mean(.clubhook_4out),mean(.shortraker_4out),mean(.magister_4out))
.skate_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Skate
.sablefish_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Sablefish
.dogfish_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Spiny Dogfish"
.grenadier_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Grenadier
.clubhook_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Clubhook Squid"
.shortraker_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Shortraker Rockfish"
.magister_5out<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$Old<-c(mean(.skate_5out),mean(.sablefish_5out),mean(.dogfish_5out),mean(.grenadier_5out),mean(.clubhook_5out),mean(.shortraker_5out),mean(.magister_5out))
write.table(allwhale.prop.output, file="allwhale.prop.output.csv", sep=",") #export it and put it in the right format; 
allwhale.prop.output2<- read.csv('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/allwhale.prop.output2.csv',sep=",",header=TRUE)

shape1<-c(16,16,24,23,16)
color1<-c("purple", "orange", "black", "black","blue")
ggplot(allwhale.prop.output2, aes(Species, Proportion, shape=Group, color=Group)) + geom_point(size=4) + xlab("Species") + ylab("Proportion of Diet") +
  scale_shape_manual(values=shape1)+
  scale_color_manual(values=color1)+
  theme_bw(base_size=24, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

.skate_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Skate
.sablefish_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Sablefish
.dogfish_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Spiny Dogfish"
.grenadier_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Grenadier
.clubhook_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Clubhook Squid"
.shortraker_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Shortraker Rockfish"
.magister_6out_early<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$EarlySumm<-c(mean(.skate_6out_early),mean(.sablefish_6out_early),mean(.dogfish_6out_early),mean(.grenadier_6out_early),mean(.clubhook_6out_early),mean(.shortraker_6out_early),mean(.magister_6out_early))
.skate_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Skate
.sablefish_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Sablefish
.dogfish_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Spiny Dogfish"
.grenadier_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Grenadier
.clubhook_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Clubhook Squid"
.shortraker_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Shortraker Rockfish"
.magister_6out_mid<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$MidSumm<-c(mean(.skate_6out_mid),mean(.sablefish_6out_mid),mean(.dogfish_6out_mid),mean(.grenadier_6out_mid),mean(.clubhook_6out_mid),mean(.shortraker_6out_mid),mean(.magister_6out_mid))
.skate_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Skate
.sablefish_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Sablefish
.dogfish_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Spiny Dogfish"
.grenadier_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$Grenadier
.clubhook_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Clubhook Squid"
.shortraker_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Shortraker Rockfish"
.magister_6out_late<-as.data.frame(as.matrix(simmr_5_out$output[[1]]))$"Magister Squid"
allwhale.prop.output$LateSumm<-c(mean(.skate_6out_late),mean(.sablefish_6out_late),mean(.dogfish_6out_late),mean(.grenadier_6out_late),mean(.clubhook_6out_late),mean(.shortraker_6out_late),mean(.magister_6out_late))

as.data.frame(as.matrix(simmr_6_out$output[[1]]))$Grenadier
grpwhale.prop.output<-data.frame(Species=c("Clubhook", "Grenadier", "Magister", "Sablefish", "Shortraker", "Skate", "Dogfish","Clubhook", "Grenadier", "Magister", "Sablefish", "Shortraker", "Skate", "Dogfish","Clubhook", "Grenadier", "Magister", "Sablefish", "Shortraker", "Skate", "Dogfish"),Proportion = c(0.135, 0.073, 0.071, 0.172, 0.139, 0.259, 0.151, 0.111,0.129,0.138,0.158,0.174,0.125,0.163,0.073,0.052,0.067,0.267,0.095,0.155,0.290), Group=c("Early","Early","Early","Early","Early","Early","Early","Mid","Mid","Mid","Mid","Mid","Mid","Mid","Late","Late","Late","Late","Late","Late","Late")) 
color2<-c("purple", "orange", "blue")
ggplot(grpwhale.prop.output, aes(Species, Proportion,color=Group)) + geom_point(size=4) + 
  xlab("Species") + ylab("Proportion of Diet") +
  #scale_color_manual(values=color2)+
  scale_color_viridis_d()+
  theme_bw(base_size=24, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

################################# 
library(boot)
#### This code used to build loop that selects 10 rand samples from 33 for Freq/Non-Freq comparisons ####
### Authors: Justin Priest & Lauren Wild
### Date: January 2019
#change mcmc.control to fewer iterations (maybe 1000), keep burn in 1000, thin to 100 or 1000
#then put simmr load and simmr out into for loop, 
#save the model output from simmr_out to an empty data frame, or just
#grab something from the output; probably mean for each species separately. 

## First here is Justin Priest's code to set up an example data set and start the loop work:
n15_nodep <- rnorm(23, 17, 2) #create nitrogen isotope data set with 23 values (non-serial depr)
n15_dep <- rnorm(10, 18, 2) #create nitrogen isotope data set with 10 values (serial depr)
c13_nodep <- rnorm(23, -16, 2) #create carbon isotope data set for non-serial depr 
c13_dep <- rnorm(10, -18, 2)  #create carbon isotope data set for serial depr.
wholedata <- data.frame(n15 = c(n15_dep, n15_nodep), c13= c(c13_dep, c13_nodep)) #merge the data frames
depdata <- data.frame(n15 = n15_dep, c13= c13_dep) #isolate just the "serial depredators" data fram
loops <- 100 #can adjust this and will run in loop so you don't have to adjust multiple times within the loop
newoutput <- data.frame(pvaln = rep(NA, loops), pvalc = rep(NA, loops)) #build an output data frame with the same number of rows as there are "loops"
for(i in 1:loops){
  .rowselect <- round(runif(10, 1, 33)) #randomly select 10 rows to use (number from 1 to 33)
  .subset10 <- wholedata[c(.rowselect),] # this selects our 10 rows (note, the '.' at the beginning just means it isn't saved into your workspace, clutterig it up)
  .subset23 <- wholedata[-c(.rowselect),] # this gives the remaining 23 rows to use
  newoutput[i,1] <- t.test(.subset10$n15, .subset23$n15)$p.value
  newoutput[i,2] <- t.test(.subset10$c13, .subset23$c13)$p.value
}
sum(newoutput$pvaln < 0.05) #number of times in bootsrapt the p-value was significant
sum(newoutput$pvaln<0.05)/ loops #no. of times p-value was sig / no. of times it was bootstrapped

# Now try to do the loop with my data
loops <- 1000
newoutput <- data.frame(pval = rep(NA, loops))
newoutput2 <- data.frame(MedianDiffSkate = rep(NA, loops), MedianDiffSablefish = rep(NA, loops), MedianDiffDogfish = rep(NA, loops), MedianDiffGrenadier = rep(NA, loops), MedianDiffClubhook = rep(NA, loops), MedianDiffShortraker = rep(NA, loops), MedianDiffMagister = rep(NA, loops))
for(i in 1:loops){
  .rowselect <- round(runif(10, 1, 33))
  .subset10 <- as.matrix(PmSimmData2[c(.rowselect),]) # this gives us 10
  .subset23 <- as.matrix(PmSimmData2[-c(.rowselect),]) # this gives the remaining 23 - somehow giving 25?!?
  simmr_10loop <- simmr_load(mixtures=.subset10,  
                             source_names=as.character(Source[,1]),
                             source_means=Source[,c(2,4)],
                             source_sds=Source[,c(3,5)],
                             correction_means=SourceTEFsd[,c(2,4)],
                             correction_sds=SourceTEFsd[,c(3,5)])
  simmr_23loop <- simmr_load(mixtures=.subset23,
                             source_names=as.character(Source[,1]),
                             source_means=Source[,c(2,4)],
                             source_sds=Source[,c(3,5)],
                             correction_means=SourceTEFsd[,c(2,4)],
                             correction_sds=SourceTEFsd[,c(3,5)])
  simmr_10loop_out = simmr_mcmc(simmr_10loop, mcmc.control=list(iter=5000,burn=1000,thin=100, n.chain=4))
  simmr_23loop_out = simmr_mcmc(simmr_23loop, mcmc.control=list(iter=5000,burn=1000,thin=100, n.chain=4))
  
  .skate10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$Skate
  .skate23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$Skate
  newoutput2[i,1]<- median(.skate10_out) - median(.skate23_out) #difference between median of random 10 and leftover 23
  .sablefish10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$Sablefish
  .sablefish23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$Sablefish
  newoutput2[i,2]<- median(.sablefish10_out) - median(.sablefish23_out)
  .dogfish10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$"Spiny Dogfish"
  .dogfish23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$"Spiny Dogfish"
  newoutput2[i,3]<- median(.dogfish10_out) - median(.dogfish23_out)
  .grenadier10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$Grenadier
  .grenadier23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$Grenadier
  newoutput2[i,4]<- median(.grenadier10_out) - median(.grenadier23_out)
  .clubhook10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$"Clubhook Squid"
  .clubhook23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$"Clubhook Squid"
  newoutput2[i,5]<- median(.clubhook10_out) - median(.clubhook23_out)
  .shortraker10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$"Shortraker Rockfish"
  .shortraker23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$"Shortraker Rockfish"
  newoutput2[i,6]<- median(.shortraker10_out) - median(.shortraker23_out)
  .magister10_out<-as.data.frame(as.matrix(simmr_10loop_out$output[[1]]))$"Magister Squid"
  .magister23_out<-as.data.frame(as.matrix(simmr_23loop_out$output[[1]]))$"Magister Squid"
  newoutput2[i,7]<- median(.magister10_out) - median(.magister23_out)
  #newoutput[i,1] <- t.test(skate10_out, skateserial_out)$p.value                        
}
#Took 15min to run 200 loops
#To run 1000 loops, took 1hr 20min

newoutput2$Max <- apply(abs(newoutput2), 1, max)
hist(newoutput2$Max) #where does my serial-nonserial difference in skate proportion () fall? 
.skate10_serial<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$Skate
.skate23_nonserial<-as.data.frame(as.matrix(simmr_3_out$output[[1]]))$Skate
median(.skate10_serial) - median(.skate23_nonserial) #0.2611831 is the difference between my 10 serial and 23 non-serial for skate consumption
0.05*loops #tells the nth highest number to compare my serial vs. non-serial with
library(Rfast)
Rfast::nth(newoutput2$Max, 50, descending=T) #should give the 50th highest max difference
#Want to see if my difference in serial skate proportion-non-serial skate prop is higher or lower than the top 5th percent of random loops
#Get 0.23985, and my diff is 0.26118 so my skate diff is significantly to if we randomly select 10 indiv. & calculate diff in proportions from max difference prey. 
Rfast::nth(newoutput2$MedianDiffSkate, 50, descending=T) #0.2398497
#If I just isolate skate max difference for 10 rand samp vs. remaining 23, the 5% cutoff is 0.239897

#summary(simmr_10loop_out, type='statistics') #gives you the means and sd of each proportion for each species from that mcmc run
#skateserial_out<-as.data.frame(as.matrix(simmr_2_out$output[[1]]))$Skate #this gives you a dataframe with just the raw skate proportions from the mcmc run in simmr_mcmc output

str(simmr_10loop_out$output) #give the structure of the r-output so you can find how to reference a p-value, or a species or something
sum(newoutput$pvaln < 0.05) #number of times in bootsrapt the p-value was significant
sum(newoutput$pvaln<0.05)/ loops #no. of times p-value was sig / no. of times it was bootstrapped


#####################################################
########### COMPARE MODELS WITH MIX-SIAR ############
#####################################################
# Authors: Brian Stock [cre, aut], Brice Semmens [aut], Eric Ward [ctb], Andrew Parnell [ctb], Andrew Jackson [ctb], Donald Phillips [ctb]
# Maintained by Brian Stock, b1stock@ucsd.edu
# Uses rjags and JAGS
# https://github.com/brianstock/MixSIAR
# cite: Stock, B.C. and B.X. Semmens (2016). MixSIAR GUI User Manual. Version 3.1. https://github.com/brianstock/MixSIAR/. doi: 10.5281/zenodo.47719
# https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf
browseVignettes("MixSIAR") #shows examples with scripts

library(MixSIAR)
mixsiar.dir <- find.package("MixSIAR")
paste0(mixsiar.dir,"/example_scripts")
#source(paste0(mixsiar.dir,"/example_scripts/mixsiar_script_killerwhale.R")) #Needs 'splancs' package to run calc_area #run sexample scripts for killer whale example! 
setwd('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/')

#Load in data for MixSIAR:
#Full data set, no factors. 
mixPm <- load_mix_data(filename="PmSimmData2.csv",
                       iso_names=c("d13C","d15N"),
                       factors=NULL,
                       fac_random=NULL,
                       fac_nested=NULL,
                       cont_effects=NULL) #Needs to be in format of just d13C and d15N columns, and "factor" column if there is one. 
#Temporal data with column of "old" or "recent" sample
mixPmTemporal <- load_mix_data(filename="PmMixSIARDataTemporal.csv",
                               iso_names=c("d13C","d15N"),
                               factors='Temporal',
                               fac_random=TRUE,
                               fac_nested=FALSE,
                               cont_effects=NULL)
# "Seasonal" with Early, Mid, and Late season factor variable
mixPmSeasonal <- load_mix_data(filename="PmMixSIARDataSeasonal.csv",
                               iso_names=c("d13C","d15N"),
                               factors='Season',
                               fac_random=TRUE,
                               fac_nested=FALSE,
                               cont_effects=NULL)
mixPmFreqNFreq <- load_mix_data(filename="PmMixSIARDataSerialNSerial.csv",
                                iso_names=c("d13C","d15N"),
                                factors='Frequent',
                                fac_random=TRUE,
                                fac_nested=FALSE,
                                cont_effects=NULL)

sourcePm <- load_source_data(filename="PmSources-GpIa2.csv",
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="means",
                             mixPm)   #Source data frame needs to be in format: "Species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", "n" (sample size each source)
sourcePmTemporal <- load_source_data(filename="PmSources-GpIa2.csv",
                                     source_factors=NULL,
                                     conc_dep=FALSE,
                                     data_type="means",
                                     mixPmTemporal)
sourcePmSeasonal <- load_source_data(filename="PmSources-GpIa2.csv",
                                     source_factors=NULL,
                                     conc_dep=FALSE,
                                     data_type="means",
                                     mixPmSeasonal)
sourcePmFreqNFreq <- load_source_data(filename="PmSources-GpIa2.csv",
                                      source_factors=NULL,
                                      conc_dep=FALSE,
                                      data_type="means",
                                      mixPmFreqNFreq)
discrPm <- load_discr_data(filename="PmSourceTEFsd-GpIa.csv", mixPm)
discrPmTemporal <- load_discr_data(filename="PmSourceTEFsd-GpIa.csv", mixPmTemporal)
discrPmSeasonal <- load_discr_data(filename="PmSourceTEFsd-GpIa.csv", mixPmSeasonal)
discrPmFreqNFreq <- load_discr_data(filename="PmSourceTEFsd-GpIa.csv", mixPmFreqNFreq)

# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPm,sourcePm,discrPm)
plot_data(filename="isospace_plot_Temporal", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmTemporal,sourcePmTemporal,discrPmTemporal)
plot_data(filename="isospace_plot_Seasonal", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmSeasonal,sourcePmSeasonal,discrPmSeasonal)
plot_data(filename="isospace_plot_FreqNFreq", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmFreqNFreq,sourcePmFreqNFreq,discrPmFreqNFreq)

# Plot uninformative prior
plot_prior(alpha.prior=1, sourcePm, filename = "prior_plot_pm_uninf", plot_save_pdf=FALSE)
plot_prior(alpha.prior=1, sourcePmTemporal, filename = "prior_plot_pmTemp_uninf", plot_save_pdf=FALSE)
plot_prior(alpha.prior=1, sourcePmSeasonal, filename = "prior_plot_pmSeason_uninf", plot_save_pdf=FALSE)
plot_prior(alpha.prior=1, sourcePmFreqNFreq, filename = "prior_plot_pmFnF_uninf", plot_save_pdf=FALSE)

# Define model structure and write JAGS model file
model_Pm <- "MixSIAR_model_Pm_uninf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_Pm, resid_err, process_err, mixPm, sourcePm)
model_Pm_Temporal <- "MixSIAR_model_Pm_Temporal_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_Temporal, resid_err, process_err, mixPmTemporal, sourcePmTemporal)
model_Pm_Seasonal <- "MixSIAR_model_Pm_Seasonal_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_Seasonal, resid_err, process_err, mixPmSeasonal, sourcePm)
model_Pm_FreqNFreq <- "MixSIAR_model_Pm_FreqNFreq_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_FreqNFreq, resid_err, process_err, mixPmFreqNFreq, sourcePm)

# Run the JAGS model ("test" took ~20sec; "short" should be good... about 15min, with chain length 50,000, burn in 25,000, and thin 25 and it seems to converge nicely) 
# Pg 14 of MixSIAR manual lists how many MCMC burn-ins and chains each one uses ("very long" is 1,000,000 chain length, 500,000 burn-in, and thin by 500, should take about 5hr to run)
jags.uninf <- run_model(run="short",mixPm,sourcePm,discrPm,model_Pm,alpha.prior = 1, resid_err, process_err) #took 14min to run
jags.uninf.Temporal <- run_model(run="short",mixPmTemporal,sourcePmTemporal,discrPmTemporal, model_Pm_Temporal,alpha.prior = 1, resid_err, process_err) #short took 24 min to run
jags.uninf.Seasonal <- run_model(run="short",mixPmSeasonal,sourcePmSeasonal,discrPmSeasonal,model_Pm_Seasonal,alpha.prior = 1, resid_err, process_err) #took 24min to run
jags.uninf.FreqNFreq <- run_model(run="short",mixPmFreqNFreq,sourcePmFreqNFreq,discrPmFreqNFreq,model_Pm_FreqNFreq,alpha.prior = 1, resid_err, process_err)

SblGrp = combine_sources(jags.uninf, mixPm, sourcePm, alpha.prior=1, groups=list(Clubhook.Squid="Clubhook Squid", Glass.Squid="Glass Squid", Grenadier="Grenadier", Magister.Squid="Magister Squid", SablefishGroup=c("Ragfish","Sablefish","Spiny Dogfish"), Shortraker.Rockfish="Shortraker Rockfish", Skate="Skate")) #Here my new aggregated priors in the plot produced show that the weighting on grups has changed, so I should group the sources before running the mixing model (a priori).

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf, mixPm, sourcePm)
output_JAGS(jags.uninf.Temporal, mixPmTemporal, sourcePm)
output_JAGS(jags.uninf.Seasonal, mixPmSeasonal, sourcePm)
output_JAGS(jags.uninf.FreqNFreq, mixPmFreqNFreq, sourcePm)

compare_models(x=list(mod.uninf = jags.uninf, mod.uninf.temporal = jags.uninf.Temporal), loo=TRUE)
compare_models(x=list(mod.uninf = jags.uninf, mod.uninf.seasonal = jags.uninf.Seasonal), loo=TRUE)
compare_models(x=list(mod.uninf = jags.uninf, mod.unif.FreqNFreq = jags.uninf.FreqNFreq), loo=TRUE) 


#Load in grouping data for MixSIAR:
#Whale Data Sets are the same: mixPm, mixPmTemporal, mixPmSeasonal, mixPmFreqNFreq 
#Source Data:
sourcePm.Comb <- load_source_data(filename="PmSources-GpIa3.csv",
                                  source_factors=NULL,
                                  conc_dep=FALSE,
                                  data_type="means",
                                  mixPm)   #Source data frame needs to be in format: "Species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", "n" (sample size each source)
sourcePm.TempComb <- load_source_data(filename="PmSources-GpIa3.csv",
                                      source_factors=NULL,
                                      conc_dep=FALSE,
                                      data_type="means",
                                      mixPmTemporal)
sourcePm.SeasonComb <- load_source_data(filename="PmSources-GpIa3.csv",
                                        source_factors=NULL,
                                        conc_dep=FALSE,
                                        data_type="means",
                                        mixPmSeasonal)
sourcePm.FreqNFreqComb <- load_source_data(filename="PmSources-GpIa3.csv",
                                           source_factors=NULL,
                                           conc_dep=FALSE,
                                           data_type="means",
                                           mixPmFreqNFreq)
discrPm.Comb <- load_discr_data(filename="PmSourceTEFsd-GpIa-AfIaSaComb.csv", mixPm)
discrPm.TempComb <- load_discr_data(filename="PmSourceTEFsd-GpIa-AfIaSaComb.csv", mixPmTemporal)
discrPm.SeasonComb <- load_discr_data(filename="PmSourceTEFsd-GpIa-AfIaSaComb.csv", mixPmSeasonal)
discrPm.FreqNFreqComb <- load_discr_data(filename="PmSourceTEFsd-GpIa-AfIaSaComb.csv", mixPmFreqNFreq)

# Make an isospace plot
plot_data(filename="isospace_plot_combined", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPm,sourcePm.Comb,discrPm.Comb)
plot_data(filename="isospace_plot_combined_Temporal", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmTemporal,sourcePm.TempComb,discrPm.TempComb)
plot_data(filename="isospace_plot_combined_Seasonal", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmSeasonal,sourcePm.SeasonComb,discrPm.SeasonComb)
plot_data(filename="isospace_plot_combined_FreqNFreq", plot_save_pdf=FALSE, plot_save_png=TRUE, mixPmFreqNFreq,sourcePm.FreqNFreqComb,discrPm.FreqNFreqComb)

# Plot uninformative prior
plot_prior(alpha.prior=1, sourcePm, filename = "prior_plot_pm_comb_uninf", plot_save_pdf=FALSE)
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
jags.uninf.Pm.Comb <- run_model(run="normal",mixPm,sourcePm.Comb,discrPm.Comb,model_Pm.Comb,alpha.prior = 1, resid_err, process_err) #took 36 min to run.
jags.uninf.Pm.TempComb <- run_model(run="normal",mixPmTemporal,sourcePm.TempComb,discrPm.TempComb, model_Pm_TempComb,alpha.prior = 1, resid_err, process_err) #'normal' took 36min 
jags.uninf.Pm.SeasonComb <- run_model(run="normal",mixPmSeasonal,sourcePm.SeasonComb,discrPm.SeasonComb,model_Pm_SeasonComb,alpha.prior = 1, resid_err, process_err) #Took 35min
jags.uninf.Pm.FreqNFreqComb <- run_model(run="normal",mixPmFreqNFreq,sourcePm.FreqNFreqComb, discrPm.FreqNFreqComb,model_Pm_FreqNFreqComb,alpha.prior = 1, resid_err, process_err)

# get posterior medians for new source groupings
apply(combined$post, 2, median)
summary_stat(jags.uninf.Pm.Comb, meanSD=FALSE, quantiles=c(.025,.5,.975), savetxt=FALSE)

# Process diagnostics, summary stats, and posterior plots
jags_output_fullmodel<- output_JAGS(jags.uninf.Pm.Comb, mixPm, sourcePm.Comb) #DIC=118.3458
jags_output_tempcomb<- output_JAGS(jags.uninf.Pm.TempComb, mixPmTemporal, sourcePm.TempComb) #DIC 117.3729 - model doesn't converge super well, might need to run more iterations! 
jags_output_seasoncomb<- output_JAGS(jags.uninf.Pm.SeasonComb, mixPmSeasonal, sourcePm.SeasonComb) #DIC 120.2877 - model did converge pretty well,
jags_output_fnfcomb<- output_JAGS(jags.uninf.Pm.FreqNFreqComb, mixPmFreqNFreq, sourcePm.FreqNFreqComb) #DIC 114.9355 - model converged pretty well

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

#Prey List: 1=Clubhook squid, 2=Glass Squid, 3=Grenadier, 4=Magister Squid, 5=Sablefish Group, 6=Shortraker rockfish, 7=Skate
###################################################################

# Compare Proportion of "Sablefish.Dogfish.Ragfish" vs "Skate"  
# This computes the probability that the proportion of "Skates" of "frequent deperedator" diets is bigger than the dietary proportion of "Sablefish.group" for "freq depr", given the data
# Find out which group is "Frequent" depredators:
freq_level = which(levels(mixPmFreqNFreq$data$Frequent) == "Frequent") #also first column
# Find out which group is "Non-Frequent" depredators:
nonfreq_level = which(levels(mixPmFreqNFreq$data$Frequent) == "Non-Frequent") #also 2nd column
# Find out which source is "Sablefish.Dogfish.Ragfish"
source_1 = which(sourcePm.FreqNFreqComb$source_names == "Sablefish.Dogfish.Ragfish")
# Find out which source is "Skates"
source_2 = which(sourcePm.FreqNFreqComb$source_names == "Skate")

prop_Sablefish.Dogfish.Ragfish_freq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,freq_level,source_1]
prop_Skate_freq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,freq_level,source_2]
prop_Sablefish.Dogfish.Ragfish_nfreq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,nonfreq_level,source_1]
prop_Skate_nfreq = jags.uninf.Pm.FreqNFreqComb$BUGSoutput$sims.list$p.fac1[,nonfreq_level,source_2]

# Create a plot of the differences
diff = prop_Skate_freq - prop_Sablefish.Dogfish.Ragfish_freq
qplot(diff, geom = 'histogram',
      xlab = 'Posterior difference in dietary proportions between Skate and Sablefish/Dogfish for Frequent depredators')
diff2 = prop_Sablefish.Dogfish.Ragfish_nfreq - prop_Skate_nfreq
qplot(diff2, geom = 'histogram',
      xlab = 'Posterior difference in dietary proportions between Skate and Sablefish/Dogfish for non-frequent depredators')

diff3 = prop_Skate_freq - prop_Skate_nfreq
# Find the probability one is bigger than the other
sum(diff>0)/length(diff) # 0.796667
# 80% chance that, in the frequent depredator group, whales are eating a higher proportion Skates than Sablefish/Dogfish.

sum(diff2>0)/length(diff2) # 0.822
# 82% chance that, in the non-frequent depredator group, whales are eating a higher proportion Sablefish group than Skates.

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


#####################################################
####### ------- ISOTOPIC NICHE SPACE ------- ########
# --------------------------------------------------#
install.packages("SIBER")
library(SIBER)
