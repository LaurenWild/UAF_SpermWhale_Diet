# CHAPTER 2, APPENDIX 1. 
# This code creates plots and data analysis used in Appendix 1 to supplement the stable isotope analysis and mixing models of sperm whale diet in the Gulf of Alaska
# Author: Lauren Wild, lawild@alaska.edu; May 2019

library(ggplot2)
library(psych)
library(devtools)
library(plotrix)  #library for the std.error function; can use mean_se() too, same answer
library(MASS)
library(stats)
library(ggpubr)
library(dplyr)
library(viridis) # color scheme for plots that is easy to read (from simmr)
library(jmv) #Has the mancova function
library(doBy)
library(gtools)
library(here)
library(rjags)
library(siar)
library(simmr)
library(MixSIAR)  ### Code starts line 853
library(patternplot) #to put patterns in the fill of boxplots, etc. 


Pm2<-read.table(here::here("data/PmIsotopes2_forR.csv"),sep=",",header=TRUE)

#Need to make Layer a numeric variable: 
Pm2["Layer2"] <- NA #Creates new column
Pm2$Layer2<- 5-as.numeric(Pm2$Layer) # Creates a new column making layer numeric
levels(Pm2$Layer)
str(Pm2$Layer2)
levels(Pm2$Month)

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
#View(Pm2Inner) #Should have dates in Date2 Column! 
doy <- strftime(Pm2Inner$Date2, format = "%j") # new vector of dates in julian format
doy<-as.numeric(doy) #make julian date numeric
Pm2Inner$doy<-doy #add julian date to data frame
View(Pm2Inner) # Check everything looks good! 

##################################################################################


### APPENDIX 1: Show how ragfish fit into species plots; 
### ----------------------------------------------------
Prey3<-read.table(here::here('data/Prey-outliers-removed-Sep2018.csv'),sep=",",header=TRUE)
View(Prey3)
Prey3$Depth.Strata<-as.factor(Prey3$Depth.Strata) 
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
Cy2<-Cy[!is.na(Cy$d15N.bulk),]
Cy2$d15N<-Cy2$d15N.bulk
Sb2<-Sb[!is.na(Sb$d15N.bulk),]
Sb2$d15N<-Sb2$d15N.bulk
Sa2<-Sa[!is.na(Sa$d15N.bulk),]
Sa2$d15N<-Sa2$d15N.bulk
Bm2<-Bm[!is.na(Bm$d15N.bulk),]
Bm2$d15N<-Bm2$d15N.bulk
Ia$d15N<-Ia$d15N.bulk
Gp$d15N<-Gp$d15N.bulk
Or$d15N<-Or$d15N.bulk
Rb$d15N<-Rb$d15N.LE
mean(Ia$d15N.bulk)
sd(Ia$d15N.bulk)
mean(Ia$d13C.LE)
sd(Ia$d13C.LE)

Prey3.2<-rbind(Or,Cy2,Af2,Sb2,Sa2,Rb,Bm2,Ia) 
Prey3.2$d13C<- Prey3.2$d13C.LE
All.Prey.Sum2<- summaryBy(d15N+d13C~Species, data=Prey3.2, FUN=myfun1) #Add Ia
write.table(All.Prey.Sum2, file="PmSourceWithRagfishs.csv", sep=",")

# Add sperm whales in to the top of it
Pm.Prey <- cbind(Prey3.2)

#use sperm whale inner layer data frame
Pm2Inner$Species<-"Sperm Whale" 
Pm2Inner2<-Pm2Inner[,c(1,3,4)]

Prey3.2<-Prey3.2[,c(1,30,31)]
Pm.Prey2<-rbind(Prey3.2, Pm2Inner2) #Top 7 and Ragfish
str(Pm.Prey2)

#set up data sheet with averages and std.error of each species
All.Sp.Sum2 <-summaryBy(d15N+d13C~Species, data=Pm.Prey2, FUN=myfun1) #All 7  
View(All.Sp.Sum2)

Ap1_AllSpRagfish<-ggplot(All.Sp.Sum2,aes(d13C.m,d15N.m, label=Species)) + 
  geom_point(size=8) +
  geom_errorbarh(aes(xmax=All.Sp.Sum2$d13C.m+All.Sp.Sum2$d13C.sd,xmin=All.Sp.Sum2$d13C.m-All.Sp.Sum2$d13C.sd, height = 0.01)) +
  geom_errorbar(aes(ymax=All.Sp.Sum2$d15N.m+All.Sp.Sum2$d15N.sd,ymin=All.Sp.Sum2$d15N.m-All.Sp.Sum2$d15N.sd, width = 0.01))+
  geom_text(color='black', hjust=-0.03,vjust = -0.7, size = 10)+
  xlab(expression(paste(delta^13, "C (\u2030)",sep="")))+ 
  ylab(expression(paste(delta^15, "N (\u2030)",sep="")))+
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=28), axis.title=element_text(size=32,face="bold"))+
  theme(legend.position="none")
ggsave(filename="Appendix1_AllSpWithRagfish_Isotopes.png", plot=Ap1_AllSpRagfish, dpi=500, width=17, height=12, units="in")


#Add humboldt squid for biplot:
#All.Prey.Sum3<-read.table(here::here('data/PmSources-IaDg.csv'),sep=",",header=TRUE)
#ggplot(All.Prey.Sum3, aes(Mean.d13C, Mean.d15N, color=Species, label=Species)) + geom_point(size=3)+
  #geom_errorbarh(aes(xmax=All.Prey.Sum3$Mean.d13C+All.Prey.Sum3$SD.d13C,xmin=All.Prey.Sum3$Mean.d13C-All.Prey.Sum3$SD.d13C, height = 0.01)) +
 # geom_errorbar(aes(ymax=All.Prey.Sum3$Mean.d15N+All.Prey.Sum3$SD.d15N,ymin=All.Prey.Sum3$Mean.d15N-All.Prey.Sum3$SD.d15N, width = 0.01))+
  #geom_text(color="black",hjust=-0.05,vjust = -0.7, size = 5)+
  #xlab(expression(paste(delta^13, "C (\u2030)")))+
  #ylab(expression(paste(delta^15, "N (\u2030)")))+
  #scale_color_viridis_d()+
  #theme_bw(base_size = 24, base_family = "Helvetica")+
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #theme(legend.position="none")

####################################################################
# SIMM Data Appendix 1: 
# Mixing model results including ragfish in models
# I ran this with sablefish, dogfish & ragfish combined. 
#      [because of similarity in isotopic values/isotope space]

mixsiar.dir <- find.package("MixSIAR")
#paste0(mixsiar.dir,"/example_scripts")

#Sperm whale predator mixture data is the same: 
PmSimmData2<-read.table(here::here("data/PmSimmData2.csv"),sep=",",header=TRUE)

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


#Source data; first for no combination of sablefish & dogfish
sourcePm.Ia.Comb <- load_source_data(filename="data/Appendix1/PmSources-Ia-SaAfIaComb.csv",
                                  source_factors=NULL,
                                  conc_dep=FALSE,
                                  data_type="means",
                                  mixPm)   #Source data frame needs to be in format: "Species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", "n" (sample size each source)
sourcePm.Ia.TempComb <- load_source_data(filename="data/Appendix1/PmSources-Ia-SaAfIaComb.csv",
                                      source_factors=NULL,
                                      conc_dep=FALSE,
                                      data_type="means",
                                      mixPmTemporal)
sourcePm.Ia.SeasonComb <- load_source_data(filename="data/Appendix1/PmSources-Ia-SaAfIaComb.csv",
                                        source_factors=NULL,
                                        conc_dep=FALSE,
                                        data_type="means",
                                        mixPmSeasonal)
sourcePm.Ia.FreqNFreqComb <- load_source_data(filename="data/Appendix1/PmSources-Ia-SaAfIaComb.csv",
                                           source_factors=NULL,
                                           conc_dep=FALSE,
                                           data_type="means",
                                           mixPmFreqNFreq)
discrPm.Ia.Comb <- load_discr_data(filename="data/Appendix1/PmSourceTEFsd-Ia-AfIaSaComb.csv", mixPm)
discrPm.Ia.TempComb <- load_discr_data(filename="data/Appendix1/PmSourceTEFsd-Ia-AfIaSaComb.csv", mixPmTemporal)
discrPm.Ia.SeasonComb <- load_discr_data(filename="data/Appendix1/PmSourceTEFsd-Ia-AfIaSaComb.csv", mixPmSeasonal)
discrPm.Ia.FreqNFreqComb <- load_discr_data(filename="data/Appendix1/PmSourceTEFsd-Ia-AfIaSaComb.csv", mixPmFreqNFreq)

# Define model structure and write JAGS model file
resid_err <- TRUE
process_err <- TRUE
model_Pm_Ia_Comb <- "MixSIAR_model_Pm_Ia_Comb_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Pm_Ia_Comb, resid_err, process_err, mixPm, sourcePm.Ia.Comb)

model_Pm_Ia_TempComb <- "MixSIAR_model_Pm_Ia_TempComb_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_Ia_TempComb, resid_err, process_err, mixPmTemporal, sourcePm.Ia.TempComb)

model_Pm_Ia_SeasonComb <- "MixSIAR_model_Pm_Ia_SeasonComb_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_Ia_SeasonComb, resid_err, process_err, mixPmSeasonal, sourcePm.Ia.SeasonComb)

model_Pm_Ia_FreqNFreqComb <- "MixSIAR_model_Pm_Ia_FreqNFreqComb_uninf.txt" # Name of JAGS model file
write_JAGS_model(model_Pm_Ia_FreqNFreqComb, resid_err, process_err, mixPmFreqNFreq, sourcePm.Ia.FreqNFreqComb)

# Run the JAGS model ("normal" is working well ... about 35min to run, with chain length 100,000; burn in 50,000, and thin 50; with 3 chains; it seems to converge nicely) 
jags.uninf.Pm.Ia.Comb <- run_model(run="normal",mixPm,sourcePm.Ia.Comb,discrPm.Ia.Comb,model_Pm_Ia_Comb,alpha.prior = 1, resid_err, process_err) #took 25 min to run.
jags.uninf.Pm.Ia.TempComb <- run_model(run="normal",mixPmTemporal,sourcePm.Ia.TempComb,discrPm.Ia.TempComb, model_Pm_Ia_TempComb,alpha.prior = 1, resid_err, process_err) #'normal' took 25min 
jags.uninf.Pm.Ia.SeasonComb <- run_model(run="normal",mixPmSeasonal,sourcePm.Ia.SeasonComb,discrPm.Ia.SeasonComb,model_Pm_Ia_SeasonComb,alpha.prior = 1, resid_err, process_err) #Took 35min
jags.uninf.Pm.Ia.FreqNFreqComb <- run_model(run="normal",mixPmFreqNFreq,sourcePm.Ia.FreqNFreqComb, discrPm.Ia.FreqNFreqComb,model_Pm_Ia_FreqNFreqComb,alpha.prior = 1, resid_err, process_err) #Took 35min to run

# get posterior medians for new source groupings
#apply(combined$post, 2, median)
#summary_stat(jags.uninf.Pm.Comb, meanSD=FALSE, quantiles=c(.025,.5,.975), savetxt=FALSE)

# Process diagnostics, summary stats, and posterior plots
jags_output_fullmodel_Ia<- output_JAGS(jags.uninf.Pm.Ia.Comb, mixPm, sourcePm.Ia.Comb) #DIC=119.4686
jags_output_tempcomb_Ia<- output_JAGS(jags.uninf.Pm.Ia.TempComb, mixPmTemporal, sourcePm.Ia.TempComb) #DIC 118.4287 - model doesn't converge super well, might need to run more iterations! 
jags_output_seasoncomb_Ia<- output_JAGS(jags.uninf.Pm.Ia.SeasonComb, mixPmSeasonal, sourcePm.Ia.SeasonComb) #DIC 120.8179 - model did converge pretty well,
jags_output_fnfcomb_Ia<- output_JAGS(jags.uninf.Pm.Ia.FreqNFreqComb, mixPmFreqNFreq, sourcePm.Ia.FreqNFreqComb) #DIC 114.8703 - model converged pretty well

which(levels(mixPmTemporal$data$Temporal)=="Old") # [1] 1
which(levels(mixPmTemporal$data$Temporal)=="Recent") # [1] 2
which(levels(mixPmFreqNFreq$data$Frequent)=="Frequent") #[1] 1
which(levels(mixPmFreqNFreq$data$Frequent)=="Non-Frequent") #[1] 2
which(levels(mixPmSeasonal$data$Season)=="Early") #[1] 1
which(levels(mixPmSeasonal$data$Season)=="Mid") #[1] 3
which(levels(mixPmSeasonal$data$Season)=="Late") #[1] 2
which(sourcePm.Ia.TempComb$source_names == "Clubhook Squid") # [1] 1
which(sourcePm.Ia.TempComb$source_names == "Grenadier") # [1] 2
which(sourcePm.Ia.TempComb$source_names == "Magister Squid") # [1] 3
which(sourcePm.Ia.TempComb$source_names == "Sablefish.Dogfish.Ragfish") #[1] 4
which(sourcePm.Ia.TempComb$source_names == "Shortraker Rockfish") #[1] 5
which(sourcePm.Ia.TempComb$source_names == "Skate") #[1] 6
which(sourcePm.Ia.SeasonComb$source_names == "Clubhook Squid") # [1] 1
which(sourcePm.Ia.SeasonComb$source_names == "Grenadier") # [1] 2
which(sourcePm.Ia.SeasonComb$source_names == "Magister Squid") # [1] 3
which(sourcePm.Ia.SeasonComb$source_names == "Sablefish.Dogfish.Ragfish") #[1] 4
which(sourcePm.Ia.SeasonComb$source_names == "Shortraker Rockfish") #[1] 5
which(sourcePm.Ia.SeasonComb$source_names == "Skate") #[1] 6

#Prey List: 1=Clubhook squid, 2=Grenadier, 3=Magister Squid, 4=Sablefish Group, 5=Shortraker rockfish, 6=Skate
#Predator Factor List: TempComb (1=Old, 2=Recent); SeasonComb (1=Early, 2=Late, 3=Mid); FreqNFreqComb (1=Freq, 2=Non-Freq, 3=Unk)
AllCombIaSummary<-jags.uninf.Pm.Ia.Comb$BUGSoutput$summary
write.table(AllCombIaSummary, file="AllCombIaSum.csv",sep=',')
TempCombIaSummary<-jags.uninf.Pm.Ia.TempComb$BUGSoutput$summary
TempCombIaSummary2<-TempCombIaSummary[c("p.fac1[1,1]","p.fac1[2,1]","p.fac1[1,2]","p.fac1[2,2]","p.fac1[1,3]","p.fac1[2,3]","p.fac1[1,4]","p.fac1[2,4]","p.fac1[1,5]","p.fac1[2,5]","p.fac1[1,6]","p.fac1[2,6]"),]
write.table(TempCombIaSummary2, file="TempCombIaSum.csv", sep=',')
SeasonCombIaSummary<-jags.uninf.Pm.Ia.SeasonComb$BUGSoutput$summary
SeasonCombIaSummary2<-SeasonCombIaSummary[c("p.fac1[1,1]","p.fac1[2,1]","p.fac1[3,1]","p.fac1[1,2]","p.fac1[2,2]","p.fac1[3,2]","p.fac1[1,3]","p.fac1[2,3]","p.fac1[3,3]","p.fac1[1,4]","p.fac1[2,4]","p.fac1[3,4]","p.fac1[1,5]","p.fac1[2,5]","p.fac1[3,5]","p.fac1[1,6]","p.fac1[2,6]","p.fac1[3,6]"),]
write.table(SeasonCombIaSummary2, file="SeasonCombIaSum.csv", sep=',')
FNFCombIaSummary<-jags.uninf.Pm.Ia.FreqNFreqComb$BUGSoutput$summary
FNFCombIaSummary2<-FNFCombIaSummary[c("p.fac1[1,1]","p.fac1[2,1]","p.fac1[3,1]","p.fac1[1,2]","p.fac1[2,2]","p.fac1[3,2]","p.fac1[1,3]","p.fac1[2,3]","p.fac1[3,3]","p.fac1[1,4]","p.fac1[2,4]","p.fac1[3,4]","p.fac1[1,5]","p.fac1[2,5]","p.fac1[3,5]","p.fac1[1,6]","p.fac1[2,6]","p.fac1[3,6]"),]
write.table(FNFCombIaSummary2, file="FNFCombIaSum2.csv", sep=',')


TempCombIaSumBox<- read.table('data/Appendix1/TempCombIaSumBox.csv',sep=",",header=TRUE)
TempCombIaSumBox$Middle<-as.numeric(TempCombIaSumBox$Middle)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')
A<-ggplot(TempCombIaSumBox, aes(x, Middle)) +
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

SeasonCombIaSumBox<- read.table('data/Appendix1/SeasonCombIaSumBox.csv',sep=",",header=TRUE)
SeasonCombIaSumBox$Group<-factor(SeasonCombIaSumBox$Group, levels=c("Early","Mid","Late"), ordered=TRUE)
B<-ggplot(SeasonCombIaSumBox, aes(x, Middle)) +
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

FNFCombIaSumBox<- read.table('data/Appendix1/FNFCombIaSumBox.csv',sep=",",header=TRUE)
library(stringr)

C<- ggplot(FNFCombIaSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax, color=Group,fill=Group), stat="identity", color="black") +
  scale_fill_manual(values=c("grey",'white'))+
  xlab("Species") +
  ylab("Proportion of Diet")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=12),axis.text.y=element_text(size=12), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))+
  theme(legend.justification=c(0.8,0), legend.position=c(0.24,0.6))+
  theme(legend.box.background=element_rect(color="black"))+
  theme(legend.title=element_blank(), legend.text=element_text(size=11))+
  theme(legend.key.size=unit(0.8,'cm'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
C

Fig2_Appendix1<-ggarrange(A,B,C, labels = c("A", "B", "C"), ncol=1, nrow=3)
tiff(filename="Fig2_Appendix1_boxplots.tiff", height = 24, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
Fig2_Appendix1
dev.off()

AllCombIaSumBox<-read.table('data/Appendix1/AllCombIaSumBox.csv', sep=",",header=TRUE)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')
AllCombIa<-ggplot(AllCombIaSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity", fill="grey", color="black") +
  #scale_fill_manual(values=c("grey"))+
  xlab("Species") +
  ylab("Proportion of Diet")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_bw()+
  theme(axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), axis.title.x=element_text(size=28), axis.title.y=element_text(size=28))+
  #theme(legend.justification=c(0.8,0), legend.position=c(0.18,0.6))+
  #theme(legend.box.background=element_rect(color="black"))+
  #theme(legend.title=element_blank(), legend.text=element_text(size=11))+
  #theme(legend.key.size=unit(0.8,'cm'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
AllCombIa
ggsave(filename="Appendix1_AllSpWithRagfish_Boxplot.png", plot=AllCombIa, dpi=500, width=17, height=12, units="in")


#################################################################################
# Mixing model results separating sablefish and dogfish
#
# Mixtures (sperm whales) are the same, just need to re-organize sources

sourcePm.All7.Comb <- load_source_data(filename="data/Appendix1/PmSources.csv",
                                     source_factors=NULL,
                                     conc_dep=FALSE,
                                     data_type="means",
                                     mixPm)   #Source data frame needs to be in format: "Species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", "n" (sample size each source)

discrPm.All7.Comb <- load_discr_data(filename="data/Appendix1/PmSourceTEF.csv", mixPm)

# Define model structure and write JAGS model file
resid_err <- TRUE
process_err <- TRUE
model_Pm_All7_Comb <- "MixSIAR_model_Pm_All7_Comb_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Pm_All7_Comb, resid_err, process_err, mixPm, sourcePm.All7.Comb)

# Run the JAGS model ("normal"; chain length 100,000; burn in 50,000, and thin 50; with 3 chains; it seems to converge nicely) 
jags.uninf.Pm.All7.Comb <- run_model(run="normal",mixPm,sourcePm.All7.Comb,discrPm.All7.Comb,model_Pm_All7_Comb,alpha.prior = 1, resid_err, process_err) #took 25 min to run.

which(sourcePm.All7.Comb$source_names == "Clubhook Squid") # [1] 1
which(sourcePm.All7.Comb$source_names == "Grenadier") # [1] 2
which(sourcePm.All7.Comb$source_names == "Magister Squid") # [1] 3
which(sourcePm.All7.Comb$source_names == "Sablefish") #[1] 4
which(sourcePm.All7.Comb$source_names == "Spiny Dogfish") #[1] 7
which(sourcePm.All7.Comb$source_names == "Shortraker Rockfish") #[1] 5
which(sourcePm.All7.Comb$source_names == "Skate") #[1] 6
#Prey List: 1=Clubhook squid, 2=Grenadier, 3=Magister Squid, 4=Sablefish, 5=Shortraker rockfish, 6=Skate, 7=Dogfish

AllCombAll7Summary<-jags.uninf.Pm.All7.Comb$BUGSoutput$summary
write.table(AllCombAll7Summary, file="AllCombAll7Sum.csv",sep=',')

AllCombAll7SumBox<-read.table('data/Appendix1/AllCombAll7SumBox.csv', sep=",",header=TRUE)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')

All7<-ggplot(AllCombAll7SumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity", fill="grey", color="black") +
  xlab("Species") +
  ylab("Proportion of Diet")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_bw()+
  theme(axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), axis.title.x=element_text(size=28), axis.title.y=element_text(size=28))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
All7

ggsave(filename="Appendix1_AllSp_AfSa_separate_Boxplot.png", plot=All7, dpi=500, width=17, height=12, units="in")


Fig2_Appendix1<-ggarrange(A,B,C, labels = c("A", "B", "C"), ncol=1, nrow=3)
tiff(filename="Fig2_Appendix1_boxplots.tiff", height = 24, width = 17, units = 'cm', 
     compression = "lzw", res = 200)
Fig2_Appendix1
dev.off()

##########################################################################################
# Mixing model results for low TEF & high TEF
# 
# Mixtures (sperm whales) are the same, just need to re-organize source TEFs
#
# #Load in source data for MixSIAR:
sourcePm.LowHighTEF.Comb <- load_source_data(filename="data/PmSources_SaAf_Comb.csv",
                                  source_factors=NULL,
                                  conc_dep=FALSE,
                                  data_type="means",
                                  mixPm)   #Source data frame needs to be in format: "Species", "Meand13C", "SDd13C", "Meand15N", "SDd15N", "n" (sample size each source)


discrPm.LowTEF.Comb <- load_discr_data(filename="data/Appendix1/PmSourceTEF1.6-AfSaComb.csv", mixPm)
discrPm.HighTEF.Comb <- load_discr_data(filename="data/Appendix1/PmSourceTEF2.8-AfSaComb.csv", mixPm)

model_Pm_LowTEF_Comb <- "MixSIAR_model_Pm_LowTEF_Comb_uninf.txt"   # Name of the JAGS model file
write_JAGS_model(model_Pm_LowTEF_Comb, resid_err, process_err, mixPm, sourcePm.LowHighTEF.Comb)
model_Pm_HighTEF_Comb <- "MixSIAR_model_Pm_HighTEF_Comb_uninf.txt"
write_JAGS_model(model_Pm_HighTEF_Comb, resid_err, process_err, mixPm, sourcePm.LowHighTEF.Comb)

# Run the JAGS model
jags.uninf.Pm.LowTEF.Comb <- run_model(run="normal",mixPm,sourcePm.LowHighTEF.Comb,discrPm.LowTEF.Comb,model_Pm_LowTEF_Comb,alpha.prior = 1, resid_err, process_err) #took 25 min to run
jags.uninf.Pm.HighTEF.Comb <- run_model(run="normal",mixPm,sourcePm.LowHighTEF.Comb,discrPm.HighTEF.Comb,model_Pm_HighTEF_Comb,alpha.prior = 1, resid_err, process_err) #took 25 min to run

which(sourcePm.LowHighTEF.Comb$source_names == "Clubhook Squid") # [1] 1
which(sourcePm.LowHighTEF.Comb$source_names == "Grenadier") # [1] 2
which(sourcePm.LowHighTEF.Comb$source_names == "Magister Squid") # [1] 3
which(sourcePm.LowHighTEF.Comb$source_names == "Sablefish.Dogfish") #[1] 4
which(sourcePm.LowHighTEF.Comb$source_names == "Shortraker Rockfish") #[1] 5
which(sourcePm.LowHighTEF.Comb$source_names == "Skate") #[1] 6
#Prey List: 1=Clubhook squid, 2=Grenadier, 3=Magister Squid, 4=Sablefish.Dogfish, 5=Shortraker rockfish, 6=Skate

AllCombLowTEFSummary<-jags.uninf.Pm.LowTEF.Comb$BUGSoutput$summary
write.table(AllCombLowTEFSummary, file="AllCombLowTEFSum.csv",sep=',')
AllCombHighTEFSummary<-jags.uninf.Pm.HighTEF.Comb$BUGSoutput$summary
write.table(AllCombHighTEFSummary, file="AllCombHighTEFSum.csv",sep=',')

AllCombLowTEFSumBox<-read.table('data/Appendix1/AllCombLowTEFSumBox.csv', sep=",",header=TRUE)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')
LowTEF<-ggplot(AllCombLowTEFSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity", fill="grey", color="black") +
  xlab("Species") +
  ylab("Proportion of Diet")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_bw()+
  theme(axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), axis.title.x=element_text(size=28), axis.title.y=element_text(size=28))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
LowTEF
ggsave(filename="Appendix1_LowTEF_Boxplot.png", plot=LowTEF, dpi=500, width=17, height=12, units="in")

AllCombHighTEFSumBox<-read.table('data/Appendix1/AllCombHighTEFSumBox.csv', sep=",",header=TRUE)
pattern.type=c('nwlines',"blank") #also 'waves', 'hdashes', 'crosshatch', 'dots', 'grid', 'hlines', 'nelines', shells', 'circles1', 'circles2', 'vdashes', 'bricks'.
pattern.color=c('black','black')
background.color=c('white', 'gray80')
HighTEF<-ggplot(AllCombHighTEFSumBox, aes(x, Middle)) +
  geom_boxplot(aes(ymin=ymin, lower=Lower, middle=Middle, upper=Upper, ymax=ymax), stat="identity", fill="grey", color="black") +
  xlab("Species") +
  ylab("Proportion of Diet")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_bw()+
  theme(axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), axis.title.x=element_text(size=28), axis.title.y=element_text(size=28))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
HighTEF
ggsave(filename="Appendix1_HighTEF_Boxplot.png", plot=HighTEF, dpi=500, width=17, height=12, units="in")

##############################################################################################
library(simmr)
PmSimmData3<-as.matrix(PmSimmData2) #[,c(2,1)]
PmSources<-read.csv(here::here("data/PmSources_SaAf_Comb.csv"), sep=",",header=TRUE)
colnames(PmSources)<-c("Species","d13C.m", "d13C.sd","d15N.m", "d15N.sd", "n")
PmSources<-PmSources[,c(1,2,3,4,5)]
SourceTEFLow<- read.csv(here:here("data/Appendix1/PmSourceTEF1.6-AfSaComb.csv"), sep=",", header=TRUE)

simmr_1_comb = simmr_load(mixtures=PmSimmData3,  
                          source_names=as.character(PmSources[,1]),
                          source_means=PmSources[,c(2,4)],
                          source_sds=PmSources[,c(3,5)],
                          correction_means=SourceTEFsd[,c(2,4)],
                          correction_sds=SourceTEFsd[,c(3,5)])


