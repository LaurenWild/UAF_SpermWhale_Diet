# Bathymetry Mapping - Code From Jordan Watson, NOAA & Cheryl Barnes, UAF; 
# 8-30-2018 & 4-16-2019
# IMPORTANT!!!!! Load the tidyverse library first. 
# I have received a funky error message from the map_data() function 
# if I load the "maps" library prior to the "tidyverse" (or "ggplot2") library.
install.packages('tidyverse')
library(tidyverse) # Contains the ggplot2 and dplyr packages
install.packages('PBSmapping')
library(PBSmapping) # From which we'll use the clipPolys function
install.packages('marmap')
library(marmap) # From which we can get bathymetry
install.packages('maps')
library(maps) # From which we get the world2 basemap
install.packages('mapproj')
library(mapproj)

########################################################
# MAP OF EASTERN GOA WITH BIOPSY & PREY LOCATIONS #
########################################################
# Load map data and set new grid boundaries (DD):
library(PBSmapping)
data(nepacLLhigh)
lonmin = -150
lonmax = -133.5
latmin = 54
latmax = 62

# Clip maps to predetermined boundaries (this will take a few moments to complete):
world = fortify(nepacLLhigh)
world2 = clipPolys(world, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

setwd("/Users/laurenwild/Desktop/UAF/Thesis/MAPS/")
Canada = raster::getData("GADM", country = "CAN", level = 0)
Canada = fortify(Canada)
names(Canada) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Canada$PID = as.numeric(Canada$PID)
Canada = clipPolys(Canada, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

#to add datapoints: 
Biopsy<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/TissueSpreadsheets/Biopsy_Details_forR_20180829.csv',sep=",",header=TRUE)
View(Biopsy)
Prey<- read.table('/Users/laurenwild/Desktop/UAF/Thesis/StableIsotopes/Data/Prey.LE.final.Isotopes.forR5.csv',sep=",",header=TRUE)
View(Prey)
Prey$Longitude4<-Prey$Longitude3-360
Prey$Longitude4<-Prey$Longitude3 * -1

#Create the map:
Fig1 = ggplot() + 
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_point(data=Prey,aes(x=Longitude4,y=Latitude3),shape=17,size=2,color="red")+
  geom_point(data=PmFull, aes(x=Longitude, y=Latitude), size=1.5, col="blue")+
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.title.y = element_text(size=16)) +
  xlab("Longitude")+
  ylab("Latitude")

setwd('/Users/laurenwild/Desktop')
ggsave(filename="Fig1_Biopsy_Prey_Map.png", plot=Fig1, dpi=500, width=12, height=8, units="in")