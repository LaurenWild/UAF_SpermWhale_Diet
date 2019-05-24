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
