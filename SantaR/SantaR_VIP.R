#install.packages("santaR")
if(!require("devtools")) install.packages("devtools")
devtools::install_github("adwolfer/santaR", ref="master")
library(ggpubr)
library(santaR)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(pheatmap)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(tidyverse)
library(gridExtra)

setwd('~/Desktop/GitHub/SantaR/VIP')
getwd()

scaled_table <- read_csv("CLR_Scaled.csv")
metadata <- read_csv("md_new.csv")
metadata <- metadata[-1:-10,-1] #exclude bedding and blanks
colnames(scaled_table)[1] <- "row_ID"


#Normal-Deficient
numbers_to_match_ND <- c('row_ID', '6817_', '7913_', '2673_', '3825_', '4813_', '3777_', 
                         '3034_', '6441_', '5670_', '6661_', '5183_', '4817_', '4083_', '3760_', '3954_', '2094_',
                         '2352_', '5823_', '2300_', '4081_', '2057_', '2254_', '3767_', '3336_', '1163_', '6528_', '6516_', 
                         '2299_', '2670_', '3784_')
scaled_table_ND<- scaled_table %>% select(matches(paste(numbers_to_match_ND, collapse = "|")))
data_merge_santaR_ND <- merge(scaled_table_ND, metadata, by.x = "row_ID", by.y = "filename",  
                                  all.x = FALSE, all.y = FALSE)

inputData     <- data.frame(data_merge_santaR_ND[,2:40])
ind           <- data_merge_santaR_ND[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR_ND[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR_ND[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

#3784-oligopeptide
santaR_plot(SANTAObj[[18]])

data_merge_santaR_ND$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_ND$ATTRIBUTE_Study_Day)
x3784 <-santaR_plot(SANTAObj[[18]], colorVect=c("blue","red","green"),  title = "Oligopeptide",xlab = "Days")
x3784 <- x3784  + scale_x_continuous(breaks = data_merge_santaR_ND$ATTRIBUTE_Study_Day)
x3784<- x3784 + theme(panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        axis.title.y = element_blank(),  
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(), 
                        axis.text.x = element_blank())  
x3784
x3784<-x3784+  theme(legend.position = "x3784")

#2299-N-acyl-alpha amino acids and derivatives
santaR_plot(SANTAObj[[3]])

data_merge_santaR_ND$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_ND$ATTRIBUTE_Study_Day)
x2299 <-santaR_plot(SANTAObj[[3]], colorVect=c("blue","red","green"),  title = "N-acyl-alpha amino acids and derivatives",xlab = "Days")
x2299 <- x2299  + scale_x_continuous(breaks = data_merge_santaR_ND$ATTRIBUTE_Study_Day)
x2299<- x2299 + theme(panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        axis.title.y = element_blank(),  
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(), 
                        axis.text.x = element_blank())  
x2299
x2299<-x2299+  theme(legend.position = "x2299")

#Normal-0.6%Fe
numbers_to_match_N06 <- c('row_ID','13115_', '15035_', '2254_', '7347_', '7373_', '3810_', '8401_', 
                          '13093_', '11829_', '6896_', '13890_', '10043_', '5680_', '6841_', '13978_', 
                          '6142_', '9040_', '5519_', '7664_', '7997_', '12251_', '7574_', '6857_', 
                          '6018_', '9305_', '12537_', '11281_', '5508_', '4660_', '7395_')
scaled_table_N06<- scaled_table %>%   select(matches(paste(numbers_to_match_N06, collapse = "|")))
data_merge_santaR_N06 <- merge(scaled_table_N06, metadata, by.x = "row_ID", by.y = "filename",  
                                  all.x = FALSE, all.y = FALSE)
inputData     <- data.frame(data_merge_santaR_N06[,2:32])
ind           <- data_merge_santaR_N06[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR_N06[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR_N06[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)


#13890-Acyclic diterpenoids
santaR_plot(SANTAObj[[21]])

data_merge_santaR_N06$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_N06$ATTRIBUTE_Study_Day)
x13890 <-santaR_plot(SANTAObj[[21]], colorVect=c("blue","red","green"),  title = "Acyclic diterpenoid",xlab = "Days")
x13890 <- x13890  + scale_x_continuous(breaks = data_merge_santaR_N06$ATTRIBUTE_Study_Day)
x13890<- x13890 + theme(panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        axis.title.y = element_blank(),  
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(), 
                        axis.text.x = element_blank())  
x13890
x138907<-x13890+  theme(legend.position = "x13890")

#6857-Fatty alcohols
 santaR_plot(SANTAObj[[18]])
 
 data_merge_santaR_N06$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_N06$ATTRIBUTE_Study_Day)
 x6857 <-santaR_plot(SANTAObj[[18]], colorVect=c("blue","red","green"),  title = "Fatty alcohol",xlab = "Days")
 x6857 <- x6857  + scale_x_continuous(breaks = data_merge_santaR_N06$ATTRIBUTE_Study_Day)
 x6857<- x6857 + theme(panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank(),
                           panel.grid.major.y = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           axis.title.y = element_blank(),  
                           axis.text.y = element_blank(),
                           axis.title.x = element_blank(), 
                           axis.text.x = element_blank())  
 x6857
 x6857<-x6857+  theme(legend.position = "x6857")
 


 
 
 
