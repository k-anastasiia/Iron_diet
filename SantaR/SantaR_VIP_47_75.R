#install.packages("santaR")
# Install devtools
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

setwd('/home/anastasiia/Desktop/mouse_studies/SantaR/SantaR_VIP')
getwd() #to see the working directory 
metadata<- read_csv("metadata_full (ak).csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/Day_75/SantaR_VIP/CLR_Scaled.csv")
colnames(scaled_table)[1] <- "row_ID"

#Day 75, Normal-Deficient
numbers_to_match_75_ND <- c('row_ID','8823_', '2352_', '8818_', '8454_', '2768_', '2299_', '2057_', '4317_', '2670_', '2523_', '8707_','8721_', '8458_',
'3784_', '8471_', '8504_', '6101_', '8714_', '1841_', '1692_', '2254_', '8748_', '1784_', '8511_', '8522_', '2202_', '3302_', '8782_', '5741_', '1819_')
scaled_table_75_ND<- scaled_table %>% select(matches(paste(numbers_to_match_75_ND, collapse = "|")))
excluded_columns <- c(15, 18, 24, 27, 30,32,33)
scaled_table_75_ND<- scaled_table_75_ND[, -excluded_columns]
data_merge_santaR_75_ND <- merge(scaled_table_75_ND, metadata, by.x = "row_ID", by.y = "...1",  
                                  all.x = FALSE, all.y = FALSE)

inputData     <- data.frame(data_merge_santaR_75_ND[,2:31])
ind           <- data_merge_santaR_75_ND[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR_75_ND[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR_75_ND[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

p1 <- santaR_plot(SANTAObj[[1]])
p2 <- santaR_plot(SANTAObj[[2]])
p3 <- santaR_plot(SANTAObj[[3]])
p4 <- santaR_plot(SANTAObj[[4]])
p5 <- santaR_plot(SANTAObj[[5]])
p6 <- santaR_plot(SANTAObj[[6]])
p7 <- santaR_plot(SANTAObj[[7]])
p8 <- santaR_plot(SANTAObj[[8]])
p9 <- santaR_plot(SANTAObj[[9]])
p10 <- santaR_plot(SANTAObj[[10]])
p11 <- santaR_plot(SANTAObj[[11]])
p12 <- santaR_plot(SANTAObj[[12]])
p13 <- santaR_plot(SANTAObj[[13]])
p14 <- santaR_plot(SANTAObj[[14]])
p15 <- santaR_plot(SANTAObj[[15]])
p16 <- santaR_plot(SANTAObj[[16]])
p17 <- santaR_plot(SANTAObj[[17]])
p18 <- santaR_plot(SANTAObj[[18]])
p19 <- santaR_plot(SANTAObj[[19]])
p20 <- santaR_plot(SANTAObj[[20]])
p21 <- santaR_plot(SANTAObj[[21]])
p22 <- santaR_plot(SANTAObj[[22]])
p23 <- santaR_plot(SANTAObj[[23]])
p24 <- santaR_plot(SANTAObj[[24]])
p25 <- santaR_plot(SANTAObj[[25]])
p26 <- santaR_plot(SANTAObj[[26]])
p27 <- santaR_plot(SANTAObj[[27]])
p28 <- santaR_plot(SANTAObj[[28]])
p29 <- santaR_plot(SANTAObj[[29]])
p30 <- santaR_plot(SANTAObj[[30]])


data_merge_santaR_75_ND$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_75_ND$ATTRIBUTE_Study_Day)
p23 <- p23 + scale_x_continuous(breaks = data_merge_santaR_75_ND$ATTRIBUTE_Study_Day)
p23 <- p23 + theme(panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())
p23

grid.arrange(p1, p2, p3, p4)
grid.arrange(p5, p6, p7, p8)
grid.arrange(p9, p10, p11, p12)
grid.arrange(p13, p14, p15, p16)
grid.arrange(p17, p18, p19, p20)
grid.arrange(p21, p22, p23, p24)
grid.arrange(p25, p26, p27, p28)
grid.arrange(p29, p30)


#Day 75, Normal-0.6%Fe
numbers_to_match_75_N06 <- c('row_ID', '7395_', '6857_', '13978_', '7347_', '5519_', '4660_', '6974_', '13093_', '12251_', '13890_', '5508_', '7997_', '2352_', '6018_', '12537_', '7364_',
                             '12337_', '11281_', '7664_', '11329_', '3810_', '10043_', '8775_', '13679_', '7372_', '1413_', '2768_', '13751_', '4471_', '5509_')
scaled_table_75_N06<- scaled_table %>%   select(matches(paste(numbers_to_match_75_N06, collapse = "|")))
excluded_columns <- c(16,25,27,30)
scaled_table_75_N06<- scaled_table_75_N06[, -excluded_columns]
data_merge_santaR_75_N06 <- merge(scaled_table_75_N06, metadata, by.x = "row_ID", by.y = "...1",  
                                  all.x = FALSE, all.y = FALSE)
inputData     <- data.frame(data_merge_santaR_75_N06[,2:31])
ind           <- data_merge_santaR_75_N06[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR_75_N06[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR_75_N06[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

p1 <- santaR_plot(SANTAObj[[1]])
p2 <- santaR_plot(SANTAObj[[2]])
p3 <- santaR_plot(SANTAObj[[3]])
p4 <- santaR_plot(SANTAObj[[4]])
p5 <- santaR_plot(SANTAObj[[5]])
p6 <- santaR_plot(SANTAObj[[6]])
p7 <- santaR_plot(SANTAObj[[7]])
p8 <- santaR_plot(SANTAObj[[8]])
p9 <- santaR_plot(SANTAObj[[9]])
p10 <- santaR_plot(SANTAObj[[10]])
p11 <- santaR_plot(SANTAObj[[11]])
p12 <- santaR_plot(SANTAObj[[12]])
p13 <- santaR_plot(SANTAObj[[13]])
p14 <- santaR_plot(SANTAObj[[14]])
p15 <- santaR_plot(SANTAObj[[15]])
p16 <- santaR_plot(SANTAObj[[16]])
p17 <- santaR_plot(SANTAObj[[17]])
p18 <- santaR_plot(SANTAObj[[18]])
p19 <- santaR_plot(SANTAObj[[19]])
p20 <- santaR_plot(SANTAObj[[20]])
p21 <- santaR_plot(SANTAObj[[21]])
p22 <- santaR_plot(SANTAObj[[22]])
p23 <- santaR_plot(SANTAObj[[23]])
p24 <- santaR_plot(SANTAObj[[24]])
p25 <- santaR_plot(SANTAObj[[25]])
p26 <- santaR_plot(SANTAObj[[26]])
p27 <- santaR_plot(SANTAObj[[27]])
p28 <- santaR_plot(SANTAObj[[28]])
p29 <- santaR_plot(SANTAObj[[29]])
p30 <- santaR_plot(SANTAObj[[30]])

grid.arrange(p1, p2, p3, p4)
grid.arrange(p5, p6, p7, p8)
grid.arrange(p9, p10, p11, p12)
grid.arrange(p13, p14, p15, p16)
grid.arrange(p17, p18, p19, p20)
grid.arrange(p21, p22, p23, p24)
grid.arrange(p25, p26, p27, p28)
grid.arrange(p29, p30)

grid.arrange(p1, p2, p3, p4)
grid.arrange(p5, p9, p12, p10)
grid.arrange(p13, p22, p27,p20, p28)

















#Day 47, Normal-0.6%Fe
numbers_to_match_47_N06 <- c('row_ID','7395_', '4660_', '1413_', '7664_', '14240_', '6857_', '6871_', '12251_', '11281_', '8559_', '9442_', '7083_', '2254_', '6018_',
                             '5352_', '5680_', '7997_', '7018_', '8868_', '10980_', '11829_', '7373_', '7655_', '6142_', '4194_', '8401_', '1014_', '10043_', '5519_', '414_')
scaled_table_47_N06<- scaled_table %>%   select(matches(paste(numbers_to_match_47_N06, collapse = "|")))
excluded_columns <- c(12,16,26,28,35,37)
scaled_table_47_N06<- scaled_table_47_N06[, -excluded_columns]
data_merge_santaR_47_N06 <- merge(scaled_table_47_N06, metadata, by.x = "row_ID", by.y = "...1",  
                                  all.x = FALSE, all.y = FALSE)

inputData     <- data.frame(data_merge_santaR_47_N06[,2:31])
ind           <- data_merge_santaR_47_N06[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR_47_N06[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR_47_N06[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

p1 <- santaR_plot(SANTAObj[[1]])
p2 <- santaR_plot(SANTAObj[[2]])
p3 <- santaR_plot(SANTAObj[[3]])
p4 <- santaR_plot(SANTAObj[[4]])
p5 <- santaR_plot(SANTAObj[[5]])
p6 <- santaR_plot(SANTAObj[[6]])
p7 <- santaR_plot(SANTAObj[[7]])
p8 <- santaR_plot(SANTAObj[[8]])
p9 <- santaR_plot(SANTAObj[[9]])
p10 <- santaR_plot(SANTAObj[[10]])
p11 <- santaR_plot(SANTAObj[[11]])
p12 <- santaR_plot(SANTAObj[[12]])
p13 <- santaR_plot(SANTAObj[[13]])
p14 <- santaR_plot(SANTAObj[[14]])
p15 <- santaR_plot(SANTAObj[[15]])
p16 <- santaR_plot(SANTAObj[[16]])
p17 <- santaR_plot(SANTAObj[[17]])
p18 <- santaR_plot(SANTAObj[[18]])
p19 <- santaR_plot(SANTAObj[[19]])
p20 <- santaR_plot(SANTAObj[[20]])
p21 <- santaR_plot(SANTAObj[[21]])
p22 <- santaR_plot(SANTAObj[[22]])
p23 <- santaR_plot(SANTAObj[[23]])
p24 <- santaR_plot(SANTAObj[[24]])
p25 <- santaR_plot(SANTAObj[[25]])
p26 <- santaR_plot(SANTAObj[[26]])
p27 <- santaR_plot(SANTAObj[[27]])
p28 <- santaR_plot(SANTAObj[[28]])
p29 <- santaR_plot(SANTAObj[[29]])
p30 <- santaR_plot(SANTAObj[[30]])

grid.arrange(p1, p2, p3, p4)
grid.arrange(p5, p6, p7, p8)
grid.arrange(p9, p10, p11, p12)
grid.arrange(p13, p14, p15, p16)
grid.arrange(p17, p18, p19, p20)
grid.arrange(p21, p22, p23, p24)
grid.arrange(p25, p26, p27, p28)
grid.arrange(p29, p30)





#6857
data_merge_santaR_47_N06$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p19 <- p19 + scale_x_continuous(breaks = data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p19 <- p19 + theme(panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank())
p19


#11829
data_merge_santaR_47_N06$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p28 <- p28 + scale_x_continuous(breaks = data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p28 <- p28 + theme(panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())
p28


data_merge_santaR_47_N06$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p6 <- p6 + scale_x_continuous(breaks = data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p6 <- p6 + theme(panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank())
p6

p3 <- p3 + scale_x_continuous(breaks = data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p3 <- p3 + theme(panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank())
p3



p9 <- p9 + scale_x_continuous(breaks = data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p9 <- p9 + theme(panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank())
p9




p30 <- p30 + scale_x_continuous(breaks = data_merge_santaR_47_N06$ATTRIBUTE_Study_Day)
p30 <- p30 + theme(panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())
p30






#Day 47, Normal-Deficient
numbers_to_match_47_ND <- c('row_ID','3784_', '2670_', '1807_', '1692_', '6516_', '1163_', '3767_', '3954_', '2300_', '3760_', '5823_', '6661_', '6638_', '4083_', '3548_', '4081_', '4772_', '5536_',
                            '5439_', '6528_', '6817_', '4771_', '5526_', '2679_', '5393_', '2402_', '6824_', '5642_', '4817_', '5670_')
scaled_table_47_ND<- scaled_table %>%   select(matches(paste(numbers_to_match_47_ND, collapse = "|")))
excluded_columns <- c(8,10,11,16,18,22)
scaled_table_47_ND<- scaled_table_47_ND[, -excluded_columns]
data_merge_santaR_47_ND <- merge(scaled_table_47_ND, metadata, by.x = "row_ID", by.y = "...1",  
                                 all.x = FALSE, all.y = FALSE)

inputData     <- data.frame(data_merge_santaR_47_ND[,2:31])
ind           <- data_merge_santaR_47_N06[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR_47_N06[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR_47_N06[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

p1 <- santaR_plot(SANTAObj[[1]])
p2 <- santaR_plot(SANTAObj[[2]])
p3 <- santaR_plot(SANTAObj[[3]])
p4 <- santaR_plot(SANTAObj[[4]])
p5 <- santaR_plot(SANTAObj[[5]])
p6 <- santaR_plot(SANTAObj[[6]])
p7 <- santaR_plot(SANTAObj[[7]])
p8 <- santaR_plot(SANTAObj[[8]])
p9 <- santaR_plot(SANTAObj[[9]])
p10 <- santaR_plot(SANTAObj[[10]])
p11 <- santaR_plot(SANTAObj[[11]])
p12 <- santaR_plot(SANTAObj[[12]])
p13 <- santaR_plot(SANTAObj[[13]])
p14 <- santaR_plot(SANTAObj[[14]])
p15 <- santaR_plot(SANTAObj[[15]])
p16 <- santaR_plot(SANTAObj[[16]])
p17 <- santaR_plot(SANTAObj[[17]])
p18 <- santaR_plot(SANTAObj[[18]])
p19 <- santaR_plot(SANTAObj[[19]])
p20 <- santaR_plot(SANTAObj[[20]])
p21 <- santaR_plot(SANTAObj[[21]])
p22 <- santaR_plot(SANTAObj[[22]])
p23 <- santaR_plot(SANTAObj[[23]])
p24 <- santaR_plot(SANTAObj[[24]])
p25 <- santaR_plot(SANTAObj[[25]])
p26 <- santaR_plot(SANTAObj[[26]])
p27 <- santaR_plot(SANTAObj[[27]])
p28 <- santaR_plot(SANTAObj[[28]])
p29 <- santaR_plot(SANTAObj[[29]])
p30 <- santaR_plot(SANTAObj[[30]])

grid.arrange(p1, p2, p3, p4)
grid.arrange(p5, p6, p7, p8)
grid.arrange(p9, p10, p11, p12)
grid.arrange(p13, p14, p15, p16)
grid.arrange(p17, p18, p19, p20)
grid.arrange(p21, p22, p23, p24)
grid.arrange(p25, p26, p27, p28)
grid.arrange(p29, p30)



grid.arrange(p9,p13, p14)


#peptides+cytoscape
grid.arrange(p13, p14, p19,p20,p21, p22, p28, p29)

data_merge_santaR_47_ND$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_47_ND$ATTRIBUTE_Study_Day)

p28 <- p28 + scale_x_continuous(breaks = data_merge_santaR_47_ND$ATTRIBUTE_Study_Day)
p28 <- p28 + theme(panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank())
p28


data_merge_santaR_47_ND$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_47_ND$ATTRIBUTE_Study_Day)
#3784
p15 <- p15 + scale_x_continuous(breaks = data_merge_santaR_47_ND$ATTRIBUTE_Study_Day)
p15 <- p15 + theme(panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())
p15

#6638
p19 <- p19 + scale_x_continuous(breaks = data_merge_santaR_47_ND$ATTRIBUTE_Study_Day)
p19 <- p19 + theme(panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank())
p19
