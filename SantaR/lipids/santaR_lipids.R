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

setwd('~/Desktop/GitHub/SantaR/lipids')
getwd()



#####All diets
lipids <- read_csv("lipids.csv")

inputDataM     <- data.frame(lipids[,2:8])
indM           <- lipids[["ATTIBUTE_mouse_no"]]
timeM          <- as.numeric(lipids[["ATTRIBUTE_Study_Day"]])
groupM         <- lipids[["ATTRIBUTE_Diet"]]
SANTAObjM  <- santaR_auto_fit(inputDataM, indM, timeM, groupM, df=5)

Met_C6 <-santaR_plot(SANTAObjM[[1]], colorVect=c("royalblue","green4",'red3'),  title = "Met-caproic acid 1",xlab = "Days")
Met_C6   <-Met_C6   + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
Met_C6   <- Met_C6  + theme(panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            axis.title.y = element_blank(),  # Remove y-axis label
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank()) 
Met_C6 <-Met_C6 + theme(legend.position = "none") 


Met_C6a <-santaR_plot(SANTAObjM[[2]], colorVect=c("royalblue","green4",'red3'),  title = "Met-caproic acid 2",xlab = "Days")
Met_C6a   <-Met_C6a   + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
Met_C6a   <- Met_C6a  + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank(),
                              axis.title.y = element_blank(),  # Remove y-axis label
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.text.x = element_blank()) 
Met_C6a <-Met_C6a + theme(legend.position = "none") 


Met_C6_all <-santaR_plot(SANTAObjM[[3]], colorVect=c("royalblue","green4",'red3'),  title = "Met-Caproic acid",xlab = "Days")
Met_C6_all   <-Met_C6_all   + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
Met_C6_all   <- Met_C6_all  + theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    axis.title.y = element_blank(),  # Remove y-axis label
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank()) 
Met_C6_all <-Met_C6_all + theme(legend.position = "none") 


PC_20_2 <-santaR_plot(SANTAObjM[[4]], colorVect=c("royalblue","green4",'red3'),  title = "PC(20:2/0:0) 1",xlab = "Days")
PC_20_2  <-PC_20_2  + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
PC_20_2  <- PC_20_2  + theme(panel.grid.major.x = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             axis.title.y = element_blank(),  # Remove y-axis label
                             axis.text.y = element_blank(),
                             axis.title.x = element_blank(),
                             axis.text.x = element_blank()) 
PC_20_2 <-PC_20_2 + theme(legend.position = "none") 


PC_20_2a <-santaR_plot(SANTAObjM[[5]], colorVect=c("royalblue","green4",'red3'),  title = "PC(20:2/0:0) 2",xlab = "Days")
PC_20_2a   <-PC_20_2a   + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
PC_20_2a   <- PC_20_2a  + theme(panel.grid.major.x = element_blank(),
                                panel.grid.minor.x = element_blank(),
                                panel.grid.major.y = element_blank(),
                                panel.grid.minor.y = element_blank(),
                                axis.title.y = element_blank(),  # Remove y-axis label
                                axis.text.y = element_blank(),
                                axis.title.x = element_blank(),
                                axis.text.x = element_blank()) 
PC_20_2a <-PC_20_2a + theme(legend.position = "none") 


PC_20_2_all <-santaR_plot(SANTAObjM[[6]], colorVect=c("royalblue","green4",'red3'),  title = "PC(20:2/0:0)",xlab = "Days")
PC_20_2_all   <-PC_20_2_all   + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
PC_20_2_all   <- PC_20_2_all  + theme(panel.grid.major.x = element_blank(),
                                      panel.grid.minor.x = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.minor.y = element_blank(),
                                      axis.title.y = element_blank(),  # Remove y-axis label
                                      axis.text.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.text.x = element_blank()) 
PC_20_2_all <-PC_20_2_all + theme(legend.position = "none") 



Stearidonic_acid <-santaR_plot(SANTAObjM[[7]], colorVect=c("royalblue","green4",'red3'),  title = "Stearidonic acid ",xlab = "Days")
Stearidonic_acid   <-Stearidonic_acid   + scale_x_continuous(breaks = lipids$ATTRIBUTE_Study_Day)
Stearidonic_acid   <- Stearidonic_acid   + theme(panel.grid.major.x = element_blank(),
                                                 panel.grid.minor.x = element_blank(),
                                                 panel.grid.major.y = element_blank(),
                                                 panel.grid.minor.y = element_blank(),
                                                 axis.title.y = element_blank(),  # Remove y-axis label
                                                 axis.text.y = element_blank(),
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()) 

Stearidonic_acid  <-Stearidonic_acid + theme(legend.position = "none") 

grid.arrange(PC_20_2_all,PC_20_2, PC_20_2a, ncol=3)
grid.arrange(Met_C6_all, Met_C6, Met_C6a, Stearidonic_acid,ncol=3)

########Normal-overload diets
lipids_N06 <- read_csv("lipids_N06.csv")

inputData     <- data.frame(lipids_N06[,2:8])
ind          <- lipids_N06[["ATTIBUTE_mouse_no"]]
time         <- as.numeric(lipids_N06[["ATTRIBUTE_Study_Day"]])
group        <- lipids_N06[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)


Met_C6_all <-santaR_plot(SANTAObj[[3]], colorVect=c("royalblue","green4"),  title = "Met-Caproic acid",xlab = "Days")
Met_C6_all   <-Met_C6_all   + scale_x_continuous(breaks = lipids_N06$ATTRIBUTE_Study_Day)
Met_C6_all   <- Met_C6_all  + theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank()) 
Met_C6_all
SANTAObj[[3]]$general$pval.dist


PC_20_2_all <-santaR_plot(SANTAObj[[6]], colorVect=c("royalblue","green4"),  title = "PC(20:2/0:0)",xlab = "Days")
PC_20_2_all   <-PC_20_2_all   + scale_x_continuous(breaks = lipids_N06$ATTRIBUTE_Study_Day)
PC_20_2_all   <- PC_20_2_all  + theme(panel.grid.major.x = element_blank(),
                                      panel.grid.minor.x = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.minor.y = element_blank()) 
PC_20_2_all
SANTAObj[[6]]$general$pval.dist


Stearidonic_acid <-santaR_plot(SANTAObj[[7]], colorVect=c("royalblue","green4"),  title = "Stearidonic acid ",xlab = "Days")
Stearidonic_acid   <-Stearidonic_acid   + scale_x_continuous(breaks = lipids_N06$ATTRIBUTE_Study_Day)
Stearidonic_acid   <- Stearidonic_acid   + theme(panel.grid.major.x = element_blank(),
                                                 panel.grid.minor.x = element_blank(),
                                                 panel.grid.major.y = element_blank(),
                                                 panel.grid.minor.y = element_blank()) 
Stearidonic_acid 
SANTAObj[[7]]$general$pval.dist


