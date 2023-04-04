#https://github.com/DaniRuizPerez/SoYouThinkYouCanPLS-DA_Public/blob/master/AllMethods/PLSDASignalAndNoise.r
library(mixOmics) 
library(lme4) 
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(caret)
library(tidyverse)
library(Biobase)

setwd("/home/anastasiia/Desktop/mouse_studies/PLS-DA")
getwd() #to see the working directory 

metadata_full_ak<- read_csv("~/Desktop/mouse_studies/PLS-DA/metadata_full (ak).csv")
metadata <- metadata_full_ak[-188,] #exclude bedding
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

# Load example data
X <- data_merge %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge$ATTRIBUTE_Diet # response variables 
plsda <- mixOmics::splsda(X, Y, scale = FALSE)
plotIndiv(plsda, col = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
plotVar(plsda)
        

#############BEFORE DIET SWITCH######################      
metadata_full_ak<- read_csv("~/Desktop/mouse_studies/PLS-DA/metadata_full (ak).csv")
metadata_before <- metadata_full_ak[1:67,] 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_before <- metadata_before %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

# Load example data
X <- data_merge_before %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_before$ATTRIBUTE_Diet # response variables 
plsda_before <- mixOmics::splsda(X, Y, scale = FALSE)
plotIndiv(plsda_before, col = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
plotVar(plsda_before)

#############AFTER DIET SWITCH######################      
metadata_full_ak<- read_csv("~/Desktop/mouse_studies/PLS-DA/metadata_full (ak).csv")
metadata_after <- metadata_full_ak[68:187,] 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_after <- metadata_after %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

# Load example data
X <- data_merge_after %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_after$ATTRIBUTE_Diet # response variables 
plsda_after <- mixOmics::splsda(X, Y, scale = FALSE)
plotIndiv(plsda_after, col = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
plotVar(plsda_after)
