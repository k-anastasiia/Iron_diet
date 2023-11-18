#https://github.com/DaniRuizPerez/SoYouThinkYouCanPLS-DA_Public/blob/master/AllMethods/PLSDASignalAndNoise.r
#https://github.com/simonezuffa/Casein_Manuscript/blob/main/Casein_Analysis.Rmd
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
library(ggrepel)

setwd("/home/anastasiia/Desktop/mouse_studies/PLS-DA/Separate_Days")
getwd() #to see the working directory 

metadata_full_ak<- read_csv("~/Desktop/mouse_studies/PLS-DA/Separate_Days/metadata_full (ak).csv")
metadata <- metadata_full_ak[-188,] #exclude bedding
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/Separate_Days/CLR_Scaled.csv")
metadata$ATTRIBUTE_Study_Day <- as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata <- metadata[order(metadata$ATTRIBUTE_Study_Day),]
metadata_0 <- metadata[1:12,] 
metadata_47 <- metadata[95:107,] 
metadata_75 <- metadata[147:161,] 
metadata_96 <- metadata[173:187,] 
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]



#############################Day_0
data_merge_0 <- metadata_0 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_0 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_0$ATTRIBUTE_Diet # response variables 

plsda_0 <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 4, scale = FALSE)
plotIndiv(plsda_0, col = c("royalblue", "red3","green4"), title = 'Day 0',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_0, group = as.character(metadata_0$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA Day 0",
          legend.title = "Diet",
          X.label = "Component 1 (13%)", Y.label = "Component 2 (11%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_0, comp = 1, title = 'Day 0', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

#06 VS DEF
metadata_D06_0 <- metadata_0[order(metadata_0$ATTRIBUTE_Diet), ]
metadata_D06_0 <- metadata_D06_0[-9:-12,] 

data_merge_D06_0 <- metadata_D06_0 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_D06_0 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_0$ATTRIBUTE_Diet # response variables 


plsda_D06_0<- mixOmics::plsda(X, Y, ncomp = 4,scale = FALSE) 
plotIndiv(plsda_D06_0,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_D06_0, group = as.character(metadata_D06_0$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA Day 0 - Fe-deficient vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (19%)", Y.label = "Component 2 (19%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)
plotLoadings(plsda_D06_0, comp = 1, title = 'Day 0', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 25, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.15, 0.15), ylim = c(-2, 2))


#Normal VS DEF
metadata_ND_0 <- metadata_0[order(metadata_0$ATTRIBUTE_Diet), ]
metadata_ND_0 <- metadata_ND_0[-1:-4,] #exclude 
data_merge_ND_0 <- metadata_ND_0 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_ND_0 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_0$ATTRIBUTE_Diet # response variables 

plsda_ND_0 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE) 
plotIndiv(plsda_ND_0, col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_ND_0, group = as.character(metadata_ND_0$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA Day 0 - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (20%)", Y.label = "Component 2 (11%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_ND_0, comp = 1, title = 'Day 0', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.15, 0.15), ylim = c(-2, 2))


#Normal VS 06
metadata_N06_0 <- metadata_0[order(metadata_0$ATTRIBUTE_Diet), ]
metadata_N06_0 <- metadata_N06_0[-5:-8,] 
data_merge_N06_0 <- metadata_N06_0 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_N06_0 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_0$ATTRIBUTE_Diet # response variables

plsda_N06_0 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE)  
plotIndiv(plsda_N06_0, col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_N06_0, group = as.character(metadata_N06_0$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA Day 0 - Normal vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (14%)", Y.label = "Component 2 (28%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_N06_0, comp = 1, title = 'Day 0', legend.color = c("royalblue","green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.15, 0.15), ylim = c(-2, 2))





#############################Day_47
data_merge_47 <- metadata_47 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_47 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_47$ATTRIBUTE_Diet # response variables 

plsda_47 <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 4, scale = FALSE)
plotIndiv(plsda_47, col = c("royalblue","red3", "green4"), title = '47',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_47, group = as.character(metadata_47$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA Day 47",
          legend.title = "Diet",
          X.label = "Component 1 (42%)", Y.label = "Component 2 (15%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)
vip_47<-vip(plsda_47)
view(vip_47)
write.csv(vip_47, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/vip_47.csv",row.names = TRUE)

plotLoadings(plsda_47, comp = 1, title = 'Day 47', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

########06 VS DEF
metadata_D06_47 <- metadata_47[order(metadata_0$ATTRIBUTE_Diet), ]
metadata_D06_47 <- metadata_D06_47[-1:-5,] 

data_merge_D06_47 <- metadata_D06_47 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_D06_47 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_47$ATTRIBUTE_Diet # response variables 


plsda_D06_47<- mixOmics::plsda(X, Y, ncomp = 4,scale = FALSE) 
plotIndiv(plsda_D06_47,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_D06_47, group = as.character(metadata_D06_47$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA Day 47 - Fe-deficient vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (63%)", Y.label = "Component 2 (10%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)
plotLoadings(plsda_D06_47, comp = 1, title = 'Day 47', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 25, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

######Normal VS DEF
metadata_ND_47 <- metadata_47[order(metadata_47$ATTRIBUTE_Diet), ]
metadata_ND_47 <- metadata_ND_47[-1:-4,] #exclude 
data_merge_ND_47 <- metadata_ND_47 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_ND_47 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_47$ATTRIBUTE_Diet # response variables 

plsda_ND_47 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE) 
plotIndiv(plsda_ND_47, col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_ND_47, group = as.character(metadata_ND_47$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA Day 47 - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (33%)", Y.label = "Component 2 (18%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_ND_47, comp = 1, title = 'Day 47', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))
#Merge
day47_ND <- read_csv("day47_ND.csv")
day47_ND <- day47_ND %>% t() %>% as.data.frame %>% rownames_to_column("Row_ID")

Merged_Hits_day47_ND<- merge(day47_ND , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                              all.x = FALSE, all.y = FALSE)
View(Merged_Hits_day47_ND)
write.csv(Merged_Hits_day47_ND , "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Hits_day47_ND.csv",row.names = TRUE)

Merged_canopus_day47_ND  <- merge(day47_ND, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                   all.x = FALSE, all.y = FALSE)
View(Merged_canopus_day47_ND)
write.csv(Merged_canopus_day47_ND, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Merged_canopus_day47_ND.csv",row.names = TRUE)

######Normal VS 06
metadata_N06_47 <- metadata_47[order(metadata_47$ATTRIBUTE_Diet), ]
metadata_N06_47 <- metadata_N06_47[-5:-8,] 
data_merge_N06_47 <- metadata_N06_47 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_N06_47 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_47$ATTRIBUTE_Diet # response variables

plsda_N06_47 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE)  
plotIndiv(plsda_N06_47, col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_N06_47, group = as.character(metadata_N06_47$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA Day 47 - Normal vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (52%)", Y.label = "Component 2 (15%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_N06_47, comp = 1, title = 'Day 47', legend.color = c("royalblue","green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))
#Merge
day47_N06 <- read_csv("day47_N06.csv")
day47_N06 <- day47_N06 %>% t() %>% as.data.frame %>% rownames_to_column("Row_ID")

Merged_Hits_day47_N06<- merge(day47_N06 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                              all.x = FALSE, all.y = FALSE)
View(Merged_Hits_day47_N06)
write.csv(Merged_Hits_day47_N06 , "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Hits_day47_N06.csv",row.names = TRUE)

Merged_canopus_day47_N06  <- merge(day47_N06, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                   all.x = FALSE, all.y = FALSE)
View(Merged_canopus_day47_N06)
write.csv(Merged_canopus_day47_N06, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Merged_canopus_day47_N06.csv",row.names = TRUE)


#############################Day_75
data_merge_75 <- metadata_75 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_75 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_75$ATTRIBUTE_Diet # response variables 

plsda_75 <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 4, scale = FALSE)
plotIndiv(plsda_75, col = c("royalblue","red3", "green4"), title = '75',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_75, group = as.character(metadata_75$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA Day 75",
          legend.title = "Diet",
          X.label = "Component 1 (47%)", Y.label = "Component 2 (14%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

vip_75<-vip(plsda_75)
view(vip_75)
write.csv(vip_75, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/vip_75.csv",row.names = TRUE)

plotLoadings(plsda_75, comp = 1, title = 'Day 75', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.2, 0.2), ylim = c(-2, 2))

#####06 VS DEF
metadata_D06_75 <- metadata_75[order(metadata_75$ATTRIBUTE_Diet), ]
metadata_D06_75 <- metadata_D06_75[-11:-15,] 

data_merge_D06_75 <- metadata_D06_75 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_D06_75 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_75$ATTRIBUTE_Diet # response variables 


plsda_D06_75<- mixOmics::plsda(X, Y, ncomp = 4,scale = FALSE) 
plotIndiv(plsda_D06_75,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_D06_75, group = as.character(metadata_D06_75$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA Day 75 - Fe-deficient vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (58%)", Y.label = "Component 2 (8%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)
plotLoadings(plsda_D06_75, comp = 1, title = 'Day 75', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 25, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#####Normal VS DEF
metadata_ND_75 <- metadata_75[order(metadata_75$ATTRIBUTE_Diet), ]
metadata_ND_75 <- metadata_ND_75[-1:-5,] #exclude 
data_merge_ND_75 <- metadata_ND_75 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_ND_75 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_75$ATTRIBUTE_Diet # response variables 

plsda_ND_75 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE) 
plotIndiv(plsda_ND_75, col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_ND_75, group = as.character(metadata_ND_75$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA Day 47 - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (40%)", Y.label = "Component 2 (14%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_ND_75, comp = 1, title = 'Day 75', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

#Merge Tables
day75_ND <- read_csv("day75_ND.csv")
day75_ND <- day75_ND %>% t() %>% as.data.frame %>% rownames_to_column("Row_ID")

Merged_Hits_day75_ND<- merge(day75_ND , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                              all.x = FALSE, all.y = FALSE)
View(Merged_Hits_day75_ND)
write.csv(Merged_Hits_day75_ND , "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Hits_day75_ND.csv",row.names = TRUE)
Merged_canopus_day75_ND  <- merge(day75_ND, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                   all.x = FALSE, all.y = FALSE)
View(Merged_canopus_day75_ND)
write.csv(Merged_canopus_day75_ND, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Merged_canopus_day75_ND.csv",row.names = TRUE)


#####Normal VS 06
metadata_N06_75 <- metadata_75[order(metadata_75$ATTRIBUTE_Diet), ]
metadata_N06_75 <- metadata_N06_75[-6:-10,] 
data_merge_N06_75 <- metadata_N06_75 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_N06_75 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_75$ATTRIBUTE_Diet # response variables

plsda_N06_75 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE)  
plotIndiv(plsda_N06_75, col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_N06_75, group = as.character(metadata_N06_75$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA Day 75 - Normal vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (57%)", Y.label = "Component 2 (12%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_N06_75, comp = 1, title = 'Day 75', legend.color = c("royalblue","green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))
#Merge
day75_N06 <- read_csv("day75_N06.csv")
day75_N06 <- day75_N06 %>% t() %>% as.data.frame %>% rownames_to_column("Row_ID")

Merged_Hits_day75_N06<- merge(day75_N06 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                             all.x = FALSE, all.y = FALSE)
View(Merged_Hits_day75_N06)
write.csv(Merged_Hits_day75_N06 , "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Hits_day75_N06.csv",row.names = TRUE)

Merged_canopus_day75_N06  <- merge(day75_N06, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_canopus_day75_N06)
write.csv(Merged_canopus_day75_N06, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/Merged_canopus_day75_N06.csv",row.names = TRUE)

#############################Day_96
data_merge_96 <- metadata_96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_96 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_96$ATTRIBUTE_Diet # response variables 

plsda_96 <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 4, scale = FALSE)
plotIndiv(plsda_96, col = c("royalblue","red3", "green4"), title = '96',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_96, group = as.character(metadata_96$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA Day 96",
          legend.title = "Diet",
          X.label = "Component 1 (16%)", Y.label = "Component 2 (11%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

vip_96<-vip(plsda_96)
view(vip_96)
write.csv(vip_96, "~/Desktop/mouse_studies/PLS-DA/Separate_Days/vip_96.csv",row.names = TRUE)


plotLoadings(plsda_96, comp = 1, title = 'Day 96', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.2, 0.2), ylim = c(-2, 2))

#06 VS DEF
metadata_D06_96 <- metadata_96[order(metadata_75$ATTRIBUTE_Diet), ]
metadata_D06_96 <- metadata_D06_96[-11:-15,] 

data_merge_D06_96 <- metadata_D06_96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_D06_96 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_96$ATTRIBUTE_Diet # response variables 


plsda_D06_96<- mixOmics::plsda(X, Y, ncomp = 4,scale = FALSE) 
plotIndiv(plsda_D06_96,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_D06_96, group = as.character(metadata_D06_96$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA Day 96 - Fe-deficient vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (21%)", Y.label = "Component 2 (16%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)
plotLoadings(plsda_D06_96, comp = 1, title = 'Day 96', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 25, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Normal VS DEF
metadata_ND_96 <- metadata_96[order(metadata_96$ATTRIBUTE_Diet), ]
metadata_ND_96 <- metadata_ND_96[-1:-5,] #exclude 
data_merge_ND_96 <- metadata_ND_96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_ND_96%>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_96$ATTRIBUTE_Diet # response variables 

plsda_ND_96 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE) 
plotIndiv(plsda_ND_96, col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_ND_96, group = as.character(metadata_ND_96$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA Day 96 - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (26%)", Y.label = "Component 2 (16%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_ND_96, comp = 1, title = 'Day 96', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Normal VS 06
metadata_N06_96 <- metadata_96[order(metadata_96$ATTRIBUTE_Diet), ]
metadata_N06_96 <- metadata_N06_96[-6:-10,] 
data_merge_N06_96 <- metadata_N06_96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_N06_96 %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_965$ATTRIBUTE_Diet # response variables

plsda_N06_96 <- mixOmics::plsda(X, Y, ncomp = 4, scale = FALSE)  
plotIndiv(plsda_N06_96, col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_N06_96, group = as.character(metadata_N06_96$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA Day 96 - Normal vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (18%)", Y.label = "Component 2 (21%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_N06_96, comp = 1, title = 'Day 96', legend.color = c("royalblue","green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


