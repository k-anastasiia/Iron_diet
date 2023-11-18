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
library(ggrepel)

setwd("/home/anastasiia/Desktop/mouse_studies/PLS-DA")
getwd() #to see the working directory 

metadata_full_ak<- read_csv("~/Desktop/mouse_studies/PLS-DA/metadata_full (ak).csv")
metadata <- metadata_full_ak[-188,] #exclude bedding
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")




#############AFTER DIET SWITCH######################      

metadata$ATTRIBUTE_Study_Day <- as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata <- metadata[order(metadata$ATTRIBUTE_Study_Day),]
metadata_after <- metadata[108:187,] 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_after <- metadata_after %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_after %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_after$ATTRIBUTE_Diet # response variables 
plsda_after <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_after,comp = c(1,2), col = c("royalblue",'green4', "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_after, group = as.character(metadata_after$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA",
          legend.title = "Diet",
          X.label = "Component 1 (27%)", Y.label = "Component 2 (9%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_after)
plotLoadings(plsda_after, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))



#Merge Tables
after_all_comp1 <- read_csv("after_all_comp1.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_after_all<- merge(after_all_comp1 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                all.x = FALSE, all.y = FALSE)
View(Merged_Hits_after_all)
#write.csv(Merged_Hits_after_all , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_after_all.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_after_all  <- merge(after_all_comp1, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_after_all)
#write.csv(Merged_Hits2_after_all , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_after_all.csv",row.names = TRUE)


#06 VS DEF

metadata_D06_a <- metadata_after[order(metadata_after$ATTRIBUTE_Diet), ]
metadata_D06_a <- metadata_D06_a[-54:-80,] #exclude 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_D06_a <- metadata_D06_a %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_D06_a %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_a$ATTRIBUTE_Diet # response variables 
plsda_D06_a <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_D06_a,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_D06_a, group = as.character(metadata_D06_a$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA - Fe-deficient vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (37%)", Y.label = "Component 2 (8%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_D06_a)
plotLoadings(plsda_D06_a, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Merge Tables
after_D06 <- read_csv("after_D06.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_after_D06 <- merge(after_D06 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                all.x = FALSE, all.y = FALSE)
View(Merged_Hits_after_D06)
#write.csv(Merged_Hits_after_D06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_after_D06.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_after_D06  <- merge(after_D06, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_after_D06 )
#write.csv(Merged_Hits2_after_D06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_after_D06.csv",row.names = TRUE)

#Normal VS DEF

metadata_ND_a <- metadata_after[order(metadata_after$ATTRIBUTE_Diet), ]
metadata_ND_a <- metadata_ND_a[-1:-26,] #exclud
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_ND_a <- metadata_ND_a %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_ND_a %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_a$ATTRIBUTE_Diet # response variables 
plsda_ND_a <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_ND_a,comp = c(1,2), col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_ND_a, group = as.character(metadata_ND_a$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (20%)", Y.label = "Component 2 (16%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_ND_a)
plotLoadings(plsda_ND_a, comp = 1, title = 'Loadings on comp 1', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Merge Tables
after_ND <- read_csv("after_ND.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_after_ND <- merge(after_ND  , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                all.x = FALSE, all.y = FALSE)
View(Merged_Hits_after_ND )
#write.csv(Merged_Hits_after_ND  , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_after_ND.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_after_ND   <- merge(after_ND , canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_after_ND)
#write.csv(Merged_Hits2_after_ND , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_after_ND.csv",row.names = TRUE)



#Normal VS 06
metadata_N06_a <- metadata_after[order(metadata_after$ATTRIBUTE_Diet), ]
metadata_N06_a <- metadata_N06_a[-27:-53,] #exclude 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_N06_a <- metadata_N06_a %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_N06_a %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_a$ATTRIBUTE_Diet # response variables 
plsda_N06_a <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_N06_a,comp = c(1,2), col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_N06_a, group = as.character(metadata_N06_a$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (35%)", Y.label = "Component 2 (9%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_N06_a)
plotLoadings(plsda_N06_a, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))



#Merge Tables
after_N06 <- read_csv("after_N06.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_after_N06 <- merge(after_N06 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                all.x = FALSE, all.y = FALSE)
View(Merged_Hits_after_N06)
#write.csv(Merged_Hits_after_N06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_after_N06.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_after_N06  <- merge(after_N06, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_after_N06 )
#write.csv(Merged_Hits2_after_N06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_after_N06.csv",row.names = TRUE)

