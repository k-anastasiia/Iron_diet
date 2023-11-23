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

setwd('~/Desktop/GitHub/PLS_DA')
getwd()


metadata<- read_csv("~/Desktop/GitHub/PLS_DA/md_new.csv")
metadata <- metadata[-1:-10,-1] #exclude bedding and blanks
colnames(metadata)[1] <- "...1"
scaled_table <- read_csv("~/Desktop/GitHub/PLS_DA/CLR_scaled.csv")

metadata$ATTRIBUTE_Study_Day <- as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata <- metadata[order(metadata$ATTRIBUTE_Study_Day),]
metadata_before <- metadata[13:107,] #Exclude day 0

data_merge_before <- metadata_before %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_before %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_before$ATTRIBUTE_Diet # response variables 


plsda_before <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 7, scale = FALSE)
plotIndiv(plsda_before,comp = c(1,2), col = c("royalblue",'green4', "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_before, group = as.character(metadata_before$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA",
          legend.title = "Diet",
          X.label = "Component 1 (23%)", Y.label = "Component 2 (6%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)


#undergo performance evaluation in order to tune the number of components to use
perf.plsda <- perf(plsda_before, validation = "Mfold", 
                   folds = 5, nrepeat = 10, # use repeated cross-validation
                   progressBar = FALSE, auc = TRUE) # include AUC values
plot(perf.plsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")


plotLoadings(plsda_before, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#06 VS DEF
metadata_D06_b <- metadata_before[order(metadata_before$ATTRIBUTE_Diet), ]
metadata_D06_b <- metadata_D06_b[-64:-95,] #exclude 

data_merge_D06_b <- metadata_D06_b %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_D06_b %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_b$ATTRIBUTE_Diet # response variables 

plsda_D06_b <- mixOmics::plsda(X, Y, ncomp = 7,scale = FALSE) 
plotIndiv(plsda_D06_b,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_D06_b, group = as.character(metadata_D06_b$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA - Fe-deficient vs Fe-overload",
          legend.title = "Diet",
          X.label = "Component 1 (27%)", Y.label = "Component 2 (15%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_D06_b, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

#Normal VS DEF
metadata_ND_b <- metadata_before[order(metadata_before$ATTRIBUTE_Diet), ]
metadata_ND_b <- metadata_ND_b[-1:-32,] #exclude 
data_merge_ND_b <- metadata_ND_b %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_ND_b %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_b$ATTRIBUTE_Diet # response variables 

plsda_ND_b <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_ND_b,comp = c(1,2), col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_ND_b, group = as.character(metadata_ND_b$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (9%)", Y.label = "Component 2 (20%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_ND_b, comp = 1, title = 'Loadings on comp 1', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Normal VS 06
metadata_N06_b <- metadata_before[order(metadata_before$ATTRIBUTE_Diet), ]
metadata_N06_b <- metadata_N06_b[-33:-63,] #exclude 

data_merge_N06_b <- metadata_N06_b %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_N06_b %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_b$ATTRIBUTE_Diet # response variables 

plsda_N06_b <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE)  
plotIndiv(plsda_N06_b,comp = c(1,2), col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_N06_b, group = as.character(metadata_N06_b$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (27%)", Y.label = "Component 2 (13%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_N06_b, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue","green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#############AFTER DIET SWITCH######################     
metadata_after <- metadata[108:187,] 
data_merge_after <- metadata_after %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_after %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_after$ATTRIBUTE_Diet # response variables 
plsda_after <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_after,comp = c(1,2), col = c("royalblue",'green4', "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotIndiv(plsda_after, group = as.character(metadata_after$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA",
          legend.title = "Diet",
          X.label = "Component 1 (29%)", Y.label = "Component 2 (10%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)


plotLoadings(plsda_after, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#06 VS DEF
metadata_D06_a <- metadata_after[order(metadata_after$ATTRIBUTE_Diet), ]
metadata_D06_a <- metadata_D06_a[-54:-80,] #exclude 

data_merge_D06_a <- metadata_D06_a %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_D06_a %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_a$ATTRIBUTE_Diet # response variables 
plsda_D06_a <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_D06_a,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_D06_a, group = as.character(metadata_D06_a$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA - Fe-deficient vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (38%)", Y.label = "Component 2 (8%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_D06_a, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Normal VS DEF
metadata_ND_a <- metadata_after[order(metadata_after$ATTRIBUTE_Diet), ]
metadata_ND_a <- metadata_ND_a[-1:-26,] 
data_merge_ND_a <- metadata_ND_a %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_ND_a %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_a$ATTRIBUTE_Diet # response variables 
plsda_ND_a <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_ND_a,comp = c(1,2), col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_ND_a, group = as.character(metadata_ND_a$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (20%)", Y.label = "Component 2 (16%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

plotLoadings(plsda_ND_a, comp = 1, title = 'Loadings on comp 1', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

#Normal VS 06
metadata_N06_a <- metadata_after[order(metadata_after$ATTRIBUTE_Diet), ]
metadata_N06_a <- metadata_N06_a[-27:-53,] #exclude 
data_merge_N06_a <- metadata_N06_a %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_N06_a %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_a$ATTRIBUTE_Diet # response variables 
plsda_N06_a <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE) 
plotIndiv(plsda_N06_a,comp = c(1,2), col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")

plotIndiv(plsda_N06_a, group = as.character(metadata_N06_a$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (36%)", Y.label = "Component 2 (8%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)


plotLoadings(plsda_N06_a, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

