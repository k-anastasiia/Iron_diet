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

# Load example data
X <- data_merge %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge$ATTRIBUTE_Diet # response variables 
plsda <- mixOmics::plsda(X, Y, ncomp = 8, scale = FALSE) 
plotIndiv(plsda, col = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotVar(plsda)
l <- plotLoadings(plsda, contrib = "max")


# http://mixomics.org/case-studies/splsda-srbct-case-study/
#undergo performance evaluation in order to tune the number of components to use
perf.plsda <- perf(plsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.plsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

# what is the optimal value of components according to perf()
perf.plsda$choice.ncomp 


library(caret)

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.plsda <- tune.plsda(X, Y, ncomp = 8, # calculate for first 8 components
                         validation = 'Mfold',
                         folds = 5, nrepeat = 10, # use repeated cross-validation
                         dist = 'max.dist', # use max.dist measure
                         measure = "BER", # use balanced error rate of dist measure
                         test.keepX = list.keepX,
                         cpus = 2) # allow for paralleliation to decrease runtime

plot(tune.plsda, col = color.jet(4)) # plot output of variable number tuning






#select the variables that have the highest correlation with a given component of the PLS-DA model
dr1<-as.list(selectVar(plsda, comp=1)) 
View(dr1)
dr2<-as.list(selectVar(plsda, comp=2)) 
View(dr2)

# select variables based on their VIP scores 
vip(plsda) 
View(vip(plsda))


















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
plsda_before <- mixOmics::plsda(X, Y, ncomp = 10, scale = FALSE)
plotIndiv(plsda_before, col = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotVar(plsda_before)
selectVar(plsda_before)


# undergo performance evaluation in order to tune the number of components to use
perf.plsda_before <- perf(plsda_before, validation = "Mfold", 
                   folds = 5, nrepeat = 10, # use repeated cross-validation
                   progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.plsda_before, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

# what is the optimal value of components according to perf()
perf.plsda_before$choice.ncomp 

#select the variables that have the highest correlation with a given component of the PLS-DA model
dr1b<-as.list(selectVar(plsda_before, comp=1)) 
View(dr1b)
dr2b<-as.list(selectVar(plsda_before, comp=2)) 
View(dr2b)

# select variables based on their VIP scores 
vip(plsda_before) 
View(vip(plsda_before))




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
plsda_after <- mixOmics::plsda(X, Y, ncomp = 10, scale = FALSE)
plotIndiv(plsda_after, col = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
plotVar(plsda_after)
selectVar(plsda_after)


# undergo performance evaluation in order to tune the number of components to use
perf.plsda_after <- perf(plsda_after, validation = "Mfold", 
                   folds = 5, nrepeat = 10, # use repeated cross-validation
                   progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.plsda_after, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

# what is the optimal value of components according to perf()
perf.plsda_after$choice.ncomp 








#select the variables that have the highest correlation with a given component of the PLS-DA model
dr1a<-as.list(selectVar(plsda_after, comp=1)) 
View(dr1a)
dr2a<-as.list(selectVar(plsda_after, comp=2)) 
View(dr2a)

# select variables based on their VIP scores 
vip(plsda_after) 
View(vip(plsda_after))




#################10-fold cross-validation
#http://mixomics.org/wp-content/uploads/2014/08/Running_perf_function4.pdf

tune.plsda = perf(plsda, dist = "max.dist", validation = "Mfold", folds = 10,
                  progressBar = FALSE)
tune.plsda$error.rate
view(summary(tune.plsda))

