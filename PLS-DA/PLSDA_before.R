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

setwd("/home/anastasiia/Desktop/mouse_studies/PLS-DA")
getwd() #to see the working directory 

metadata_full_ak<- read_csv("~/Desktop/mouse_studies/PLS-DA/metadata_full (ak).csv")
metadata <- metadata_full_ak[-188,] #exclude bedding
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


#############BEFORE DIET SWITCH######################      

metadata$ATTRIBUTE_Study_Day <- as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata <- metadata[order(metadata$ATTRIBUTE_Study_Day),]
metadata_before <- metadata[13:107,] 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_before <- metadata_before %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

X <- data_merge_before %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_before$ATTRIBUTE_Diet # response variables 


###PCA
design <- data.frame(samp = metadata_before$ATTIBUTE_mouse_no)
print(design)

result.pca.multi <- pca(X, ncomp = 2,  scale = FALSE)   # run the method
plotIndiv(result.pca.multi, group = as.character(metadata_before$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PCA",
          legend.title = "Diet",
          X.label = "Component 1 (23%)", Y.label = "Component 2 (6%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)  # plot the samples






#plsda_before <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 2, multilevel=metadata_before$ATTRIBUTE_Study_Day, scale = FALSE)
plsda_before <- mixOmics::plsda(X, Y, max.iter = 500, ncomp = 2, scale = FALSE)
plotIndiv(plsda_before,comp = c(1,2), col = c("royalblue",'green4', "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
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

# plot the outcome of performance evaluation across all ten components
plot(perf.plsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")


plsda_before <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE, multilevel = design) 
plotIndiv(plsda_before,comp = c(1,2), col = c("royalblue",'green4', "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_before, group = as.character(metadata_before$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue","red3",'green4'), 
          pch = 20, title = "PLS-DA",
          legend.title = "Diet",
          X.label = "Component 1 (23%)", Y.label = "Component 2 (6%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_before)
plotLoadings(plsda_before, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3",'green4'),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

#Merge Tables
before_all_comp1 <- read_csv("before_all_comp1.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_before_all <- merge(before_all_comp1 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits_before_all)
#write.csv(Merged_Hits_before_all , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_before_all.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_before_all  <- merge(before_all_comp1, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                    all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_before_all)
#write.csv(Merged_Hits2_before_all, "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_before_all.csv",row.names = TRUE)


#06 VS DEF
metadata_D06_b <- metadata_before[order(metadata_before$ATTRIBUTE_Diet), ]
metadata_D06_b <- metadata_D06_b[-64:-95,] #exclude 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_D06_b <- metadata_D06_b %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_D06_b %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_D06_b$ATTRIBUTE_Diet # response variables 

design <- data.frame(samp = metadata_D06_b$ATTIBUTE_mouse_no)

plsda_D06_b <- mixOmics::plsda(X, Y, ncomp = 7,scale = FALSE, multilevel = design) 
plotIndiv(plsda_D06_b,comp = c(1,2), col = c("royalblue", "red3"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")


# Plot PLS-DA
plotIndiv(plsda_D06_b, group = as.character(metadata_D06_b$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c("royalblue", "red3"), 
          pch = 20, title = "PLS-DA - Fe-deficient vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (27%)", Y.label = "Component 2 (16%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)


#plotVar(plsda_D06_b)
plotLoadings(plsda_D06_b, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))

#Merge Tables
before_D06 <- read_csv("before_D06.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_before_D06 <- merge(before_D06 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits_before_D06)
#write.csv(Merged_Hits_before_D06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_before_D06.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_before_D06  <- merge(before_D06, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                    all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_before_D06 )
#write.csv(Merged_Hits2_before_D06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_before_D06.csv",row.names = TRUE)






#Normal VS DEF

metadata_ND_b <- metadata_before[order(metadata_before$ATTRIBUTE_Diet), ]
metadata_ND_b <- metadata_ND_b[-1:-32,] #exclude 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_ND_b <- metadata_ND_b %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_ND_b %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_ND_b$ATTRIBUTE_Diet # response variables 
design <- data.frame(samp = metadata_ND_b$ATTIBUTE_mouse_no)
print(design)


plsda_ND_b <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE,multilevel = design) 

plotIndiv(plsda_ND_b,comp = c(1,2), col = c("red3", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_ND_b, group = as.character(metadata_ND_b$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "red3", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-deficient",
          legend.title = "Diet",
          X.label = "Component 1 (9%)", Y.label = "Component 2 (21%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_ND_b)
plotLoadings(plsda_ND_b, comp = 1, title = 'Loadings on comp 1', legend.color = c("red3", "green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))


#Merge Tables
before_ND <- read_csv("before_ND.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_before_ND <- merge(before_ND , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                all.x = FALSE, all.y = FALSE)
View(Merged_Hits_before_ND)
#write.csv(Merged_Hits_before_ND , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_before_ND.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_before_ND  <- merge(before_ND, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_before_ND )
#write.csv(Merged_Hits2_before_ND , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_before_ND.csv",row.names = TRUE)



#Normal VS 06

metadata_N06_b <- metadata_before[order(metadata_before$ATTRIBUTE_Diet), ]
metadata_N06_b <- metadata_N06_b[-33:-63,] #exclude 
scaled_table <- read_csv("~/Desktop/mouse_studies/PLS-DA/CLR_Scaled.csv")
data_merge_N06_b <- metadata_N06_b %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


X <- data_merge_N06_b %>% dplyr::select(-ATTRIBUTE_Diet)# predictor variable
Y <- data_merge_N06_b$ATTRIBUTE_Diet # response variables 

design <- data.frame(samp = metadata_N06_b $ATTIBUTE_mouse_no)
print(design)


plsda_N06_b <- mixOmics::plsda(X, Y, ncomp = 7, scale = FALSE,multilevel = design)  
plotIndiv(plsda_N06_b,comp = c(1,2), col = c("royalblue", "green4"), title = 'PLSDA results',ellipse = TRUE, legend = TRUE, legend.title = "Diet")
# Plot PLS-DA
plotIndiv(plsda_N06_b, group = as.character(metadata_N06_b$ATTRIBUTE_Diet), 
          legend = TRUE, centroid = TRUE, style = "lattice", ind.names = FALSE, 
          ellipse = TRUE, col.per.group = c( "royalblue", "green4"), 
          pch = 20, title = "PLS-DA - Normal vs Fe-supplemented",
          legend.title = "Diet",
          X.label = "Component 1 (27%)", Y.label = "Component 2 (13%)",
          size.title = 0.5, size.xlabel = 1, size.ylabel = 1, size.axis = 1,
          size.legend.title = 1, size.legend = 1, alpha = 0.8)

#plotVar(plsda_N06_b)
plotLoadings(plsda_N06_b, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue","green4"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))



#Merge Tables
before_N06 <- read_csv("before_N06.csv")
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_before_N06 <- merge(before_N06 , Library_Hits_Refined_AK, by.x = "Row_ID", by.y = "Row_ID",
                                all.x = FALSE, all.y = FALSE)
View(Merged_Hits_before_N06)
#write.csv(Merged_Hits_before_N06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits_before_N06.csv",row.names = TRUE)


canopus <- read_delim("canopus_formula_summary_adducts.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
canopus[c('c_ID', 'date', 'b', 'c', 'd', 'a','row_ID')] <- str_split_fixed(canopus$'id', '_', 7)
canopus1 <- canopus [,c(11, 19,17,15,28)]

Merged_Hits2_before_N06  <- merge(before_N06, canopus1, by.x = "Row_ID", by.y = "row_ID",
                                  all.x = FALSE, all.y = FALSE)
View(Merged_Hits2_before_N06 )
#write.csv(Merged_Hits2_before_N06 , "~/Desktop/mouse_studies/PLS-DA/hits/Hits2_before_N06.csv",row.names = TRUE)

