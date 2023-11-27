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

getwd() #to see the working directory 
setwd("/home/anastasiia/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance")
#scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/CLR_Scaled.csv")
metadata<- read_csv("metadata_full (ak).csv")
AA_BA <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/AA_BA.csv")
AA_FA <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/AA_FA.csv")
Normalised <- read_csv("Normalised_Quant_table.csv")
dummy <- read_delim("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/dummy.csv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)


#BA
AA_BA$Unique_Compound_Name <- AA_BA$Compound_Name
duplicates <- duplicated(AA_BA$Compound_Name)
counter <- 1
for (i in 2:length(AA_BA$Compound_Name)) {
  if (duplicates[i]) {
    counter <- counter + 1
    AA_BA$Unique_Compound_Name[i] <- paste0(AA_BA$Compound_Name[i], "_", counter)
  }
}
view(AA_BA)


#NORMALISED TABLE
# Split name column into several names
colnames(Normalised)[1] ="row_ID"
Normalised[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(Normalised$'row_ID', '_', 3)

#Merge Tables
Merged_Hits_AA_BA <- merge(Normalised, AA_BA, by.x = "row_ID", by.y = "Row_ID",  
                        all.x = FALSE, all.y = FALSE)
Merged_Hits_AA_BA_cut <- Merged_Hits_AA_BA [,-189:-201]
Merged_Hits_AA_BA_cut_f <- Merged_Hits_AA_BA_cut  %>% column_to_rownames ('Unique_Compound_Name') %>% t() %>% as.data.frame %>% rownames_to_column
Merged_Hits_AA_BA_cut_f <- Merged_Hits_AA_BA_cut_f[-1,]
metadata <- metadata[-188,-34:-37]
metadata <- metadata[,-2:-30]
data_merge_santaR <- merge(Merged_Hits_AA_BA_cut_f, metadata, by.x = "rowname", by.y = "...1",  
                           all.x = FALSE, all.y = FALSE)

#all diets
inputData     <- data.frame(data_merge_santaR[,2:38])
#inputData     <- data.frame(data_merge_santaR[,22:36])
ind           <- data_merge_santaR[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

santaR_plot(SANTAObj[[1]])
santaR_plot(SANTAObj[[2]])
santaR_plot(SANTAObj[[3]])
santaR_plot(SANTAObj[[4]])
santaR_plot(SANTAObj[[5]])
santaR_plot(SANTAObj[[6]])
santaR_plot(SANTAObj[[7]])
santaR_plot(SANTAObj[[8]])
santaR_plot(SANTAObj[[9]])
santaR_plot(SANTAObj[[10]])
santaR_plot(SANTAObj[[11]])
santaR_plot(SANTAObj[[12]])
santaR_plot(SANTAObj[[13]])
santaR_plot(SANTAObj[[14]])
santaR_plot(SANTAObj[[15]])
santaR_plot(SANTAObj[[16]])
santaR_plot(SANTAObj[[17]])
santaR_plot(SANTAObj[[18]])
santaR_plot(SANTAObj[[19]])
santaR_plot(SANTAObj[[20]])
santaR_plot(SANTAObj[[21]])
santaR_plot(SANTAObj[[22]])
santaR_plot(SANTAObj[[23]])
santaR_plot(SANTAObj[[24]])
santaR_plot(SANTAObj[[25]])
santaR_plot(SANTAObj[[26]])
santaR_plot(SANTAObj[[27]])
santaR_plot(SANTAObj[[28]])
santaR_plot(SANTAObj[[29]])
santaR_plot(SANTAObj[[30]])
santaR_plot(SANTAObj[[31]])
santaR_plot(SANTAObj[[32]])
santaR_plot(SANTAObj[[33]])
santaR_plot(SANTAObj[[34]])
santaR_plot(SANTAObj[[35]])
santaR_plot(SANTAObj[[36]])
santaR_plot(SANTAObj[[37]])

p1
p2
p3
#grid.arrange(p1, p2, ncol=2) # force them side by side




#FA
AA_FA$Unique_Compound_Name <- AA_FA$Compound_Name
duplicates <- duplicated(AA_FA$Compound_Name)
counter <- 1
for (i in 2:length(AA_FA$Compound_Name)) {
  if (duplicates[i]) {
    counter <- counter + 1
    AA_FA$Unique_Compound_Name[i] <- paste0(AA_FA$Compound_Name[i], "_", counter)
  }
}
view(AA_FA)


#NORMALISED TABLE
# Split name column into several names
colnames(Normalised)[1] ="row_ID"
Normalised[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(Normalised$'row_ID', '_', 3)

#Merge Tables
Merged_Hits_AA_FA <- merge(Normalised, AA_FA, by.x = "row_ID", by.y = "Row_ID",  
                           all.x = FALSE, all.y = FALSE)
Merged_Hits_AA_FA_cut <- Merged_Hits_AA_FA [,-189:-201]
Merged_Hits_AA_FA_cut_f <- Merged_Hits_AA_FA_cut  %>% column_to_rownames ('Unique_Compound_Name') %>% t() %>% as.data.frame %>% rownames_to_column
Merged_Hits_AA_FA_cut_f <- Merged_Hits_AA_FA_cut_f[-1,]
metadata<- read_csv("metadata_full (ak).csv")
metadata <- metadata[-188,-34:-37]
metadata <- metadata[,-2:-30]
data_merge_santaR <- merge(Merged_Hits_AA_FA_cut_f, metadata, by.x = "rowname", by.y = "...1",  
                           all.x = FALSE, all.y = FALSE)

#all diets
inputData     <- data.frame(data_merge_santaR[,2:7])
#inputData     <- data.frame(data_merge_santaR[,22:36])
ind           <- data_merge_santaR[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santaR[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santaR[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

santaR_plot(SANTAObj[[1]])
santaR_plot(SANTAObj[[2]])
santaR_plot(SANTAObj[[3]])
santaR_plot(SANTAObj[[4]])
santaR_plot(SANTAObj[[5]])
santaR_plot(SANTAObj[[6]])
santaR_plot(SANTAObj[[7]])
santaR_plot(SANTAObj[[8]])
santaR_plot(SANTAObj[[9]])
santaR_plot(SANTAObj[[10]])
santaR_plot(SANTAObj[[11]])
santaR_plot(SANTAObj[[12]])
santaR_plot(SANTAObj[[13]])
santaR_plot(SANTAObj[[14]])
santaR_plot(SANTAObj[[15]])
santaR_plot(SANTAObj[[16]])
santaR_plot(SANTAObj[[17]])
santaR_plot(SANTAObj[[18]])
santaR_plot(SANTAObj[[19]])
santaR_plot(SANTAObj[[20]])
santaR_plot(SANTAObj[[21]])
santaR_plot(SANTAObj[[22]])
santaR_plot(SANTAObj[[23]])
santaR_plot(SANTAObj[[24]])
santaR_plot(SANTAObj[[25]])
santaR_plot(SANTAObj[[26]])
santaR_plot(SANTAObj[[27]])
santaR_plot(SANTAObj[[28]])
santaR_plot(SANTAObj[[29]])
santaR_plot(SANTAObj[[30]])
santaR_plot(SANTAObj[[31]])
santaR_plot(SANTAObj[[32]])
santaR_plot(SANTAObj[[33]])
santaR_plot(SANTAObj[[34]])
santaR_plot(SANTAObj[[35]])
santaR_plot(SANTAObj[[36]])
santaR_plot(SANTAObj[[37]])



