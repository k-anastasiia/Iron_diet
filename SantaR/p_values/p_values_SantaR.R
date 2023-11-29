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

setwd('~/Desktop/GitHub/SantaR/main_figure_p_values')
getwd()

metadata <- read_csv("md_new.csv")
metadata <- metadata[-1:-10,-1] #exclude bedding and blanks
metadata <- metadata[,-34:-37]
metadata <- metadata[,-2:-30]
hits <- read_delim("hits_all.csv",
                   delim = "\t", escape_double = FALSE,
                   trim_ws = TRUE)
Normalised <- read_csv("Normalised_Quant_table.csv")

hits$Unique_Compound_Name <- hits$Compound_Name
duplicates <- duplicated(hits$Compound_Name)
counter <- 1
for (i in 2:length(hits$Compound_Name)) {
  if (duplicates[i]) {
    counter <- counter + 1
    hits$Unique_Compound_Name[i] <- paste0(hits$Compound_Name[i], "_", counter)
  }
}

colnames(Normalised)[1] ="row_ID"
Normalised[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(Normalised$'row_ID', '_', 3)
Normalised$"row_ID"<-as.numeric(Normalised$"row_ID")
hits$"Row_ID"<-as.numeric(hits$"Row_ID")
Merged_Hits_result <- merge(Normalised, hits, by.x = "row_ID", by.y = "Row_ID",  
                            all.x = FALSE, all.y = FALSE)
Merged_Hits_result_cut <- Merged_Hits_result [,-189:-195]
Merged_Hits_result_cut_f <- Merged_Hits_result_cut  %>% column_to_rownames ('Unique_Compound_Name') %>% t() %>% as.data.frame %>% rownames_to_column
Merged_Hits_result_cut_f <- Merged_Hits_result_cut_f[-1,]

data_merge_santaR <- merge(Merged_Hits_result_cut_f, metadata, by.x = "rowname", by.y = "filename",  
                           all.x = FALSE, all.y = FALSE)
write.csv(data_merge_santaR, file = "library_hits_new_merged.csv", row.names = FALSE)

data_merge_santaR <- data_merge_santaR[order(data_merge_santaR$ATTRIBUTE_Diet), ]
data_merge_santaR_ND <- data_merge_santaR[-1:-62,]
inputData_ND     <- data.frame(data_merge_santaR_ND[,2:368])
ind_ND           <- data_merge_santaR_ND[["ATTIBUTE_mouse_no"]]
time_ND          <- as.numeric(data_merge_santaR_ND[["ATTRIBUTE_Study_Day"]])
group_ND         <- data_merge_santaR_ND[["ATTRIBUTE_Diet"]]
SANTAObj_ND  <- santaR_auto_fit(inputData_ND, ind_ND, time_ND, group_ND, df=5)
SANTAObj_ND[[15]]$general$pval.dist #ser-leu
SANTAObj_ND[[40]]$general$pval.dist #ile-gly-ile
SANTAObj_ND[[50]]$general$pval.dist #leu-ek
SANTAObj_ND[[77]]$general$pval.dist #arg-bmca


data_merge_santaR_N06 <- data_merge_santaR[-63:-124,]
inputData_N06     <- data.frame(data_merge_santaR_N06 [,2:368])
ind_N06           <- data_merge_santaR_N06 [["ATTIBUTE_mouse_no"]]
time_N06          <- as.numeric(data_merge_santaR_N06 [["ATTRIBUTE_Study_Day"]])
group_N06          <- data_merge_santaR_N06 [["ATTRIBUTE_Diet"]]
SANTAObj_N06   <- santaR_auto_fit(inputData_N06 , ind_N06 , time_N06 , group_N06 , df=5)



SANTAObj_N06[[10]]$general$pval.dist #inosine
SANTAObj_N06[[9]]$general$pval.dist #guanosine
SANTAObj_N06[[12]]$general$pval.dist #deoxyguanosine
SANTAObj_N06[[120]]$general$pval.dist #stearidonic
SANTAObj_N06[[326]]$general$pval.dist #pc_20:2 1
SANTAObj_N06[[340]]$general$pval.dist #pc_20:2 2
SANTAObj_N06[[86]]$general$pval.dist #met-c:0 1
SANTAObj_N06[[99]]$general$pval.dist #met-c:0 2
SANTAObj_N06[[214]]$general$pval.dist #glu-cdca
SANTAObj_N06[[197]]$general$pval.dist #phe-gmca

