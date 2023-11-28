library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(pheatmap)
library(heatmap3)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(tidyverse)

setwd("/Users/Anastasiia/Desktop/Github/Heatmaps/Peptides")
getwd() 

scaled_table <- read_csv("CLR_Scaled.csv")
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")

# Split name column into several names
colnames(scaled_table_t)[1] ="row_ID"
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'row_ID', '_', 3)

#remove day 0
columns_to_exclude <- c(2,9,18,27,36,87,120,129,139,160,170,178) 
scaled_table_t <- scaled_table_t[, -columns_to_exclude]

Library_Hits_Refined_AK <- read_delim("Peptides_no_blank.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)

Merged_Hits_peptides <- merge(scaled_table_t, Library_Hits_Refined_AK, by.x = "row_ID", by.y = "Row_ID",  
                              all.x = FALSE, all.y = FALSE)
Merged_Hits_peptides_cut <- Merged_Hits_peptides [,-182:-189]
Merged_Hits_peptides_cut <- Merged_Hits_peptides_cut [,-177:-180]
data <- Merged_Hits_peptides_cut %>% column_to_rownames("Compound_Name") 
data2 <- data[, -1]
gLogCpmData = as.matrix(data2)

gAnnotationData = read.table("md_full_noday0.txt", sep = '\t', header = TRUE,
                             fill = TRUE)
gAnnotationData

# Make helper function to map metadata category to color
mapDietToColor<-function(annotations){
  colorsVector = ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="0.6% Fe", 
                        "blue" , ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="Normal", 
                                        "green", "red"))
  return(colorsVector)
}

testHeatmap2 <- function(logCPM, annotations) {    
  sampleColors <- mapDietToColor(annotations)
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)



#CANOPUS
Peptides <- read_csv("canopus_peptides.csv")
colnames(Peptides)[1] ="row_ID"
Peptides[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(Peptides$'row_ID', '_', 7)
columns_to_keep <- c(5,14) 
Peptides <- Peptides[, columns_to_keep]
Peptides <- Peptides[order(Peptides[,1]),]
Peptides <- Peptides[-1:-9,]
Peptides <- Peptides[-293:-317,]
Peptides <- Peptides[-242,]
Peptides <- Peptides[-241:-242,]
Merged_Hits_peptides <- merge(scaled_table_t, Peptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-179]
data2 <- data[, -1]
data2 <- data2[, -176]
gLogCpmData = as.matrix(data2)

gAnnotationData = read.table("md_full_noday0.txt", sep = '\t', header = TRUE,
                             fill = TRUE)
gAnnotationData

mapDietToColor<-function(annotations){
  colorsVector = ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="0.6% Fe", 
                        "blue", ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="Normal", 
                                       "green", "red"))
  return(colorsVector)
}

testHeatmap2 <- function(logCPM, annotations) {    
  sampleColors <- mapDietToColor(annotations)
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

peptides<- testHeatmap2(gLogCpmData, gAnnotationData) 

#Dipeptides
Dipeptides <- Peptides[61:239,]
Merged_Hits_peptides <- merge(scaled_table_t, Dipeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-179]
data2 <- data[, -1]
data2 <- data2[, -176]
gLogCpmData = as.matrix(data2)

mapDietToColor<-function(annotations){
  colorsVector = ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="0.6% Fe", 
                        "blue", ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="Normal", 
                                       "green", "red"))
  return(colorsVector)
}

testHeatmap2 <- function(logCPM, annotations) {    
  sampleColors <- mapDietToColor(annotations)
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

dipeptides<- testHeatmap2(gLogCpmData, gAnnotationData) 

#Tripeptides
Tripeptides <- Peptides[290:706,]
Merged_Hits_peptides <- merge(scaled_table_t, Tripeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-179]
data2 <- data[, -1]
data2 <- data2[, -176]
gLogCpmData = as.matrix(data2)
mapDietToColor<-function(annotations){
  colorsVector = ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="0.6% Fe", 
                        "blue", ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="Normal", 
                                       "green", "red"))
  return(colorsVector)
}
testHeatmap2 <- function(logCPM, annotations) {    
  sampleColors <- mapDietToColor(annotations)
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
tripeptides<- testHeatmap2(gLogCpmData, gAnnotationData) 
