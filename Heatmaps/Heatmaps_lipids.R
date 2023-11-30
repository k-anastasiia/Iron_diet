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


setwd('~/Desktop/GitHub/Heatmaps/Lipids')
getwd()

scaled_table <- read_csv("CLR_Scaled.csv")
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'SampleID', '_', 3)

#remove day 0
columns_to_exclude <- c(2,9,18,27,36,87,120,129,139,160,170,178) 
scaled_table_t <- scaled_table_t[, -columns_to_exclude]

lipids <- read_csv("lipids_canopus.csv")
colnames(lipids)[1] ="row_ID"
lipids[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(lipids$'row_ID', '_', 7)

columns_to_keep <- c(1,6) 
lipids <- lipids[, columns_to_keep]
colnames(lipids)[2] <- "most_specific_class"
lipids <- lipids[order(lipids$most_specific_class), ]

#####ALL lipids
Merged_Hits_lipids <- merge(scaled_table_t, lipids, by.x = "row_ID", by.y = "row_ID",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_lipids [,-178:-180]
data2 <- data[, -1:-2]
data2 <- data2[, -177]

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

testHeatmap2(gLogCpmData, gAnnotationData) 


#####long-chain fa
lc <- lipids[296:403,]
Merged_Hits_lc <- merge(scaled_table_t, lc, by.x = "row_ID", by.y = "row_ID",  
                            all.x = FALSE, all.y = FALSE)

data <- Merged_Hits_lc [,-177:-179]
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

testHeatmap2(gLogCpmData, gAnnotationData) 


#####linoleic acid
la <- lipids[193:253,]

Merged_Hits_la <- merge(scaled_table_t, la, by.x = "row_ID", by.y = "row_ID",  
                        all.x = FALSE, all.y = FALSE)

data <- Merged_Hits_la [,-177:-179]
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

testHeatmap2(gLogCpmData, gAnnotationData) 



#####medium-chain fa
mc <- lipids[427:470,]
Merged_Hits_mc <- merge(scaled_table_t, mc, by.x = "row_ID", by.y = "row_ID",  
                        all.x = FALSE, all.y = FALSE)

data <- Merged_Hits_mc [,-177:-179]
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

testHeatmap2(gLogCpmData, gAnnotationData) 



#####prostagl
x <- lipids[477:498,]
Merged_Hits_x <- merge(scaled_table_t, x, by.x = "row_ID", by.y = "row_ID",  
                        all.x = FALSE, all.y = FALSE)

data <- Merged_Hits_x [,-177:-179]
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

testHeatmap2(gLogCpmData, gAnnotationData) 

