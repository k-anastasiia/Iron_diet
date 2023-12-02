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

setwd('~/Desktop/GitHub/Heatmaps/Purine_nucleosides')
getwd()

scaled_table <- read_csv("CLR_Scaled.csv")
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'SampleID', '_', 3)

#remove day 0
columns_to_exclude <- c(2,9,18,27,36,87,120,129,139,160,170,178) 
scaled_table_t <- scaled_table_t[, -columns_to_exclude]

#Merge Tables
purine <- read_csv("purine.csv")
colnames(purine)[1] ="row_ID"
purine[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(purine$'row_ID', '_', 7)
columns_to_keep <- c(2,3,4,10) 
purine <- purine[, columns_to_keep]

Merged_Hits_purine <- merge(scaled_table_t, purine, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_purine [,-178:-183]
data <- data[, -2]
data2 <- data[, -1]
gLogCpmData = as.matrix(data2)

# Load the example annotation/metadata
gAnnotationData = read.table("md_full_noday0.txt", sep = '\t', header = TRUE,
                             fill = TRUE)
gAnnotationData

# Make helper function to map metadata category to color
mapDietToColor<-function(annotations){
  colorsVector = ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="0.6% Fe", 
                        "blue", ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="Normal", 
                        "green", "red"))
                        return(colorsVector)
}

# Test heatmap with larger font size for column labels
testHeatmap2 <- function(logCPM, annotations) {    
  sampleColors <- mapDietToColor(annotations)
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
