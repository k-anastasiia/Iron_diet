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

getwd() #to see the working directory 
#Setting a working directory in Google Colab:
setwd("/home/anastasiia/Desktop/mouse_studies/Heatmaps/Bile_acids")


library(tidyverse)
scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/Bile_acids/CLR_Scaled.csv")
View(scaled_table)

#change the column name
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")

# Split name column into several names
colnames(scaled_table_t)[1] ="row_ID"
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'row_ID', '_', 3)
View(scaled_table_t)

#remove day 0
columns_to_exclude <- c(2,9,18,27,36,87,120,129,139,160,170,178) 
scaled_table_t <- scaled_table_t[, -columns_to_exclude]


#Merge Tables
BA <- read_csv("~/Desktop/mouse_studies/Heatmaps/Bile_acids/canopus_BA.csv")
colnames(BA)[1] ="row_ID"
BA[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(BA$'row_ID', '_', 7)
columns_to_keep <- c(2,8) 
BA <- BA[, columns_to_keep]



Merged_Hits_BA <- merge(scaled_table_t, BA, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)

      
data <- Merged_Hits_BA [,-177:-178]
#colnames(Merged_Hits_peptides_cut)[177] <- "Compound_Name"

#data <- Merged_Hits_peptides_cut %>% column_to_rownames("Compound_Name") 
data2 <- data[, -1]
data2 <- data2[, -176]
#View(data2)

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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

 testHeatmap2(gLogCpmData, gAnnotationData) 
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)



 

 #Hydroxy
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
 HBA <- BA[123:215,]
 
 Merged_Hits_HBA <- merge(scaled_table_t, HBA, by.x = "row_ID", by.y = "row_ID2",  
                          all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_HBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
   heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
             key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
             symbreaks = FALSE, density.info = "none", trace = "none",
             cexCol = 0.5, cexRow = 1.35) 
 }
 
 testHeatmap2(gLogCpmData, gAnnotationData) 
 
 #Dihydroxy
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
DBA <- BA[55:71,]

 Merged_Hits_DBA <- merge(scaled_table_t, DBA, by.x = "row_ID", by.y = "row_ID2",  
                          all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_DBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
   heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
             key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
             symbreaks = FALSE, density.info = "none", trace = "none",
             cexCol = 0.5, cexRow = 1.35) 
 }
 
 testHeatmap2(gLogCpmData, gAnnotationData) 
 
 
 #Glycinated
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
 GBA <- BA[72:122,]
 
 Merged_Hits_GBA <- merge(scaled_table_t, GBA, by.x = "row_ID", by.y = "row_ID2",  
                          all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_GBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
   heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
             key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
             symbreaks = FALSE, density.info = "none", trace = "none",
             cexCol = 0.5, cexRow = 1.35) 
 }
 
 testHeatmap2(gLogCpmData, gAnnotationData) 
 
 
 
 
 #Monohydroxy
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
 MBA <- BA[216:264,]
 
 Merged_Hits_MBA <- merge(scaled_table_t, MBA, by.x = "row_ID", by.y = "row_ID2",  
                          all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_MBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
   heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
             key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
             symbreaks = FALSE, density.info = "none", trace = "none",
             cexCol = 0.5, cexRow = 1.35) 
 }
 
 testHeatmap2(gLogCpmData, gAnnotationData) 
 
 
 
 
 
 
 #Taurinated
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
 TauBA <- BA[265:285,]
 
 Merged_Hits_TauBA <- merge(scaled_table_t, TauBA, by.x = "row_ID", by.y = "row_ID2",  
                          all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_TauBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
 
 
 
 
 #Trihydroxy
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
 TBA <- BA[286:325,]
 
 Merged_Hits_TBA <- merge(scaled_table_t, TBA, by.x = "row_ID", by.y = "row_ID2",  
                            all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_TBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
   heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
             key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
             symbreaks = FALSE, density.info = "none", trace = "none",
             cexCol = 0.5, cexRow = 1.35) 
 }
 
 testHeatmap2(gLogCpmData, gAnnotationData)  
 
 
 
 
 
 
 #simple Bile acids
 colnames(BA)[1] <- "Compound_Name"
 BA <- BA[order(BA$Compound_Name), ]
 SBA <- BA[1:54,]

 Merged_Hits_SBA <- merge(scaled_table_t, SBA, by.x = "row_ID", by.y = "row_ID2",  
                         all.x = FALSE, all.y = FALSE)
 
 
 data <- Merged_Hits_SBA [,-177:-178]
 data2 <- data[, -1]
 data2 <- data2[, -176]
 
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
   heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-12, 12, 0.24), col = bluered(100),  
             key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
             symbreaks = FALSE, density.info = "none", trace = "none",
             cexCol = 0.5, cexRow = 1.35) 
 }
 
 testHeatmap2(gLogCpmData, gAnnotationData) 
 
 