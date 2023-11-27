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
getwd() #to see the working directory 
#Setting a working directory in Google Colab:
setwd("/home/anastasiia/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance")


scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/CLR_Scaled.csv")
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
purine <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/purine.csv")
colnames(purine)[1] ="row_ID"
purine[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(purine$'row_ID', '_', 7)
columns_to_keep <- c(2,3,4,10) 
purine <- purine[, columns_to_keep]



Merged_Hits_purine <- merge(scaled_table_t, purine, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_purine [,-177:-180]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)



############purine_NPC_class

purine <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/purine_NPC_class.csv")
colnames(purine)[1] ="row_ID"
purine[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(purine$'row_ID', '_', 7)
columns_to_keep <- c(2,3,4,10) 
purine <- purine[, columns_to_keep]



Merged_Hits_purine <- merge(scaled_table_t, purine, by.x = "row_ID", by.y = "row_ID2",  
                            all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_purine [,-177:-180]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)


############purine_classyFire_most_specific_class.csv

purine <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/purine_classyFire_most_specific_class.csv")
colnames(purine)[1] ="row_ID"
purine[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(purine$'row_ID', '_', 7)
columns_to_keep <- c(2,3,4,10) 
purine <- purine[, columns_to_keep]



Merged_Hits_purine <- merge(scaled_table_t, purine, by.x = "row_ID", by.y = "row_ID2",  
                            all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_purine [,-177:-180]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)



############purine_ClassyFire_class.csv

purine <- read_csv("~/Desktop/mouse_studies/Heatmaps/BA_conjugates/AA_BA_abundance/purine_ClassyFire_class.csv")
colnames(purine)[1] ="row_ID"
purine[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(purine$'row_ID', '_', 7)
columns_to_keep <- c(2,3,4,10) 
purine <- purine[, columns_to_keep]



Merged_Hits_purine <- merge(scaled_table_t, purine, by.x = "row_ID", by.y = "row_ID2",  
                            all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_purine [,-177:-180]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)

































###2
library(tidyverse)
scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/CLR_Scaled.csv")


#change the column name
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")

# Split name column into several names
colnames(scaled_table_t)[1] ="row_ID"
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'row_ID', '_', 3)
#View(scaled_table_t)

#remove day 0
columns_to_exclude <- c(2,9,18,27,36,87,120,129,139,160,170,178) 
scaled_table_t <- scaled_table_t[, -columns_to_exclude]


#Merge Tables
Peptides <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/canopus_peptides.csv")
colnames(Peptides)[1] ="row_ID"
Peptides[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(Peptides$'row_ID', '_', 7)
columns_to_keep <- c(5,14) 
Peptides <- Peptides[, columns_to_keep]
colnames(Peptides)[1] <- "Compound_Name"
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Dipeptides <- Peptides[70:249,]

Merged_Hits_peptides <- merge(scaled_table_t, Dipeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_peptides [,-177:-178]
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
#legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)

#3
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Cpeptides <- Peptides[10:64,]
#Cpeptides <- Peptides[-291:-743,] #linear+cyclic
#Cpeptides <- Cpeptides[-65:-252,]
#Cpeptides <- Cpeptides[-1:-10,]
Merged_Hits_peptides <- merge(scaled_table_t, Cpeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_peptides [,-177:-178]
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


#4
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Lpeptides <- Peptides[253:290,]

Merged_Hits_peptides <- merge(scaled_table_t, Lpeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_peptides [,-177:-178]
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

#5
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Tripeptides <- Peptides[-1:-326,]

Merged_Hits_peptides <- merge(scaled_table_t, Tripeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)


data <- Merged_Hits_peptides [,-177:-178]
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


####ClassyFire#most specific class
Peptides <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/canopus_peptides.csv")
colnames(Peptides)[1] ="row_ID"
Peptides[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(Peptides$'row_ID', '_', 7)
columns_to_keep2 <- c(6,14) 
Peptides <- Peptides[, columns_to_keep2]
colnames(Peptides)[1] <- "Compound_Name"
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Dipeptides <- Peptides[1:132,]

Merged_Hits_peptides <- merge(scaled_table_t, Dipeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-178]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 

#Oligo
Opeptides <- Peptides[145:549,]

Merged_Hits_peptides <- merge(scaled_table_t, Opeptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-178]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 
#Peptides-dont work
Peptides <- Peptides[-1:-549,]

Merged_Hits_peptides <- merge(scaled_table_t, Peptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-178]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 

###tripeptides
#Merge Tables
Peptides <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/canopus_peptides.csv")
colnames(Peptides)[1] ="row_ID"
Peptides[c('row_ID', 'm/z', 'RT','a', 'c', 'b','row_ID2')] <- str_split_fixed(Peptides$'row_ID', '_', 7)
columns_to_keep <- c(5,6,14) 
Peptides <- Peptides[, columns_to_keep]
colnames(Peptides)[2] <- "Compound_Name"
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Peptides <- Peptides[1:549,]

columns_to_keep <- c(1,3) 
Peptides <- Peptides[, columns_to_keep]
colnames(Peptides)[1] <- "Compound_Name"
Peptides <- Peptides[order(Peptides$Compound_Name), ]
Peptides <- Peptides[229:549,]


Merged_Hits_peptides <- merge(scaled_table_t, Peptides, by.x = "row_ID", by.y = "row_ID2",  
                              all.x = FALSE, all.y = FALSE)
data <- Merged_Hits_peptides [,-177:-178]
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
  heatmap.2(logCPM, margins = c(6,12), ColSideColors = sampleColors, breaks = seq(-10, 10, 0.2), col = bluered(100),  
            key.xlab = "log(feature abundance)", keysize = 1.2, scale = "none", key = TRUE, symkey = FALSE,
            symbreaks = FALSE, density.info = "none", trace = "none",
            cexCol = 0.5, cexRow = 1.35) 
}

testHeatmap2(gLogCpmData, gAnnotationData) 

