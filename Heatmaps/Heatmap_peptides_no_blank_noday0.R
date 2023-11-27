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
setwd("/home/anastasiia/Desktop/mouse_studies/Heatmaps/Peptides")


library(tidyverse)
scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/CLR_Scaled.csv")
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
Peptides <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/Peptides_no_blank.csv")
Library_Hits_Refined_AK <- read_delim("Peptides_no_blank.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Library_Hits_Refined_AK <- Library_Hits_Refined_AK [-26:-27,]
Library_Hits_Refined_AK <- Library_Hits_Refined_AK [-1,]
View(Library_Hits_Refined_AK)

Merged_Hits_peptides <- merge(scaled_table_t, Library_Hits_Refined_AK, by.x = "row_ID", by.y = "Row_ID",  
                              all.x = FALSE, all.y = FALSE)
Merged_Hits_peptides_cut <- Merged_Hits_peptides [,-181:-188]
Merged_Hits_peptides_cut <- Merged_Hits_peptides_cut [,-177:-179]
#View(Merged_Hits_peptides_cut)
data <- Merged_Hits_peptides_cut %>% column_to_rownames("Compound_Name") 
data2 <- data[, -1]
View(data2)

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

ggsave("heatmap_peptides.svg", width = 8, height = 6)





