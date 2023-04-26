#calling the necessary packages:
library(ggplot2)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(gplots)
library(pheatmap)
library(heatmap3)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(clr)

#set working directory
Directory <- setwd('/home/anastasiia/Desktop/mouse_studies/Heatmaps/File_processing')

Blank_removed <- read_csv("Blank_removed_0.3.csv")
colnames(Blank_removed)[1] <- "feature"
Blank_removed[c('row ID', 'row m/z', 'row retention time')] <- str_split_fixed(Blank_removed$'feature', '_', 3)
Blank_removed <- Blank_removed[, c((ncol(Blank_removed)-2):ncol(Blank_removed), 1:(ncol(Blank_removed)-3))]
Blank_removed_ft <- Blank_removed[,-4]
# Specify columns to modify
cols_to_modify <- c(4:190)
# Add letter "a" at the end of specified column names
names(Blank_removed_ft)[cols_to_modify] <- paste0(names(Blank_removed_ft)[cols_to_modify], " Peak area")
view(Blank_removed_ft)
write.csv(Blank_removed_ft, "~/Desktop/mouse_studies/Blank_removed_ft.csv",row.names = TRUE) #delete first numerical column after saving 

metadata <- read_csv("~/Desktop/mouse_studies/File_processing/gnps_metadata_for_mouse_studies_01112023.csv")
metadata_no_blanks <- metadata[-1:-10,] 
write.csv(metadata_no_blanks, "~/Desktop/mouse_studies/metadata_no_blanks.csv",row.names = TRUE) #delete first numerical column after saving 
