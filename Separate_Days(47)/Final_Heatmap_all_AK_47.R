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
library("RColorBrewer")
library(paletteer)

getwd() #to see the working directory 
#Setting a working directory in Google Colab:
setwd("/home/anastasiia/Desktop/mouse_studies/Heatmaps/Day_47")

#Imputed_Scaled_QuantTable_day47only <- read_csv("~/Desktop/mouse_studies/Heatmaps/Day_47/Imputed_Scaled_QuantTable_day47only.csv")
Imputed_Scaled_QuantTable_day47only <- read_csv("~/Desktop/mouse_studies/PCAs/Days/Day_47/CLR_Scaled.csv")
View(Imputed_Scaled_QuantTable_day47only)

#change the column name
colnames(Imputed_Scaled_QuantTable_day47only)[1] <- "row_ID"
Imputed_Scaled_QuantTable_day47only[,1]
View(Imputed_Scaled_QuantTable_day47only)
#other other way of transposing
#Imputed_Scaled_QuantTable_day75_t = setNames(data.frame(t(Imputed_Scaled_QuantTable_day75[,-1])),Imputed_Scaled_QuantTable_day75[,1])

#other way of transposing
names <- Imputed_Scaled_QuantTable_day47only[,1]
Imputed_Scaled_QuantTable_day47only_t <- as.data.frame(as.matrix(t(Imputed_Scaled_QuantTable_day47only[,-1])))
names_t <- transpose(names)
#View(names_t)
colnames(Imputed_Scaled_QuantTable_day47only_t) <- unlist(names_t)
write.csv(Imputed_Scaled_QuantTable_day47only_t, "~/Desktop/mouse_studies/Heatmaps/Day_47/Imputed_Scaled_Quant_47_T.csv",row.names = TRUE)

Imputed_Scaled_Quant_47_T <- read_csv("Imputed_Scaled_Quant_47_T.csv")

# Split name column into several names
colnames(Imputed_Scaled_Quant_47_T)[1] ="row_ID"
Imputed_Scaled_Quant_47_T[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(Imputed_Scaled_Quant_47_T$row_ID, '_', 3)
View(Imputed_Scaled_Quant_47_T)

#Merge Tables
Library_Hits_Refined_AK <- read_delim("Library_Hits_Refined_AK_new_master.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
View(Library_Hits_Refined_AK)
Merged_Hits_47 <- merge(Imputed_Scaled_Quant_47_T, Library_Hits_Refined_AK, by.x = "row_ID", by.y = "Row_ID",
                        all.x = FALSE, all.y = FALSE)
View(Merged_Hits_47)

write.csv(Merged_Hits_47, "~/Desktop/mouse_studies/Heatmaps/Day_47/Merged_Hits_47.csv",row.names = TRUE)

#if columns have different names, then 
#merge(df1, df2, by.x = "df1ColName", by.y = "df2ColName").
#if no need to keep unmatched rows from the table y, but keep from the table x, then 
#all.x = TRUE, all.y = FALSE

#cut the table
Merged_Hits_47_cut <- Merged_Hits_47[, (names(Merged_Hits_47) %in% c("Compound_Name", "121.mzML", "122.mzML", "123.mzML","124.mzML","125.mzML","126.mzML", "128.mzML","129.mzML","130.mzML","131.mzML","132.mzML","133.mzML","135.mzML"))]
write.csv(Merged_Hits_47_cut, "~/Desktop/mouse_studies/Heatmaps/Day_47/Merged_Hits_47_cut.csv",row.names = TRUE)
#rename the rownames
data <- read_csv("Merged_Hits_47_cut2_no121.csv") %>% column_to_rownames("Compound_Name") #cut2-same as table cut but with no repeating compound names. Make sure column names are still 1st row
View(data)

#delete the first column
data2 <- data[, -1]
View(data2)

# Load the example "data"
gLogCpmData = as.matrix(data2)
print(gLogCpmData)


# Load the example annotation/metadata
gAnnotationData = read.table("md_47_no121.txt", sep = '\t', header = TRUE,
                             fill = TRUE)
gAnnotationData

# Make helper function to map metadata category to color
mapDietToColor<-function(annotations){
  colorsVector = ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="0.6% Fe", 
                        "blue", ifelse(gAnnotationData["ATTRIBUTE_Diet"]=="Normal", 
                        "green", "red"))
                        return(colorsVector)
}
# Test heatmap with column annotations
testHeatmap2<-function(logCPM, annotations) {    
  sampleColors = mapDietToColor(annotations)
  heatmap.2(logCPM, margins=c(8,12), ColSideColors=sampleColors, breaks=seq(-6,10,0.16), col=bluered(100),  
            key.xlab="feature abundance", keysize=1.2, scale="none", key=TRUE, symkey=FALSE, symbreaks=FALSE, density.info="none", trace="none") 
}        
testHeatmap2(gLogCpmData, gAnnotationData)
legend(x=locator(1), legend=c("0.6% Fe","Normal", "Deficient"), fill=c("blue", "green", "red"), cex=.7)

