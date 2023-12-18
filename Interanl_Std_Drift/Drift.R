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

setwd("Desktop/GitHub/Internal_standard_drift")


i_s <- read_csv("Sulfadimethoxine.csv")
i_s <- column_to_rownames(i_s, var = "row ID")
i_s<- i_s [, -1:-12]

# Remove "Peak area" from column names
colnames(i_s) <- gsub("area", "", colnames(i_s))
colnames(i_s) <- gsub("Peak", "", colnames(i_s))
colnames(i_s) <- gsub(".mzML", "", colnames(i_s))
colnames(i_s) <- gsub("   ", "", colnames(i_s))
colnames(i_s)[(ncol(i_s)-1):ncol(i_s)] <- sub("...", "", colnames(i_s)[(ncol(i_s)-1):ncol(i_s)])
#Remove bedding?
i_s<- i_s [,-26]
i_s1<- i_s [-1,]
i_s1<- i_s1 [-2,]
# Scatter plot-1
plot(x = seq_along(i_s1), y = i_s1[1, ], type = "p", xlab = "Samples", ylab = "Intensity", xaxt = "n", ylim = c(0, 3.5e+7), main = "Intensity of the feature 6567")
axis(side = 1, at = seq_along(colnames), labels = sorted_colnames, tick = TRUE, las = 2)


# Scatter plot-2
i_s2<- i_s [-2:-3,]
plot(x = seq_along(i_s2), y = i_s2[1, ], type = "p", xlab = "Samples", ylab = "Intensity", xaxt = "n", ylim = c(0, 3.5e+7), main = "Intensity of the feature 6051")
axis(side = 1, at = seq_along(colnames), labels = sorted_colnames, tick = TRUE, las = 2)

# Scatter plot-all
i_s_all<- i_s [-1:-2,]
plot(x = seq_along(i_s_all), y = i_s_all[1, ], type = "p", xlab = "Samples", ylab = "Intensity", xaxt = "n", ylim = c(0, 3.5e+7), main = "Intensity of added features")
axis(side = 1, at = seq_along(colnames), labels = sorted_colnames, tick = TRUE, las = 2)



