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
setwd("/home/anastasiia/Desktop/mouse_studies/Internal_standard_drift")


i_s <- read_csv("~/Desktop/mouse_studies/Internal_standard_drift/Sulfadimethoxine.csv")
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
# Line plot
plot(x = seq_along(i_s1), y = i_s1[1, ], type = "l", xlab = "Samles", ylab = "Intensity", xaxt = "n")
axis(side = 1, at = seq_along(colnames), labels = colnames, tick = TRUE, las = 2)

# Scatter plot
plot(x = seq_along(i_s1), y = i_s1[1, ], type = "p", xlab = "Samples", ylab = "Intensity", xaxt = "n")
axis(side = 1, at = seq_along(colnames), labels = colnames, tick = TRUE, las = 2)



i_s2<- i_s [-2:-3,]
# Line plot
plot(x = seq_along(i_s2), y = i_s2[1, ], type = "l", xlab = "Samles", ylab = "Intensity", xaxt = "n")
axis(side = 1, at = seq_along(sorted_colnames), labels = sorted_colnames, tick = TRUE, las = 2)

# Scatter plot
plot(x = seq_along(i_s2), y = i_s2[1, ], type = "p", xlab = "Samples", ylab = "Intensity", xaxt = "n")
axis(side = 1, at = seq_along(sorted_colnames), labels = sorted_colnames, tick = TRUE, las = 2)



i_s_all<- i_s [-1:-2,]
# Line plot
plot(x = seq_along(i_s_all), y = i_s_all[1, ], type = "l", xlab = "Samles", ylab = "Intensity", xaxt = "n")
axis(side = 1, at = seq_along(sorted_colnames), labels = sorted_colnames, tick = TRUE, las = 2)

# Scatter plot
plot(x = seq_along(i_s_all), y = i_s_all[1, ], type = "p", xlab = "Samples", ylab = "Intensity", xaxt = "n")
axis(side = 1, at = seq_along(sorted_colnames), labels = sorted_colnames, tick = TRUE, las = 2)





df <- data.frame(x = colnames(i_s), y1 = i_s[,1], y2 = i_s[,2], y3 = i_s[,3])

# Convert the dataframe to a tidy format
df_tidy <- tidyr::gather(df, key = "group", value = "value", -x)

# Plot the data using ggplot2
ggplot(df_tidy, aes(x = x, y = value, color = group)) +
  geom_point()




