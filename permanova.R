library(readr)
library(factoextra)
library(tidyverse) 
library(KODAMA)
library(ggrepel) 
library(svglite) #to save the plots in support vector graphics (svg) format
library(ggsci)
library(NbClust) 
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)

setwd("/home/anastasiia/Desktop/mouse_studies/PCAs/PERMANOVA")
getwd() 

permanova<- read_csv("~/Desktop/mouse_studies/PCAs/PERMANOVA/Permanova.csv")
permanova_b<-permanova[1:8, ]
permanova_a<-permanova[-1:-8, ]

#ALL
plot(permanova$Day, permanova$F, 
     xlab = "Day", 
     ylab = "F",
     xaxt = "n",  # This removes both labels and ticks from the x-axis
     type = "n"  , # This prevents plotting points initially
     ylim = c(1, 10)
)

# Specify custom x-axis tick positions
custom_ticks <- c(0, 4, 8, 12, 19, 26, 33, 47, 54, 61, 68, 75, 82, 96)  # Specify the specific tick positions you want

# Set the custom x-axis ticks
axis(1, at = custom_ticks, labels = custom_ticks)
# Specify custom y-axis tick positions
custom_ticks_y <- 1:10  # Specify y-axis tick positions from 1 to 10

# Customize the data points (color and thickness)
points(permanova$Day, permanova$F, pch = 19, col = "black", cex = 1.5)

# Connect the data points with a dashed line
lines(permanova$Day, permanova$F, type = "l", lty = 2, col = "black")





#BEFORE
plot(permanova_b$Day, permanova_b$F, 
     xlab = "Day", 
     ylab = "F",
     xaxt = "n",  # This removes both labels and ticks from the x-axis
     type = "n"   # This prevents plotting points initially
)

# Specify custom x-axis tick positions
custom_ticks <- c(0, 4, 8, 12, 19, 26, 33, 47)  # Specify the specific tick positions you want

# Set the custom x-axis ticks
axis(1, at = custom_ticks, labels = custom_ticks)

# Customize the data points (color and thickness)
points(permanova_b$Day, permanova_b$F, pch = 19, col = "black", cex = 1.5)

# Connect the data points with a dashed line
lines(permanova_b$Day, permanova_b$F, type = "l", lty = 2, col = "black")




#AFTER
plot(permanova_a$Day, permanova_a$F, 
     xlab = "Day", 
     ylab = "F",
     yaxt = "n",
     xaxt = "n",  # This removes both labels and ticks from the x-axis
     type = "n",  # This prevents plotting points initially
     ylim = c(1, 10)  # Set the y-axis limits to 1 to 10
)


# Specify custom x-axis tick positions
custom_ticks <- c(54, 61, 68, 75, 82, 96)  # Specify the specific tick positions you want

# Set the custom x-axis ticks
axis(1, at = custom_ticks, labels = custom_ticks)

# Customize the data points (color and thickness)
points(permanova_a$Day, permanova_a$F, pch = 19, col = "black", cex = 1.5)

# Connect the data points with a dashed line
lines(permanova_a$Day, permanova_a$F, type = "l", lty = 2, col = "black")
# Specify custom y-axis tick positions
custom_ticks_y <- 1:10  # Specify y-axis tick positions from 1 to 10

# Set the custom y-axis ticks
axis(2, at = custom_ticks_y)



#BEFORE

# Create a line plot without x-axis labels and ticks, and set y-axis limits
plot(permanova_b$Day, permanova_b$F, 
     xlab = "Day", 
     ylab = "F",
     yaxt = "n",
     xaxt = "n",  # This removes both labels and ticks from the x-axis
     type = "n",  # This prevents plotting points initially
     ylim = c(1, 10)  # Set the y-axis limits to 1 to 10
)

# Specify custom x-axis tick positions
custom_ticks_x <- c(0, 4, 8, 12, 19, 26, 33, 47)  # Specify the specific tick positions you want

# Set the custom x-axis ticks
axis(1, at = custom_ticks_x, labels = custom_ticks_x)

# Customize the data points (color and thickness)
points(permanova_b$Day, permanova_b$F, pch = 19, col = "black", cex = 1.5)

# Connect the data points with a dashed line
lines(permanova_b$Day, permanova_b$F, type = "l", lty = 2, col = "black")

# Specify custom y-axis tick positions
custom_ticks_y <- 1:10  # Specify y-axis tick positions from 1 to 10

# Set the custom y-axis ticks
axis(2, at = custom_ticks_y)













# Create a scatter plot
plot(permanova$Day, permanova$F, 
     xlab = "Day", 
     ylab = "F",
     xaxt = "n"  # This removes both labels and ticks from the x-axis
)

axis(1, at = permanova$Day, labels = permanova$Day)


# Create a line plot
plot(permanova$Day, permanova$F, 
     xlab = "Day", 
     ylab = "F",
     xaxt = "n"  # This removes both labels and ticks from the x-axis
)

# Set custom x-axis ticks
axis(1, at = permanova$Day, labels = permanova$Day)




