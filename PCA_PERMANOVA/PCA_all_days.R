library(readr)
library(factoextra)
library(tidyverse) #used for data science. The eight core packages inside this library are: ggplot2 (data visualisation), dplyr (data manipulation), tidyr, readr, purrr, tibble, stringr, and forcats
library(KODAMA) # to use the normalisation function
library(ggrepel) #mainly used to repel overlapping text labels in ggplots
library(vegan) #popular library for analysing ecological diversity and for multivariate analysis of community data. Here, we use it for PCoA
library(IRdisplay) #better display of output cells in Jupyter Notebooks running iwth IRKernel. Library not needed when running the script in RStudio
library(svglite) #to save the plots in support vector graphics (svg) format
library(factoextra) #for extracting and visualizing outputs of multivariate analyses such as PCA, k-means
library(ggsci) #provides color palettes for ggplot2 that can be used for scientific journals
library(matrixStats) #contains highly optimized functions to perform statistics on matrix data
library(cowplot) #efficient functions to arrange several plots
library(ComplexHeatmap) #for visualising heatmaps
library(dendextend) # for getting dendograms
library(NbClust) # for finding the optimum no.of clusters to be used for a clustering method
library(tidyverse)

setwd('~/Desktop/GitHub/PCA')
getwd()


metadata<- read_csv("~/Desktop/GitHub/PCA/md_new.csv")
metadata <- metadata[-1:-10,-1] #exclude bedding and blanks
colnames(metadata)[1] <- "...1"
scaled_table <- read_csv("~/Desktop/GitHub/PCA/CLR_scaled.csv")
data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
   column_to_rownames("...1")

res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
#fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE,repel = TRUE, legend.title = "Diet",
             #palette = c("Normal"="green4", "Deficient"="red3","0.6%_Fe"="royalblue1"))
fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Diet", palette = c("Normal"="green4", "Deficient"="red3", "0.6%_Fe"="royalblue1"),
             geom.ind = "point", pointsize = 2)


#COLORING BY MOUSE NUMBER

library(factoextra)
data_merge <- metadata %>% dplyr::select("...1", "ATTIBUTE_mouse_no") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")
data_merge$ATTIBUTE_mouse_no<-as.numeric(data_merge$ATTIBUTE_mouse_no)
res.pca <- data_merge %>% dplyr::select(-ATTIBUTE_mouse_no) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)


color_palette <-c("1"="#498915", "2"="#88bd4b", "3"="#a6cd78", "4"="#c3dea5", "5"="#e1eed2",
                  "6"="#990000", "7"="#cc0000", "8"="#ff0000", "9"="#ff6666", "10"="#ffb2b2",
                  "11"='#2B3A66', "12"= '#3E5A9F', "13"="#a3bbec", "14"="#c2d2f2", "15"="#e0e8f9")


shape_mapping <- c("1"=19, "2"=19, "3"=19, "4"=19, "5"=19,
                   "6"=19, "7"=19, "8"=19, "9"=19, "10"=19,
                   "11"=19, "12"=19, "13"=19, "14"=19, "15"=19)

# Create a data frame with PCA results
pca_data <- data.frame(PC1 = res.pca$x[,1], PC2 = res.pca$x[,2], ATTIBUTE_mouse_no = data_merge$ATTIBUTE_mouse_no)

# Create a PCA plot with colored points and different shapes (no labels)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(ATTIBUTE_mouse_no), shape = factor(ATTIBUTE_mouse_no))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_mapping) +
  labs(x = "PC1 (25.6%)", y = "PC2 (15.3%)", title = "PCA") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    panel.grid.major = element_blank(),          # Remove major grid lines
    panel.grid.minor = element_blank(),          # Remove minor grid lines
    axis.text = element_text(size = 12),         # Increase the size of axis tick labels
    axis.ticks = element_line(color = "black")   # Add axis ticks
  )

# Display the PCA plot with ticks
print(pca_plot)



#COLORING BY DAY


data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Study_Day") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")
data_merge$ATTRIBUTE_Study_Day<- as.numeric(data_merge$ATTRIBUTE_Study_Day)
res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_Study_Day) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)

# Create a custom color palette with 15 shades of red to blue
color_palette <- colorRampPalette(c("red", "blue"))(15)

# Create a shape mapping with shape 19 for all values
shape_mapping <- rep(19, length(unique(data_merge$ATTRIBUTE_Study_Day)))

# Create a data frame with PCA results
pca_data <- data.frame(PC1 = res.pca$x[,1], PC2 = res.pca$x[,2], ATTRIBUTE_Study_Day = data_merge$ATTRIBUTE_Study_Day)

# Create a PCA plot with colored labels and shape 19
pca_plot1 <- ggplot(pca_data, aes(x = PC1, y = PC2, label = ATTRIBUTE_Study_Day, color = factor(ATTRIBUTE_Study_Day), shape = factor(ATTRIBUTE_Study_Day))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_mapping) +
  labs(x = "PC1 (25.6%)", y = "PC2 (15.3%)", title = "PCA") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    panel.grid.major = element_blank(),          # Remove major grid lines
    panel.grid.minor = element_blank(),          # Remove minor grid lines
    axis.text = element_text(size = 12),         # Increase the size of axis tick labels
    axis.ticks = element_line(color = "black")   # Add axis ticks
  )

# Display the PCA plot with ticks
print(pca_plot1)

pca_plot <- pca_plot + ggtitle("Mouse ID")+  theme(legend.position = "none")
pca_plot1 <- pca_plot1 + ggtitle("Study day")+  theme(legend.position = "none")
grid.arrange(pca_plot, pca_plot1,  ncol = 2)

#PERMANOVA
metadata<- read_csv("~/Desktop/GitHub/PCA/md_new.csv")
metadata <- metadata[-1:-10,-1] #exclude bedding and blanks
colnames(metadata)[1] <- "...1"
scaled_table <- read_csv("~/Desktop/GitHub/PCA/CLR_scaled.csv")

data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

distm <- dist(data_merge, method = "euclidean")# compute distance
adonres <- adonis2(distm ~ data_merge[,colnames(data_merge) == 'ATTRIBUTE_Diet'])
View(adonres)



