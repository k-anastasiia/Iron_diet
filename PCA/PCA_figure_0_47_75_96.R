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
library(factoextra)
library(gridExtra)
library(factoextra)
library(ggplot2)

getwd() #to see the working directory 
#Setting a working directory:
setwd("/home/anastasiia/Desktop/mouse_studies/PCAs")

#Day 75
metadata_day75_ak <- read_csv("~/Desktop/mouse_studies/PCAs/Day_75/metadata_day75_ak.csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Day_75/CLR_Scaled.csv")
data_merge <- metadata_day75_ak %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  dplyr::filter(`...1` != "186.mzML") %>%
  column_to_rownames("...1")



data_merge <- metadata_day75_ak %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
   column_to_rownames("...1")


res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) 
pca_plot_75 <- fviz_pca_ind(
  res.pca,
  col.ind = data_merge$ATTRIBUTE_Diet,
  addEllipses = T,
  repel = F,
  #legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6% Fe" = "royalblue1"),
  geom.ind = "point",
  ellipse.alpha = 0,
  pointsize = 3,
  ggtheme = theme_minimal()+
    theme(
      axis.line = element_line(color = "black"),  # Add axis lines
      panel.grid.major = element_blank(),          # Remove major grid lines
      panel.grid.minor = element_blank(),          # Remove minor grid lines
      axis.text = element_text(size = 12),         # Increase the size of axis tick labels
      axis.ticks = element_line(color = "black")   # Add axis ticks
  ))

pca_plot_75 +
  labs(title = "PCA (Day 75)") +  # Set the title to an empty string
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")  # Remove legends

#Day 0
metadata_day_0 <- read_csv("~/Desktop/mouse_studies/PCAs/Day_0/metadata_day_0.csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Day_0/CLR_Scaled.csv")
data_merge_0 <- metadata_day_0 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_0 <- data_merge_0 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_0 <- fviz_pca_ind(
  res.pca_0,
  col.ind = data_merge_0$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6% Fe" = "royalblue1"),
  geom.ind = "point",
  ellipse.alpha = 0,
  pointsize = 3,
  ggtheme = theme_minimal()+
    theme(
      axis.line = element_line(color = "black"),  # Add axis lines
      panel.grid.major = element_blank(),          # Remove major grid lines
      panel.grid.minor = element_blank(),          # Remove minor grid lines
      axis.text = element_text(size = 12),         # Increase the size of axis tick labels
      axis.ticks = element_line(color = "black")   # Add axis ticks
    ))

pca_plot_0 + ggtitle("PCA (Day 0)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines




#Day 47
metadata_day_47 <- read_csv("~/Desktop/mouse_studies/PCAs/Day_47/metadata_day47only.csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Day_47/CLR_Scaled.csv")


data_merge_47 <- metadata_day_47 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_47 <- data_merge_47 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_47 <- fviz_pca_ind(
  res.pca_47,
  col.ind = data_merge_47$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6% Fe" = "royalblue1"),
  geom.ind = "point",
  ellipse.alpha = 0,
  pointsize = 3,
  ggtheme = theme_minimal()+
    theme(
      axis.line = element_line(color = "black"),  # Add axis lines
      panel.grid.major = element_blank(),          # Remove major grid lines
      panel.grid.minor = element_blank(),          # Remove minor grid lines
      axis.text = element_text(size = 12),         # Increase the size of axis tick labels
      axis.ticks = element_line(color = "black")   # Add axis ticks
    ))

pca_plot_47 + ggtitle("PCA (Day 47)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none") 
#Day 96
metadata_day96 <- read_csv("~/Desktop/mouse_studies/PCAs/Day_96/metadata_day96.csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Day_96/CLR_Scaled.csv")
#data_merge_96 <- metadata_day96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  #left_join(scaled_table) %>% 
#  dplyr::filter(`...1` != "222.mzML") %>%
 # column_to_rownames("...1")

data_merge_96 <- metadata_day96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


res.pca_96 <- data_merge_96 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) 

pca_plot_96 <- fviz_pca_ind(
  res.pca_96,
  col.ind = data_merge_96$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6% Fe" = "royalblue1"),
  geom.ind = "point",
  ellipse.alpha = 0,
  pointsize = 3,
  ggtheme = theme_minimal()+
    theme(
      axis.line = element_line(color = "black"),  # Add axis lines
      panel.grid.major = element_blank(),          # Remove major grid lines
      panel.grid.minor = element_blank(),          # Remove minor grid lines
      axis.text = element_text(size = 12),         # Increase the size of axis tick labels
      axis.ticks = element_line(color = "black")   # Add axis ticks
    ))

pca_plot_96 + ggtitle("PCA (Day 96)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none") 







pca_plot_0 <- pca_plot_0 + ggtitle("PCA (Day 0)")
pca_plot_47 <- pca_plot_47 + ggtitle("PCA (Day 47)")
pca_plot_75 <- pca_plot_75 + ggtitle("PCA (Day 75)")
pca_plot_96 <- pca_plot_96 + ggtitle("PCA (Day 96)")

grid.arrange(pca_plot_0, pca_plot_47, pca_plot_75, pca_plot_96, ncol = 2)
grid.arrange(pca_plot_0, pca_plot_47, pca_plot_75, pca_plot_96, ncol = 4)


