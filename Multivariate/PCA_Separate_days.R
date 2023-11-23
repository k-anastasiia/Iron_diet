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

setwd('~/Desktop/GitHub/PCA')
getwd()

metadata <- read_csv("md_new.csv")
metadata <- metadata[-1:-10,-1] #exclude bedding and blanks
colnames(metadata)[1] <- "...1"
metadata$ATTRIBUTE_Study_Day<-as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata<- metadata[order(metadata$ATTRIBUTE_Study_Day), ]
scaled_table <- read_csv("CLR_Scaled.csv")
#Day 0
metadata_day_0 <-metadata[1:12,]
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
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

#Day 4
metadata_day_4 <-metadata[13:26,]
data_merge_4 <- metadata_day_4 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_4 <- data_merge_4 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_4 <- fviz_pca_ind(
  res.pca_4,
  col.ind = data_merge_4$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_4 + ggtitle("PCA (Day 4)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines

#Day 8
metadata_day_8 <-metadata[27:38,]
data_merge_8 <- metadata_day_8 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_8 <- data_merge_8 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_8 <- fviz_pca_ind(
  res.pca_8,
  col.ind = data_merge_8$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_8 + ggtitle("PCA (Day 8)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines

#Day 12
metadata_day_12 <-metadata[39:53,]
data_merge_12 <- metadata_day_12 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_12 <- data_merge_12 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_12 <- fviz_pca_ind(
  res.pca_12,
  col.ind = data_merge_12$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_12 + ggtitle("PCA (Day 12)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines

#Day 19
metadata_day_19 <-metadata[54:68,]
data_merge_19 <- metadata_day_19 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_19 <- data_merge_19 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_19<- fviz_pca_ind(
  res.pca_19,
  col.ind = data_merge_19$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_19 + ggtitle("PCA (Day 19)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines

#Day 26
metadata_day_26 <-metadata[69:79,]
data_merge_26 <- metadata_day_26 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_26 <- data_merge_26 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_26 <- fviz_pca_ind(
  res.pca_26,
  col.ind = data_merge_26$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_26 + ggtitle("PCA (Day 26)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines

#Day 33
metadata_day_33 <-metadata[80:94,]
data_merge_33 <- metadata_day_33 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_33 <- data_merge_33 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_33 <- fviz_pca_ind(
  res.pca_33,
  col.ind = data_merge_33$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_33 + ggtitle("PCA (Day 33)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines


#Day 47
metadata_day_47 <-metadata[95:107,]
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
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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


#Day 54
metadata_day_54 <-metadata[108:119,]
data_merge_54 <- metadata_day_54 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_54 <- data_merge_54 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_54 <- fviz_pca_ind(
  res.pca_54,
  col.ind = data_merge_54$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_54 + ggtitle("PCA (Day 54)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none")   # Customize the axis lines


#Day 61
metadata_day_61 <-metadata[120:134,]
data_merge_61 <- metadata_day_61 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_61 <- data_merge_61 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_61 <- fviz_pca_ind(
  res.pca_61,
  col.ind = data_merge_61$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_61 + ggtitle("PCA (Day 61)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none") 


#Day 68
metadata_day_68 <-metadata[135:146,]
data_merge_68 <- metadata_day_68 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_68 <- data_merge_68 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_68 <- fviz_pca_ind(
  res.pca_68,
  col.ind = data_merge_68$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_68 + ggtitle("PCA (Day 68)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none") 

#Day 75
metadata_day_75 <-metadata[147:161,]
data_merge_75 <- metadata_day_75 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca <- data_merge_75 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) 
pca_plot_75 <- fviz_pca_ind(
  res.pca,
  col.ind = data_merge_75$ATTRIBUTE_Diet,
  addEllipses = T,
  repel = F,
  #legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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


#Day 82
metadata_day_82 <-metadata[162:172,]
data_merge_82 <- metadata_day_82 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca_82 <- data_merge_82 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before

pca_plot_82 <- fviz_pca_ind(
  res.pca_82,
  col.ind = data_merge_82$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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

pca_plot_82 + ggtitle("PCA (Day 82)") +
  theme(panel.grid = element_blank(),  # Remove all grid lines
        axis.line = element_line(color = "black"),  # Customize the axis lines
        legend.position = "none") 



#Day 96
metadata_day_96 <-metadata[173:187,]
data_merge_96 <- metadata_day_96 %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")


res.pca_96 <- data_merge_96 %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) 

pca_plot_96 <- fviz_pca_ind(
  res.pca_96,
  col.ind = data_merge_96$ATTRIBUTE_Diet,
  addEllipses = TRUE,
  repel = TRUE,
  legend.title = "Diet",
  palette = c("Normal" = "green4", "Deficient" = "red3", "0.6%_Fe" = "royalblue1"),
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


#PLOT
pca_plot_0 <- pca_plot_0 + ggtitle("Day 0")+  theme(legend.position = "none")
pca_plot_4 <- pca_plot_4 + ggtitle("Day 4")+  theme(legend.position = "none")
pca_plot_8 <- pca_plot_8 + ggtitle("Day 8")+  theme(legend.position = "none")
pca_plot_12 <- pca_plot_12 + ggtitle("Day 12")+  theme(legend.position = "none")
pca_plot_19 <- pca_plot_19 + ggtitle("Day 19")+  theme(legend.position = "none")
pca_plot_26 <- pca_plot_26 + ggtitle("Day 26")+  theme(legend.position = "none")
pca_plot_33 <- pca_plot_33 + ggtitle("Day 33")+  theme(legend.position = "none")
pca_plot_47 <- pca_plot_47 + ggtitle("Day 47")+  theme(legend.position = "none")
pca_plot_54 <- pca_plot_54 + ggtitle("Day 54")+  theme(legend.position = "none")
pca_plot_61 <- pca_plot_61 + ggtitle("Day 61")+  theme(legend.position = "none")
pca_plot_68 <- pca_plot_68 + ggtitle("Day 68")+  theme(legend.position = "none")
pca_plot_75 <- pca_plot_75 + ggtitle("Day 75")+  theme(legend.position = "none")
pca_plot_82 <- pca_plot_82 + ggtitle("Day 82")+  theme(legend.position = "none")
pca_plot_96 <- pca_plot_96 + ggtitle("Day 96")+  theme(legend.position = "none")


grid.arrange(pca_plot_0, pca_plot_4, pca_plot_8, pca_plot_12,pca_plot_19, pca_plot_26, 
             pca_plot_33, pca_plot_47, pca_plot_54, pca_plot_61,pca_plot_68, pca_plot_75, pca_plot_82, pca_plot_96,  ncol = 3)


####################PERMANOVA############################
library(vegan)

distm_0 <- dist(data_merge_0, method = "euclidean")# compute distance
adonres_0 <- adonis2(distm_0 ~ data_merge_0[,colnames(data_merge_0) == 'ATTRIBUTE_Diet'])
View(adonres_0)

distm_4 <- dist(data_merge_4, method = "euclidean")# compute distance
adonres_4 <- adonis2(distm_4 ~ data_merge_4[,colnames(data_merge_4) == 'ATTRIBUTE_Diet'])
View(adonres_4)

distm_8 <- dist(data_merge_8, method = "euclidean")# compute distance
adonres_8 <- adonis2(distm_8 ~ data_merge_8[,colnames(data_merge_8) == 'ATTRIBUTE_Diet'])
View(adonres_8)

distm_12 <- dist(data_merge_12, method = "euclidean")# compute distance
adonres_12 <- adonis2(distm_12 ~ data_merge_12[,colnames(data_merge_12) == 'ATTRIBUTE_Diet'])
View(adonres_12)

distm_19 <- dist(data_merge_19, method = "euclidean")# compute distance
adonres_19 <- adonis2(distm_19 ~ data_merge_19[,colnames(data_merge_19) == 'ATTRIBUTE_Diet'])
View(adonres_19)

distm_26 <- dist(data_merge_26, method = "euclidean")# compute distance
adonres_26 <- adonis2(distm_26 ~ data_merge_26[,colnames(data_merge_26) == 'ATTRIBUTE_Diet'])
View(adonres_26)

distm_33 <- dist(data_merge_33, method = "euclidean")# compute distance
adonres_33 <- adonis2(distm_33 ~ data_merge_33[,colnames(data_merge_33) == 'ATTRIBUTE_Diet'])
View(adonres_33)

distm_47 <- dist(data_merge_47, method = "euclidean")# compute distance
adonres_47 <- adonis2(distm_47 ~ data_merge_47[,colnames(data_merge_47) == 'ATTRIBUTE_Diet'])
View(adonres_47)

distm_54 <- dist(data_merge_54, method = "euclidean")# compute distance
adonres_54 <- adonis2(distm_54 ~ data_merge_54[,colnames(data_merge_54) == 'ATTRIBUTE_Diet'])
View(adonres_54)

distm_61 <- dist(data_merge_61, method = "euclidean")# compute distance
adonres_61 <- adonis2(distm_61 ~ data_merge_61[,colnames(data_merge_61) == 'ATTRIBUTE_Diet'])
View(adonres_61)

distm_68 <- dist(data_merge_68, method = "euclidean")# compute distance
adonres_68 <- adonis2(distm_68 ~ data_merge_68[,colnames(data_merge_68) == 'ATTRIBUTE_Diet'])
View(adonres_68)

distm_75 <- dist(data_merge_75, method = "euclidean")# compute distance
adonres_75 <- adonis2(distm_75 ~ data_merge_75[,colnames(data_merge_75) == 'ATTRIBUTE_Diet'])
View(adonres_75)

distm_82 <- dist(data_merge_82, method = "euclidean")# compute distance
adonres_82 <- adonis2(distm_82 ~ data_merge_82[,colnames(data_merge_82) == 'ATTRIBUTE_Diet'])
View(adonres_82)

distm_96 <- dist(data_merge_96, method = "euclidean")# compute distance
adonres_96 <- adonis2(distm ~ data_merge_96[,colnames(data_merge_96) == 'ATTRIBUTE_Diet'])
View(adonres_96)

#PLOT
Permanova <- read_csv("Permanova.csv")
Permanova<-Permanova[,-4:-5 ]
Permanova<-Permanova[,-2]

plot(Permanova$Day, Permanova$F, 
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
points(Permanova$Day, Permanova$F, pch = 19, col = "black", cex = 1.5)

# Connect the data points with a dashed line
lines(Permanova$Day,Permanova$F, type = "l", lty = 2, col = "black")



