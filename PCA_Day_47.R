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

getwd() #to see the working directory 
#Setting a working directory:
setwd("/home/anastasiia/Desktop/mouse_studies/PCAs/Days/Day_47")


####################SIMONE'S CODE################################################
metadata_day47only <- read_csv("metadata_day47only.csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Days/Day_47/CLR_Scaled.csv")
data_merge <- metadata_day47only %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  dplyr::filter(`...1` != "121.mzML") %>%
  column_to_rownames("...1")

res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE,repel = TRUE, legend.title = "Diet",
             palette = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))


pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
pc_loadings

fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
####################PERMANOVA############################
distm <- dist(data_merge, method = "euclidean")# compute distance

adonres <- adonis2(distm ~ data_merge[,colnames(data_merge) == 'ATTRIBUTE_Diet'])
View(adonres)













####################OLD############################

#choose only day47 from the feature table
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Days/Day_47/CLR_Scaled.csv")
scaled_table_cut <- scaled_table[(scaled_table$...1 %in%c("121.mzML", "122.mzML", "123.mzML","124.mzML","125.mzML","126.mzML", "128.mzML","129.mzML","130.mzML","131.mzML","132.mzML","133.mzML","135.mzML")),]
scaled_table_cut<-as.matrix(scaled_table_cut[1:13,2:5731])

#metadata <-scaled_table_cut2[1:13,1]


res.pca <- prcomp(scaled_table_cut, scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)


pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
pc_loadings

fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
