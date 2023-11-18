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
setwd("/home/anastasiia/Desktop/mouse_studies/PCAs")


####################SIMONE'S CODE################################################
metadata<- read_csv("~/Desktop/mouse_studies/PCAs/metadata_full (ak).csv")
metadata <- metadata[-188,] #exclude bedding
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/old/Days/Day_47/CLR_Scaled.csv")
data_merge <- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
   column_to_rownames("...1")

res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
#fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE,repel = TRUE, legend.title = "Diet",
             #palette = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Diet", palette = c("Normal"="green4", "Deficient"="red3", "0.6% Fe"="royalblue1"),
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
  labs(x = "PC1 (26.3%)", y = "PC2 (14.8%)", title = "PCA") +
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
  labs(x = "PC1 (26.3%)", y = "PC2 (14.8%)", title = "PCA") +
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



















pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
view(pc_loadings)


pc_loadings[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings$'feature', '_', 3)
pc1_all <- pc_loadings[,-190:-191]
pc1_all <- pc1_all [,-3:-188]
pc1_all <- pc1_all[, c("row_ID", setdiff(names(pc1_all), "row_ID"))]
pc1_all <- pc1_all [,-2]
write.csv(pc1_all, "~/Desktop/mouse_studies/PCAs/pc1_all.csv",row.names = TRUE)

fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
####################PERMANOVA############################
distm <- dist(data_merge, method = "euclidean")# compute distance

adonres <- adonis2(distm ~ data_merge[,colnames(data_merge) == 'ATTRIBUTE_Diet'])
View(adonres)

distm <- dist(data_merge_after, method = "euclidean")# compute distance

adonres_a <- adonis2(distm ~ data_merge_after[,colnames(data_merge_after) == 'ATTRIBUTE_Diet'])
View(adonres_a)

distm <- dist(data_merge_before, method = "euclidean")# compute distance

adonres_b <- adonis2(distm ~ data_merge_before[,colnames(data_merge_before) == 'ATTRIBUTE_Diet'])
View(adonres_b)





####################BEFORE################################################
metadata$ATTRIBUTE_Study_Day <- as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata <- metadata[order(metadata$ATTRIBUTE_Study_Day),]
metadata_before <- metadata[14:107,] 
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Days/Day_47/CLR_Scaled.csv")
data_merge_before <- metadata_before %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca <- data_merge_before %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = data_merge_before$ATTRIBUTE_Diet, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Diet", palette = c("Normal"="green4", "Deficient"="red3", "0.6% Fe"="royalblue1"),
             geom.ind = "point", pointsize = 2)


pc_loadings_before <- res.pca$rotation
pc_loadings_before <- as_tibble(pc_loadings_before, rownames = "feature")
view(pc_loadings_before)

fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
####################PERMANOVA_BEFORE############################
distm_before <- dist(data_merge_before, method = "euclidean")# compute distance

adonres_before <- adonis2(distm ~ data_merge_before[,colnames(data_merge_before) == 'ATTRIBUTE_Diet'])
View(adonres_before)



####################AFTER################################################
metadata$ATTRIBUTE_Study_Day <- as.numeric(metadata$ATTRIBUTE_Study_Day)
metadata <- metadata[order(metadata$ATTRIBUTE_Study_Day),]
metadata_after <- metadata[108:187,] 
scaled_table <- read_csv("~/Desktop/mouse_studies/PCAs/Days/Day_47/CLR_Scaled.csv")
data_merge_after  <- metadata_after  %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1")

res.pca <- data_merge_after  %>% dplyr::select(-ATTRIBUTE_Diet) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = data_merge_after $ATTRIBUTE_Diet, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Diet", palette = c("Normal"="green4", "Deficient"="red3", "0.6% Fe"="royalblue1"),
             geom.ind = "point", pointsize = 2)

pc_loadings_after <- res.pca$rotation
pc_loadings_after <- as_tibble(pc_loadings_after, rownames = "feature")
view(pc_loadings_after)

fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE     # Avoid text overlapping
)


####################PERMANOVA_AFTER############################
distm_after  <- dist(data_merge_after, method = "euclidean")# compute distance

adonres_after  <- adonis2(distm ~ data_merge_after [,colnames(data_merge_after) == 'ATTRIBUTE_Diet'])
View(adonres_after )







#############EXTRA
library(ggplot2)

# Create a custom color palette
color_palette <- c("1"="green4", "2"="green4", "3"="green4", "4"="green4", "5"="green4",
                   "6"="red3", "7"="red3", "8"="red3", "9"="red3", "10"="red3",
                   "11"="royalblue1", "12"="royalblue1", "13"="royalblue1", "14"="royalblue1", "15"="royalblue1")

# Create a custom shape mapping
shape_mapping <- c("1"=19, "2"=15, "3"=17, "4"=18, "5"=5,
                   "6"=19, "7"=15, "8"=17, "9"=18, "10"=5,
                   "11"=19, "12"=15, "13"=17, "14"=18, "15"=5)


# Create a data frame with PCA results
pca_data <- data.frame(PC1 = res.pca$x[,1], PC2 = res.pca$x[,2], ATTIBUTE_mouse_no = data_merge$ATTIBUTE_mouse_no)

# Create a PCA plot with colored labels and different shapes
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = ATTIBUTE_mouse_no, color = factor(ATTIBUTE_mouse_no), shape = factor(ATTIBUTE_mouse_no))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_mapping) +
  labs(x = "PC1", y = "PC2", title = "PCA Plot") +
  theme_minimal()

# Display the PCA plot
pca_plot


library(ggplot2)



data_merge$ATTIBUTE_mouse_no<- as.numeric(data_merge$ATTIBUTE_mouse_no)
# Create a custom color palette
color_palette <-c("1"="#498915", "2"="#88bd4b", "3"="#a6cd78", "4"="#c3dea5", "5"="#e1eed2",
                  "6"="#990000", "7"="#cc0000", "8"="#ff0000", "9"="#ff6666", "10"="#ffb2b2",
                  "11"='#2B3A66', "12"= '#3E5A9F', "13"="#a3bbec", "14"="#c2d2f2", "15"="#e0e8f9")




# Create a custom shape mapping
shape_mapping <- c("1"=19, "2"=15, "3"=7, "4"=17, "5"=5,
                   "6"=19, "7"=15, "8"=7, "9"=17, "10"=5,
                   "11"=19, "12"=15, "13"=7, "14"=17, "15"=5)

# Create a data frame with PCA results
pca_data <- data.frame(PC1 = res.pca$x[,1], PC2 = res.pca$x[,2], ATTIBUTE_mouse_no = data_merge$ATTIBUTE_mouse_no)

# Create a PCA plot with colored labels and different shapes
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(ATTIBUTE_mouse_no), shape = factor(ATTIBUTE_mouse_no))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_mapping) +
  geom_text(aes(label = ATTIBUTE_mouse_no), hjust = 1, vjust = -1) + # Add labels directly to the plot
  labs(x = "PC1", y = "PC2", title = "PCA Plot") +
  theme_minimal()

# Display the PCA plot
pca_plot



# Create a custom color palette

color_palette <-c("1"="#498915", "2"="#88bd4b", "3"="#a6cd78", "4"="#c3dea5", "5"="#e1eed2",
                  "6"="#990000", "7"="#cc0000", "8"="#ff0000", "9"="#ff6666", "10"="#ffb2b2",
                  "11"='#2B3A66', "12"= '#3E5A9F', "13"="#a3bbec", "14"="#c2d2f2", "15"="#e0e8f9")


#color_palette <- c("1"="#98FB98", "2"='#478778', "3"="green4", "4"="#4CBB17", "5"="#00A36C",
#   "6"="#FF6961", "7"="#DC143C", "8"="red3", "9"="#8B0000", "10"="#E0115F",
# "11"='#6495ED', "12"= '#7B68EE', "13"="#4169E1", "14"="#000080", "15"="#0000CD")

# Create a custom shape mapping
shape_mapping <- c("1"=19, "2"=15, "3"=7, "4"=17, "5"=5,
                   "6"=19, "7"=15, "8"=7, "9"=17, "10"=5,
                   "11"=19, "12"=15, "13"=7, "14"=17, "15"=5)

# Create a data frame with PCA results
pca_data <- data.frame(PC1 = res.pca$x[,1], PC2 = res.pca$x[,2], ATTIBUTE_mouse_no = data_merge$ATTIBUTE_mouse_no)

# Create a PCA plot with colored points and different shapes (no labels)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(ATTIBUTE_mouse_no), shape = factor(ATTIBUTE_mouse_no))) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_mapping) +
  labs(x = "PC1", y = "PC2", title = "PCA Plot") +
  theme_minimal()

# Display the PCA plot without labels
pca_plot











# Define start and end colors for three different color sets
start_color_set1 <- "#498915"
  end_color_set1 <- "#88bd4b"
    
  start_color_set2 <- "#990000"
    end_color_set2 <- "#ff0000"
      
    start_color_set3 <- "#2B3A66"
      end_color_set3 <- "#a3bbec"
        
      # Number of colors per set
      n_colors_per_set <- 5
      
      # Create custom color palettes for each set
      color_palette_set1 <- colorRampPalette(c(start_color_set1, end_color_set1))(n_colors_per_set)
      color_palette_set2 <- colorRampPalette(c(start_color_set2, end_color_set2))(n_colors_per_set)
      color_palette_set3 <- colorRampPalette(c(start_color_set3, end_color_set3))(n_colors_per_set)
      
      # Combine the color palettes from different sets
      color_palette <- c(
        color_palette_set1,
        color_palette_set2,
        color_palette_set3
      )
      
      # Define different shapes for each set
      shape_mapping <- c(
        rep(19, n_colors_per_set),
        rep(17, n_colors_per_set),
        rep(15, n_colors_per_set)
      )
      
      # Create a data frame with PCA results
      pca_data <- data.frame(
        PC1 = res.pca$x[, 1],
        PC2 = res.pca$x[, 2],
        ATTIBUTE_mouse_no = data_merge$ATTIBUTE_mouse_no
      )
      
      # Create a PCA plot with colored points and different shapes (no labels)
      pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(ATTIBUTE_mouse_no), shape = factor(ATTIBUTE_mouse_no))) +
        geom_point(size = 3) +
        scale_color_manual(values = color_palette) +
        scale_shape_manual(values = shape_mapping) +
        labs(x = "PC1", y = "PC2", title = "PCA Plot") +
        theme_minimal()
      
      # Display the PCA plot without labels
      pca_plot
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      