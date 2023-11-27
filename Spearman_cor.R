#https://github.com/allegra-aron/Stromatolite_analysis/blob/main/correlations_between_two_feature_tables.ipynb

#install.packages("Hmisc")
#install.packages("corrplot")
#install.packages('d3heatmap')
#install.packages("devtools")
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("rstudio/d3heatmap")
#install.packages("reshape")
#install.packages("reshape2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq")
library("DESeq")
library(reshape)
library(reshape2)
library(d3heatmap)
library(Hmisc)
library(htmlwidgets)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(santaR)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(pheatmap)
library(magrittr) 
library(tibble) 
library(paletteer)
library(tidyverse)
library(gridExtra)

setwd('/home/anastasiia/Desktop/mouse_studies/Spearman_cor/')
getwd() 

EC <- read_csv("~/Desktop/mouse_studies/Spearman_cor/EC_counts_mouse_n.csv")


sessionInfo()


##############OG


input_file <- "correlation-table-metabolotes-and-microbes.csv"
rows_ds1 <- 1:2372 # which rows (not including header) contain dataset1?
rows_ds2 <- 2373:3697 # which rows (not including header) contain dataset2?
prev_filt <- 0.1 # minimum % of samples that a feature must be observed in
norm_ds1 <- TRUE # do you want to normalize the first dataset?
norm_ds2 <- FALSE # do you want to normalize the second dataset?
scale_features <- TRUE # do you want to scale prior to correlation?
padj_meth <- "bonferroni" # method to adjust for multiple hypothesis tests (can change to "BH" to be less stringent)
p_thresh <- 0.05 #
data <- read.csv('correlation-table-metabolotes-and-microbes.csv')
if (nrow(data) == length(c(rows_ds1, rows_ds2))) {message("Looks good!")} else message("Check row numbers again")


datat <- t(data) # transpose so that samples are in rows and features are in columns
colnames(datat)<- datat[1,] # feature names are now the first row, make this row the column names
datat <- datat[-1,] # then remove it
datat <- as.data.frame(datat,stringsAsFactors=F)
datat <- as.data.frame(sapply(datat, as.numeric)) # make all values numeric
rownames(datat) <- colnames(data)[-1] # use sample names as row names



ds1 <- datat[,rows_ds1]
ds2 <- datat[,rows_ds2]
ds1_filt <- ds1[,apply(ds1, 2, function(x) {sum(x > 0) > prev_filt*nrow(ds1)})]
ds2_filt <- ds2[,apply(ds2, 2, function(x) {sum(x > 0) > prev_filt*nrow(ds2)})]

if (norm_ds1) {
  ds1_norm <- t(apply(ds1_filt, 1, function(x) {x/sum(x)})) # for each sample (row) divide each feature by the sum of all features
} else ds1_norm <- ds1_filt

if (norm_ds2) {
  ds2_norm <- t(apply(ds2_filt, 1, function(x) {x/sum(x)})) 
} else ds2_norm <- ds2_filt

if (scale_features) {
  ds1_norm <- scale(ds1_norm)
  ds2_norm <- scale(ds2_norm)
}

cor_mat <- Hmisc::rcorr(x = as(ds1_norm, "matrix"), y = as(ds2_norm, "matrix"), type = "spearman") # spearman recommended for microbiome data
# get correlations
cor_r <- cor_mat$r[1:ncol(ds1_norm), -c(1:ncol(ds1_norm))] # removing ds1~ds1 and ds2~ds2 correlations, so we only keep ds1~ds2 correlations
# get pvalues
cor_p <- cor_mat$P[1:ncol(ds1_norm), -c(1:ncol(ds1_norm))]

head(cor_r)
head(cor_p)
pheatmap(t(cor_r),  clustering_method="ward.D", clustering_distance_cols="canberra",
         show_colnames = TRUE,show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 1, 
         filename = "corr_coef_all_AT.pdf")

map <- d3heatmap(t(cor_r), distfun=function(x) dist(x, method="canberra"), 
                 hclustfun=function(x) hclust(x, method="ward.D"),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
saveWidget(map, "corr_coef_all_AT.html")


hold_r <- melt(cor_r) # putting all values in one column
hold_p <- melt(cor_p)
hold_p$padj<-p.adjust(hold_p$value, method = padj_meth) # adjust pvals for multiple hypothesis tests

hold_r_sig <- hold_r[hold_p$padj < p_thresh,] 
hold_p_sig <- hold_p[hold_p$padj < p_thresh,]

# new matrix with only significant correlations (insignificant correlations are set to 0)
corr_r_sig <- dcast(hold_r_sig, Var1~Var2, fill = 0)
row.names(corr_r_sig) <- corr_r_sig$Var1 # re-assign row names as feature names then remove that column
corr_r_sig <- corr_r_sig[,-1]


pheatmap(t(corr_r_sig),  clustering_method="ward.D", clustering_distance_cols="canberra",
         show_colnames = TRUE,show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 1, 
         filename = "corr_coef_sig_AT.pdf")

map <- d3heatmap(t(corr_r_sig), distfun=function(x) dist(x, method="canberra"), 
                 hclustfun=function(x) hclust(x, method="ward.D"),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
saveWidget(map, "corr_coef_sig.html") 

colnames(hold_r)[3] <- "corr_coef" # fix column names
colnames(hold_p)[3] <- "pval"
to_save <- cbind(hold_r, hold_p[,3:4]) # join coefs and pvals
write.csv(to_save, "corr_coef_all.csv", row.names = FALSE) # save table

# save only correlations with significant pval
write.csv(subset(to_save, to_save$pval < p_thresh), "corr_coef_sig.csv", row.names = FALSE) 

hist(-log10(to_save$pval),breaks=100)

hist(to_save$corr_coef,breaks=100)
dev.copy(png,'histogram_sig.png')
dev.off()

hist(to_save$corr_coef,-1:1, breaks=100)
h_1 <- hist(to_save$corr_coef,-1:1, breaks=100)
h_1$breaks
h_1$counts

Cstack_info()



















##############Normalised

Normalised_day96 <- read_csv("Normalised_day96.csv")
Normalised <- Normalised_day96 %>% column_to_rownames ('ID') %>% t() %>% as.data.frame %>% rownames_to_column("mouse_number")
merge<- merge(EC, Normalised, by = "mouse_number")
data <- merge %>% column_to_rownames ('mouse_number') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")
rows_ds1 <- 1:1023# which rows (not including header) contain dataset1?
rows_ds2 <- 1024:6753 # which rows (not including header) contain dataset2?
prev_filt <- 0.1 # minimum % of samples that a feature must be observed in
norm_ds1 <- TRUE # do you want to normalize the first dataset?
norm_ds2 <- FALSE # do you want to normalize the second dataset?
scale_features <- TRUE # do you want to scale prior to correlation?
padj_meth <- "BH" # method to adjust for multiple hypothesis tests (can change to "BH" to be less stringent)
p_thresh <- 0.05 # alpha for adjusted pvalues

if (nrow(data) == length(c(rows_ds1, rows_ds2))) {message("Looks good!")} else message("Check row numbers again")

datat <- t(data) # transpose so that samples are in rows and features are in columns
colnames(datat)<- datat[1,] # feature names are now the first row, make this row the column names
datat <- datat[-1,] # then remove it
datat <- as.data.frame(datat,stringsAsFactors=F)
datat <- as.data.frame(sapply(datat, as.numeric)) # make all values numeric
rownames(datat) <- colnames(data)[-1] # use sample names as row names

ds1 <- datat[,rows_ds1]
ds2 <- datat[,rows_ds2]

ds1_filt <- ds1[,apply(ds1, 2, function(x) {sum(x > 0) > prev_filt*nrow(ds1)})]
ds2_filt <- ds2[,apply(ds2, 2, function(x) {sum(x > 0) > prev_filt*nrow(ds2)})]

if (norm_ds1) {
  ds1_norm <- t(apply(ds1_filt, 1, function(x) {x/sum(x)})) # for each sample (row) divide each feature by the sum of all features
} else ds1_norm <- ds1_filt

if (norm_ds2) {
  ds2_norm <- t(apply(ds2_filt, 1, function(x) {x/sum(x)})) 
} else ds2_norm <- ds2_filt

if (scale_features) {
  ds1_norm <- scale(ds1_norm)
  ds2_norm <- scale(ds2_norm)
}

cor_mat <- Hmisc::rcorr(x = as(ds1_norm, "matrix"), y = as(ds2_norm, "matrix"), type = "pearson") # spearman recommended for microbiome data
# get correlations
cor_r <- cor_mat$r[1:ncol(ds1_norm), -c(1:ncol(ds1_norm))] # removing ds1~ds1 and ds2~ds2 correlations, so we only keep ds1~ds2 correlations
# get pvalues
cor_p <- cor_mat$P[1:ncol(ds1_norm), -c(1:ncol(ds1_norm))]

head(cor_r)
head(cor_p)
pheatmap(t(cor_r),  clustering_method="ward.D", clustering_distance_cols="canberra",
         show_colnames = TRUE,show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 1, 
         filename = "corr_coef_all.pdf")

map <- d3heatmap(t(cor_r), distfun=function(x) dist(x, method="canberra"), 
                 hclustfun=function(x) hclust(x, method="ward.D"),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
saveWidget(map, "corr_coef_all.html")


hold_r <- melt(cor_r) # putting all values in one column
#reshape2::melt(cor_r)
hold_p <- melt(cor_p)
hold_p$padj<-p.adjust(hold_p$value, method = padj_meth) # adjust pvals for multiple hypothesis tests

hold_r_sig <- hold_r[hold_p$padj < p_thresh,] 
hold_p_sig <- hold_p[hold_p$padj < p_thresh,]

# new matrix with only significant correlations (insignificant correlations are set to 0)
corr_r_sig <- dcast(hold_r_sig, Var1~Var2, fill = 0)
row.names(corr_r_sig) <- corr_r_sig$Var1 # re-assign row names as feature names then remove that column
corr_r_sig <- corr_r_sig[,-1]


pheatmap(t(corr_r_sig),  clustering_method="ward.D", clustering_distance_cols="canberra",
         show_colnames = TRUE,show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 1, 
         filename = "corr_coef_sig.pdf")

map <- d3heatmap(t(corr_r_sig), distfun=function(x) dist(x, method="canberra"), 
                 hclustfun=function(x) hclust(x, method="ward.D"),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
saveWidget(map, "corr_coef_sig.html") 

colnames(hold_r)[3] <- "corr_coef" # fix column names
colnames(hold_p)[3] <- "pval"
to_save <- cbind(hold_r, hold_p[,3:4]) # join coefs and pvals
write.csv(to_save, "corr_coef_all.csv", row.names = FALSE) # save table

# save only correlations with significant pval
write.csv(subset(to_save, to_save$pval < p_thresh), "corr_coef_sig.csv", row.names = FALSE) 

hist(-log10(to_save$pval),breaks=100)

hist(to_save$corr_coef,breaks=100)
dev.copy(png,'histogram_sig.png')
dev.off()

hist(to_save$corr_coef,-1:1, breaks=100)
h_1 <- hist(to_save$corr_coef,-1:1, breaks=100)
h_1$breaks
h_1$counts

Cstack_info()

###############CLR scaled
CLR <- read_csv("~/Desktop/mouse_studies/Spearman_cor/CLR_Scaled_day96.csv")
merge<- merge(EC, CLR, by = "mouse_number")
data <- merge %>% column_to_rownames ('mouse_number') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")

rows_ds1 <- 1:1023# which rows (not including header) contain dataset1?
rows_ds2 <- 1024:2046 # which rows (not including header) contain dataset2?
prev_filt <- 0.1 # minimum % of samples that a feature must be observed in
norm_ds1 <- TRUE # do you want to normalize the first dataset?
norm_ds2 <- FALSE # do you want to normalize the second dataset?
scale_features <- FALSE # do you want to scale prior to correlation?
padj_meth <- "bonferroni" # method to adjust for multiple hypothesis tests (can change to "BH" to be less stringent)
p_thresh <- 0.05 # alpha for adjusted pvalues

if (nrow(data) == length(c(rows_ds1, rows_ds2))) {message("Looks good!")} else message("Check row numbers again")

datat <- t(data) # transpose so that samples are in rows and features are in columns
colnames(datat)<- datat[1,] # feature names are now the first row, make this row the column names
datat <- datat[-1,] # then remove it
datat <- as.data.frame(datat,stringsAsFactors=F)
datat <- as.data.frame(sapply(datat, as.numeric)) # make all values numeric
rownames(datat) <- colnames(data)[-1] # use sample names as row names

ds1 <- datat[,rows_ds1]
ds2 <- datat[,rows_ds2]

ds1_filt <- ds1[,apply(ds1, 2, function(x) {sum(x > 0) > prev_filt*nrow(ds1)})]
ds2_filt <- ds2[,apply(ds2, 2, function(x) {sum(x > 0) > prev_filt*nrow(ds2)})]
  
if (norm_ds1) {
  ds1_norm <- t(apply(ds1_filt, 1, function(x) {x/sum(x)})) # for each sample (row) divide each feature by the sum of all features
} else ds1_norm <- ds1_filt

if (norm_ds2) {
  ds2_norm <- t(apply(ds2_filt, 1, function(x) {x/sum(x)})) 
} else ds2_norm <- ds2_filt

if (scale_features) {
  ds1_norm <- scale(ds1_norm)
  ds2_norm <- scale(ds2_norm)
}

cor_mat <- Hmisc::rcorr(x = as(ds1_norm, "matrix"), y = as(ds2_norm, "matrix"), type = "spearman") # spearman recommended for microbiome data
# get correlations
cor_r <- cor_mat$r[1:ncol(ds1_norm), -c(1:ncol(ds1_norm))] # removing ds1~ds1 and ds2~ds2 correlations, so we only keep ds1~ds2 correlations
# get pvalues
cor_p <- cor_mat$P[1:ncol(ds1_norm), -c(1:ncol(ds1_norm))]

head(cor_r)
head(cor_p)
pheatmap(t(cor_r),  clustering_method="ward.D", clustering_distance_cols="canberra",
         show_colnames = TRUE,show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 1, 
         filename = "corr_coef_all.pdf")

map <- d3heatmap(t(cor_r), distfun=function(x) dist(x, method="canberra"), 
                 hclustfun=function(x) hclust(x, method="ward.D"),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
saveWidget(map, "corr_coef_all.html")

hold_r <- melt(cor_r) # putting all values in one column
hold_p <- melt(cor_p)
hold_p$padj<-p.adjust(hold_p$value, method = padj_meth) # adjust pvals for multiple hypothesis tests

hold_r_sig <- hold_r[hold_p$padj < p_thresh,] 
hold_p_sig <- hold_p[hold_p$padj < p_thresh,]

# new matrix with only significant correlations (insignificant correlations are set to 0)
corr_r_sig <- dcast(hold_r_sig, Var1~Var2, fill = 0)
row.names(corr_r_sig) <- corr_r_sig$Var1 # re-assign row names as feature names then remove that column
corr_r_sig <- corr_r_sig[,-1]

map <- d3heatmap(t(corr_r_sig), distfun=function(x) dist(x, method="canberra"), 
                 hclustfun=function(x) hclust(x, method="ward.D"),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
saveWidget(map, "corr_coef_sig.html") 

colnames(hold_r)[3] <- "corr_coef" # fix column names
colnames(hold_p)[3] <- "pval"
to_save <- cbind(hold_r, hold_p[,3:4]) # join coefs and pvals
write.csv(to_save, "corr_coef_all.csv", row.names = FALSE) # save table

# save only correlations with significant pval
write.csv(subset(to_save, to_save$pval < p_thresh), "corr_coef_sig.csv", row.names = FALSE) 

hist(-log10(to_save$pval),breaks=100)

hist(to_save$corr_coef,breaks=100)
dev.copy(png,'histogram_sig.png')
dev.off()

hist(to_save$corr_coef,-1:1, breaks=100)
h_1 <- hist(to_save$corr_coef,-1:1, breaks=100)
h_1$breaks
h_1$counts





