# https://github.com/Functional-Metabolomics-Lab/Statistical-analysis-of-non-targeted-LC-MSMS-data/blob/main/Stats_Untargeted_Metabolomics.ipynb

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
library(ggrepel)
options(install.packages.compile.from.source="never")
if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("tidyverse","factoextra","KODAMA","vegan","IRdisplay","svglite")
#Global settings for plot size in the output cell:
options(repr.plot.width=10, repr.plot.height=10,res=600) # the parameters: width, height & resolution can be changed
pacman::p_load("matrixStats","ggsci","FSA","cowplot") #installing & loading the necessary libraries


setwd("/home/anastasiia/Desktop/mouse_studies/Univariate/day_47")
getwd() #to see the working directory 

metadata<- read_csv("~/Desktop/mouse_studies/Univariate/day_47/metadata_4_75.csv")
metadata<- metadata[83:95,]
scaled_table <- read_csv("~/Desktop/mouse_studies/Univariate/day_47/CLR_Scaled.csv")
data_merge1<- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1") 

data_merge <- data_merge1[,-1]

###################Test for normality #####################
colnames(data_merge)
norm_test <- data_merge %>% dplyr::select(.,1)  #gets the 1st column 

head(norm_test)
hist(norm_test[,1], xlab=names(norm_test[1]), main="Histogram")

qqnorm(norm_test[,1])
qqline(norm_test[,1], col = "blue") #add the line with theoretic normal distribution

shapiro.test(norm_test[,1])


data_merge2 <- data_merge1[order(metadata$ATTRIBUTE_Diet),]
#NORMAL
data_merge_normal <- data_merge2[9:13,] 
data_merge_normal <- data_merge_normal[,-1] 
op_shapiro <- data_merge_normal  %>% colnames() %>% as.data.frame #getting all the column names into a dataframe called "op_shapiro"
colnames(op_shapiro)[1] <- "Metabolites" #naming the column as "Metabolites"
op_shapiro$p <- sapply(1:ncol(data_merge_normal ), function(n){shapiro.test(data_merge_normal[,n])$p}) #performing shapiro test on each feature and getti
op_shapiro$p_adj <- p.adjust(op_shapiro$p, method="fdr") #adding a column "p_adj" with corrected p-values. The method used is
op_shapiro$distribution <- ifelse(op_shapiro$p_adj<0.05,"Non-normal","Normal")
op_shapiro %>% arrange(p) %>% head(n=5) #getting the top-5 p-value metabolites
paste("No.of features with normal distribution:", sum(op_shapiro$distribution=="Normal"))
paste("No.of features with non-normal distribution:", sum(op_shapiro$distribution=="Non-normal"))

#DEFICIENT
data_merge_deficient <- data_merge2[5:8,] 
data_merge_deficient <- data_merge_deficient[,-1] 
op_shapiro <- data_merge_deficient  %>% colnames() %>% as.data.frame #getting all the column names into a dataframe called "op_shapiro"
colnames(op_shapiro)[1] <- "Metabolites" #naming the column as "Metabolites"
op_shapiro$p <- sapply(1:ncol(data_merge_deficient ), function(n){shapiro.test(data_merge_deficient[,n])$p}) #performing shapiro test on each feature and getti
op_shapiro$p_adj <- p.adjust(op_shapiro$p, method="fdr") #adding a column "p_adj" with corrected p-values. The method used is
op_shapiro$distribution <- ifelse(op_shapiro$p_adj<0.05,"Non-normal","Normal")
op_shapiro %>% arrange(p) %>% head(n=5) #getting the top-5 p-value metabolites
paste("No.of features with normal distribution:", sum(op_shapiro$distribution=="Normal"))
paste("No.of features with non-normal distribution:", sum(op_shapiro$distribution=="Non-normal"))

#0.6%FE
data_merge_0.6 <- data_merge2[1:4,] 
data_merge_0.6 <- data_merge_0.6[,-1] 
op_shapiro <- data_merge_0.6 %>% colnames() %>% as.data.frame #getting all the column names into a dataframe called "op_shapiro"
colnames(op_shapiro)[1] <- "Metabolites" #naming the column as "Metabolites"
op_shapiro$p <- sapply(1:ncol(data_merge_0.6 ), function(n){shapiro.test(data_merge_0.6[,n])$p}) #performing shapiro test on each feature and getti
op_shapiro$p_adj <- p.adjust(op_shapiro$p, method="fdr") #adding a column "p_adj" with corrected p-values. The method used is
op_shapiro$distribution <- ifelse(op_shapiro$p_adj<0.05,"Non-normal","Normal")
op_shapiro %>% arrange(p) %>% head(n=5) #getting the top-5 p-value metabolites
paste("No.of features with normal distribution:", sum(op_shapiro$distribution=="Normal"))
paste("No.of features with non-normal distribution:", sum(op_shapiro$distribution=="Non-normal"))

#ALL
op_shapiro <- data_merge %>% colnames() %>% as.data.frame #getting all the column names into a dataframe called "op_shapiro"
colnames(op_shapiro)[1] <- "Metabolites" #naming the column as "Metabolites"
op_shapiro$p <- sapply(1:ncol(data_merge), function(n){shapiro.test(data_merge[,n])$p}) #performing shapiro test on each feature and getti
op_shapiro$p_adj <- p.adjust(op_shapiro$p, method="fdr") #adding a column "p_adj" with corrected p-values. The method used is
op_shapiro$distribution <- ifelse(op_shapiro$p_adj<0.05,"Non-normal","Normal")
op_shapiro %>% arrange(p) %>% head(n=5) #getting the top-5 p-value metabolites
paste("No.of features with normal distribution:", sum(op_shapiro$distribution=="Normal"))
paste("No.of features with non-normal distribution:", sum(op_shapiro$distribution=="Non-normal"))


###################ANOVA#####################

head(data_merge)

data_merge1$'ATTRIBUTE_Diet' <- as.factor(data_merge1$'ATTRIBUTE_Diet') # Convert "ATTRIBUTE_Diet" to factor
broom::tidy(aov(data_merge[,1]~data_merge1$'ATTRIBUTE_Diet')) #tidy summarizes the anova output in a  tibble

anova_out=list() # creates an empty list
anova_model=list() # creates an empty list

for(i in 1:ncol(data_merge)) {
  anova_model[[i]] <- aov(data_merge[,i] ~ data_merge1$'ATTRIBUTE_Diet') #performing ANOVA and storing the model as a list element in anova_model
  anova_tidy <- broom::tidy(anova_model[[i]]) # summarizing the model output as a tibble using 'tidy'
  result <- anova_tidy[1, ] # getting only the 1st row the tibble
  anova_out <- bind_rows(anova_out, result) # combining all anova outputs
}
anova_out$term <- colnames(data_merge) ##for each row in column'term', give the corresponding name of its feature

anova_out <- arrange(anova_out, p.value) # arranging anova_out by p values
anova_out["p_bonferroni"] <- p.adjust(anova_out$p.value,method="bonferroni") #adding another with bonferroni corrected p-values
anova_out["significant"]  <- ifelse(anova_out$p_bonferroni<0.05,"Significant","Non-significant")

names(anova_model) <- colnames(data_merge) #providing names for each list element

view(head(anova_out,20))
dim(anova_out)

anova_sig <- filter(anova_out, significant=="Significant")
anova_sig_names <-anova_sig$term #getting the names of all significant metabolites

print(paste("Total no.of features on which ANOVA test was performed:",nrow(anova_out)))
print(paste("No.of features that showed significant difference:",sum(anova_out$significant=="Significant")))
print(paste("No.of features that did not show significant difference:",sum(anova_out$significant=="Non-significant")))

#plot ANOVA results
plot_anova <- ggplot(anova_out,aes(x = log(statistic,base=10), y = -log(p.value,base=10), color = significant)) + 
  geom_point() + 
  theme_classic() + 
  scale_color_jama() + 
  ylab("-log(p)") + 
  xlab("log(F)") +
  ggrepel::geom_text_repel(data=head(anova_out), # showing the top 6 features
                           aes(label=term),
                           size=4,
                           show.legend = FALSE, max.overlaps = 100) +
  theme(panel.background = element_blank(), #making the ggplot background white
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) # giving a rectangle black border

plot_anova

###################Tukey's post-hoc test#####################

model_1 <- anova_model[[anova_sig_names[1]]] #looking at one of the anova model that showed significant difference
broom::tidy(model_1)

model_1 <- anova_model[[anova_sig_names[1]]] #looking at one of the anova model that showed significant difference
broom::tidy(TukeyHSD(model_1)) # Perform Tukey HSD test on the model_1 and summarize the result

# Normal-Deficient
output_tukey <- NULL
for(i in anova_sig_names) {
  tukey_tidy <- broom::tidy(TukeyHSD(anova_model[[i]])) # Perform Tukey HSD test on the model & summarize the result
  
  result <- filter(tukey_tidy, contrast == 'Normal-Deficient') # subsetting that particular comparison
  result$term <- i #giving the rowname as the name of its feature
  output_tukey <- bind_rows(output_tukey, result) #combine outputs from all features
}
output_tukey["p_bonferroni"] <- p.adjust(output_tukey$adj.p.value,method="bonferroni")
output_tukey["significant"] <- ifelse(output_tukey$p_bonferroni<0.05,"Significant","Nonsignificant")
output_tukey <- arrange(output_tukey, adj.p.value) # arranging anova_out by p_bonferroni values
view(head(output_tukey,399))

print(paste("Total no.of features on which tukey test was performed:",nrow(output_tukey)))
print(paste("No.of features that showed significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Significant")))
print(paste("No.of features that did not show significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Nonsignificant")))


filtered_data <- output_tukey %>% 
  filter(grepl("^(3784|2670|1807|1692|6516|1163|3767|3954|2300|3760|5823|6661|6638|4083|3548|4081|4772|5536|
  5439|6528|6817|4771|5526|2679|5393|2402|6824|5642|4817|5670)", term))

# EXTRA
plot_tukey <- ggplot(output_tukey, aes(x = estimate, y = -log(adj.p.value, base = 10), color = significant)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("black", "red3"), labels = c("Not significant", "Significant")) +
  ylab("-log(p)") +
  ggrepel::geom_text_repel(data = filtered_data, aes(x = estimate, y = -log(adj.p.value, base = 10), label = term),
                           size = 4, show.legend = F, max.overlaps = 100, force = 1, box.padding = 0.5,
                           hjust = 1.2) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black")
  ) +
  labs(x = "estimate", y = "-log(P)") + 
  geom_point(colour = "grey", alpha = 0.5) +
  geom_point(data = filtered_data, # New layer containing data subset il_genes       
             size = 2,
             shape = 21,
             fill = "firebrick",
             colour = "black")

print(plot_tukey)

# MAIN
  plot_tukey <- ggplot(output_tukey, aes(x = estimate, y = -log(adj.p.value, base = 10), color = significant)) +
  geom_point(data = output_tukey, aes(color = significant, fill = significant), size = 4, shape = 21) +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("darkgrey", "#d8abab")) +
  theme_minimal() +
  ylab("-log(p)") +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black")
  ) +
  labs(x = "estimate", y = "-log(P)") + 
  geom_point(data = filtered_data, size = 4, shape = 21,fill = "firebrick", colour = "black")

print(plot_tukey)




















# Normal-0.6% Fe
output_tukey <- NULL
for(i in anova_sig_names) {
  tukey_tidy <- broom::tidy(TukeyHSD(anova_model[[i]])) # Perform Tukey HSD test on the model & summarize the result
  
  result <- filter(tukey_tidy, contrast == 'Normal-0.6% Fe') # subsetting that particular comparison
  result$term <- i #giving the rowname as the name of its feature
  output_tukey <- bind_rows(output_tukey, result) #combine outputs from all features
}
output_tukey["p_bonferroni"] <- p.adjust(output_tukey$adj.p.value,method="bonferroni")
output_tukey["significant"] <- ifelse(output_tukey$p_bonferroni<0.05,"Significant","Nonsignificant")
output_tukey <- arrange(output_tukey, adj.p.value) # arranging anova_out by p_bonferroni values
view(head(output_tukey,399))

print(paste("Total no.of features on which tukey test was performed:",nrow(output_tukey)))
print(paste("No.of features that showed significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Significant")))
print(paste("No.of features that did not show significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Nonsignificant")))


filtered_data <- output_tukey %>% 
  filter(grepl("^(7395|4660|1413|7664|14240|6857|6871|12251|11281|8559|9442|7083|2254|6018|
  5352|5680|7997|7018|8868|10980|11829|7373|7655|6142|4194|8401|1014|10043|5519|414)", term))

# extra
plot_tukey <- ggplot(output_tukey, aes(x = estimate, y = -log(adj.p.value, base = 10), color = significant)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("black", "royalblue1"), labels = c("Not significant", "Significant")) +
  ylab("-log(p)") +
  ggrepel::geom_text_repel(data = filtered_data, aes(x = estimate, y = -log(adj.p.value, base = 10), label = term),
                           size = 4, show.legend = F, max.overlaps = 100, force = 1, box.padding = 0.5,
                           hjust = 1.2) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black")
  ) +
  labs(x = "estimate", y = "-log(P)")

print(plot_tukey)



# MAIN
plot_tukey <- ggplot(output_tukey, aes(x = estimate, y = -log(adj.p.value, base = 10), color = significant)) +
  geom_point(data = output_tukey, aes(color = significant, fill = significant), size = 4, shape = 21) +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("darkgrey", "#99ccff")) +
  theme_minimal() +
  ylab("-log(p)") +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black")
  ) +
  labs(x = "estimate", y = "-log(P)") + 
  geom_point(data = filtered_data, size = 4, shape = 21,fill = "#0066cc", colour = "black")

print(plot_tukey)









