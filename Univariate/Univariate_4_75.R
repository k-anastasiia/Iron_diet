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

options(install.packages.compile.from.source="never")
if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("tidyverse","factoextra","KODAMA","vegan","IRdisplay","svglite")
#Global settings for plot size in the output cell:
options(repr.plot.width=10, repr.plot.height=10,res=600) # the parameters: width, height & resolution can be changed
pacman::p_load("matrixStats","ggsci","FSA","cowplot") #installing & loading the necessary libraries


setwd("/home/anastasiia/Desktop/mouse_studies/Univariate/day4_75")
getwd() #to see the working directory 

metadata<- read_csv("~/Desktop/mouse_studies/Univariate/day4_75/metadata_4_75.csv")

scaled_table <- read_csv("~/Desktop/mouse_studies/Univariate/day4_75/CLR_Scaled.csv")
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
data_merge_normal <- data_merge2[100:149,] 
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
data_merge_deficient <- data_merge2[51:99,] 
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
data_merge_0.6 <- data_merge2[1:50,] 
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

head(anova_out,20)
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

# to export as svg
svglite("plot_ANOVA.svg", width=10, height=8, bg='white')
plot_anova
dev.off()


b <- list() # creating an empty list to store the boxplot of each metabolite as an element [[i]]
for(i in 1:4){
  b[[i]] <- ggplot(data_merge1,aes(x=ATTRIBUTE_Diet,
                            y=data_merge1[,(anova_sig_names[i])], # this will get the corresponding column from the dataframe 'Data'
                            color=ATTRIBUTE_Diet))+
    geom_boxplot() + 
    scale_color_jama() +  # specifying the 'journal-friendly' color palette
    labs(x='Diet',y=as.character(anova_sig_names[i])) + #labelling the axes
    theme_classic() + # white background theme
    geom_jitter() + #adds some random variation to each point location to avoid overfitting
    theme(axis.text.x= element_blank()) #removing the x axis labels for each boxplot in a plot
}
combined_boxplots <- plot_grid(b[[1]],b[[2]],b[[3]],b[[4]], labels = c('A','B','C','D'))
combined_boxplots

ggsave("combined_boxplots.pdf",combined_boxplots,height = 9,width = 16)

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
head(output_tukey,10)

print(paste("Total no.of features on which tukey test was performed:",nrow(output_tukey)))
print(paste("No.of features that showed significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Significant")))
print(paste("No.of features that did not show significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Nonsignificant")))

plot_tukey <- ggplot(output_tukey,aes(x=estimate,y=-log(adj.p.value,base=10),color=significant)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_jama() + 
  ylab("-log(p)") +
  ggrepel::geom_text_repel(data=head(output_tukey),aes(label=term),size=3, show.legend = F, max.overlaps = 100) +
  theme(legend.title = element_blank())            

plot_tukey

svglite("TukeyHSD_N_D.svg", width=10, height=10, bg='white')
plot_tukey
dev.off()



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
head(output_tukey,10)

print(paste("Total no.of features on which tukey test was performed:",nrow(output_tukey)))
print(paste("No.of features that showed significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Significant")))
print(paste("No.of features that did not show significant difference between",output_tukey$contrast[1],":",sum(output_tukey$significant=="Nonsignificant")))

plot_tukey <- ggplot(output_tukey,aes(x=estimate,y=-log(adj.p.value,base=10),color=significant)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_jama() + 
  ylab("-log(p)") +
  ggrepel::geom_text_repel(data=head(output_tukey),aes(label=term),size=3, show.legend = F, max.overlaps = 100) +
  theme(legend.title = element_blank())            

plot_tukey

###################T-tests#####################
data_merge2 <- data_merge1[order(metadata$ATTRIBUTE_Diet),]
data_merge_nd1<- data_merge2[51:149,] 
data_merge_nd<- data_merge_nd1[,-1] 
t.test(data_merge_nd[,1]~data_merge_nd1$ATTRIBUTE_Diet) #showing t-test result for 1st feature

rain_out <- data.frame()
for (i in 1:ncol(data_merge_nd)){
  model <- t.test(data_merge_nd[,i] ~ data_merge_nd1$ATTRIBUTE_Diet) #perform t-test for each feature against rainfall
  t_tidy <- broom::tidy(model)
  rain_out <- bind_rows(rain_out, t_tidy) #combine outputs from all features
}

rownames(rain_out) <- colnames(data_merge_nd) # naming each row with feature name

rain_out <- arrange(rain_out, p.value) # arranging rain_out by p values
rain_out["p_bonferroni"] <- p.adjust(rain_out$p.value,method="bonferroni")
rain_out["significant"] <- ifelse(rain_out$p_bonferroni<0.05,"Significant","Nonsignificant")

head(rain_out,3)

volcano_rainfall <- ggplot(rain_out,aes(x=estimate,y=-log(p.value,base=10),color=significant)) +
  geom_point() +
  theme_classic() +
  scale_color_jama() +
  geom_vline(xintercept = 0, linetype="dashed") +
  xlab("Beta") +
  ylab("-log(p)") +
  theme(legend.title = element_blank()) +
  ggrepel::geom_text_repel(data=head(rain_out), aes(label=rownames(rain_out)[1:6]), size=3, show.legend = F, max.overlaps = 100)

volcano_rainfall
ggsave("volcano_plot.pdf",volcano_rainfall)



###################Kruskal-Wallis#####################
broom::tidy(kruskal.test(data_merge[,1] ~ data_merge1$'ATTRIBUTE_Diet'))
kruskall_out <- data.frame()
for (i in 1:ncol(data_merge)){
  kw_model <- broom::tidy(kruskal.test(data_merge[,i] ~ data_merge1$'ATTRIBUTE_Diet')) #perform Kruskall Wallis for each feature against sample area
  kruskall_out <- bind_rows(kruskall_out, kw_model) #combine outputs from all features
}

rownames(kruskall_out) <- colnames(data_merge) # naming each row with feature name

kruskall_out <- arrange(kruskall_out, p.value) # arranging rain_out by p values
kruskall_out["p_bonferroni"] <- p.adjust(kruskall_out$p.value,method="bonferroni")
kruskall_out["significant"] <- ifelse(kruskall_out$p_bonferroni<0.05,"Significant","Nonsignificant")

head(kruskall_out,20)
print(paste("Total no.of features on which Kruskal test was performed:",nrow(kruskall_out)))
print(paste("No.of features that showed significant difference:",sum(kruskall_out$significant=="Significant")))
print(paste("No.of features that did not show significant difference:",sum(kruskall_out$significant=="Nonsignificant")))

kw_sig_names <- kruskall_out%>% filter(significant=="Significant") %>% rownames(.) #getting the rownames of significant features from KW output
length(kw_sig_names)


#combining anova and kruskal-wallis outputs:
AOV_KW <- kruskall_out #taking krusall-wallis output in a new variable "AOV_KW"
colnames(AOV_KW) <- paste0('KW_',colnames(AOV_KW)) #renaming the column names with a prefix 'KW'
AOV_KW <- merge(anova_out,AOV_KW, by.x="term", by.y=0) # merging anova output with AOV_KW 
AOV_KW <- arrange(AOV_KW,p.value) # arranging by p value


head(AOV_KW,10)
dim(AOV_KW)

ggplot(AOV_KW,aes(x=statistic, y=KW_statistic)) + 
  labs(x="ANOVA F-statistic", y="Kruskal-Wallis H-statistic", title="Correlation between ANOVA and Kruskal-Wallis statistics") +
  geom_point() +
  ggrepel::geom_text_repel(data=head(AOV_KW), aes(label=term), size=3, show.legend = F, max.overlaps = 100)


cor(AOV_KW$statistic, AOV_KW$KW_statistic, method="spearman")

###################Dunn's post hoc test#####################
dunnTest_model <- dunnTest(data_merge[,kw_sig_names[1]], 
                           data_merge1$`ATTRIBUTE_Diet`, 
                           method="bonferroni")
dunnTest_model$res

dunn_out <- NULL
for (i in kw_sig_names){
  dunn_out[[i]] <- dunnTest(data_merge[,i] ~  data_merge1$`ATTRIBUTE_Diet`,method="bonferroni") #perform dunn-test on Kruskal test's significant features against sample area
}
names(dunn_out) <- kw_sig_names
length(dunn_out)

#combining all dunn test results of 'La_Jolla Reefs - Mission_Bay':

dunn_ML <- data.frame() #empty dataframe for dunn test between Mission-Bay and La Jolla Reefs

for (i in 1:length(dunn_out)){
  result <- filter(dunn_out[[i]]$res, Comparison == 'Deficient - Normal') # subsetting that particular comparison
  rownames(result) <- kw_sig_names[i] # naming each row with feature name
  dunn_ML <- bind_rows(dunn_ML, result) #combine outputs from all features
}


dunn_ML <- arrange(dunn_ML, P.adj) # arranging dunn_out by p values
dunn_ML["p_bonferroni"] <- p.adjust(dunn_ML$P.adj, method="bonferroni") #p adj for bonferroni
dunn_ML["significant"] <- ifelse(dunn_ML$p_bonferroni<0.05,"Significant","Nonsignificant")

head(dunn_ML)

print(paste("Total no.of features on which Dunn test was performed:",length(dunn_out)))
print(paste("No.of features that showed significant difference between",dunn_ML$Comparison[1],":",sum(dunn_ML$significant=="Significant")))
print(paste("No.of features that did not show significant difference between",dunn_ML$Comparison[1],":",sum(dunn_ML$significant=="Nonsignificant")))



plot_dunn <- ggplot(dunn_ML, aes(x=Z, y=-log(P.adj, base=10), color=significant)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_jama()+ ylab("-log(p)") +
  ggrepel::geom_text_repel(data = head(dunn_ML), #labelling the top 6 features
                           aes(label = rownames(head(dunn_ML))), size = 3, show.legend = FALSE, max.overlaps = 100) +
  theme(legend.title = element_blank())

plot_dunn

ggsave("volcano_plot_dunnTest.pdf",plot_dunn)



##PLOT ANOVA RESULTS
data_merge_lp<- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet", 'ATTRIBUTE_Study_Day') %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1") 

lp_new <- data_merge_lp %>% 
  select(matches("11763_|2779_|1447_|2048_|6550_|11993_|5680_|4131_|5479_|1995_|6845_|7655_|
                 3048_|1336_|3332_|4096_|5430_|6896_|8006_|8625_|ATTRIBUTE_Study_Day|ATTRIBUTE_Diet"))  #select all columns that start with 11763 or 45678


lp_new %>%
  select(ATTRIBUTE_Study_Day, 3:7, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

lp_new %>%
  select(ATTRIBUTE_Study_Day, 8:12, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

lp_new %>%
  select(ATTRIBUTE_Study_Day, 13:17, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

lp_new %>%
  select(ATTRIBUTE_Study_Day, 18:22, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

##PLOT KW RESULTS
data_merge_kw<- metadata %>% dplyr::select("...1", "ATTRIBUTE_Diet", 'ATTRIBUTE_Study_Day') %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("...1") 

kw_new <- data_merge_kw %>% 
  select(matches("7979_|7908_|7116_|7618_|6809_|2048_|7609_|8006_|7620_|6285_|7327_|1555_|8406_|1568_|5678_|1666_|6556_|6056_|12541_|7606_|ATTRIBUTE_Study_Day|ATTRIBUTE_Diet"))  #select all columns that start with 11763 or 45678

kw_new %>%
  select(ATTRIBUTE_Study_Day, 3:7, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

kw_new %>%
  select(ATTRIBUTE_Study_Day, 8:12, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

kw_new %>%
  select(ATTRIBUTE_Study_Day, 13:17, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

kw_new %>%
  select(ATTRIBUTE_Study_Day, 18:22, ATTRIBUTE_Diet) %>%
  pivot_longer(cols = -c(ATTRIBUTE_Study_Day, ATTRIBUTE_Diet), names_to = "Variable") %>%
  ggplot(aes(x = ATTRIBUTE_Study_Day, y = value, color = ATTRIBUTE_Diet)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  labs(x = "Study Day", y = "Data", color = "Diet")

