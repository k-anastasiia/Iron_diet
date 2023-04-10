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

setwd("/home/anastasiia/Desktop/mouse_studies/Heme_breakdown")
getwd() 

metadata <- read_csv("~/Desktop/mouse_studies/Heme_breakdown/metadata_full (ak).csv")
scaled_table <- read_csv("~/Desktop/mouse_studies/Heme_breakdown/CLR_Scaled.csv")
master_table <- read_csv("~/Desktop/mouse_studies/Heme_breakdown/Library_Hits_Heme.csv")

#change the column name
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")

# Split name column into several names
colnames(scaled_table_t)[1] ="row_ID"
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'row_ID', '_', 3)

#Merge Feature and Library Hits Tables
Merged_Hits<- merge(scaled_table_t, master_table, by.x = "row_ID", by.y = "Row_ID",  
                              all.x = FALSE, all.y = FALSE)
Merged_Hits_cut <- Merged_Hits [,-193:-253]
Merged_Hits_cut2 <- Merged_Hits_cut [,-189:-191]
Merged_Hits_cut2[5,] <-colnames(Merged_Hits_cut2) 
Merged_Hits_cut2_t <- Merged_Hits_cut2 %>% column_to_rownames ('Compound_Name')   %>% t() %>% as.data.frame 
colnames(Merged_Hits_cut2_t)[5] <- "Row_ID"

#Merge with metadata
colnames(metadata)[1] <- "row_ID"
metadata_t <- metadata %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")
Merged_Hits_Heme<- merge(Merged_Hits_cut2_t, metadata, by.x = "Row_ID", by.y = "row_ID",  
                    all.x = TRUE, all.y = TRUE)
Merged_Hits_Heme <-  Merged_Hits_Heme %>% column_to_rownames ('Row_ID')  
Merged_Hits_Heme_cut <- Merged_Hits_Heme [,-37:-40]
Merged_Hits_Heme_cut2 <- Merged_Hits_Heme_cut [-187:-189,-5:-33]


#BILIRUBIN
 dat <- Merged_Hits_Heme_cut2[order(Merged_Hits_Heme_cut2$ATTRIBUTE_Diet),] 
 group <- c("0.6% Fe","Deficient", "Normal" )
 dat$Diet <- ifelse(dat$ATTRIBUTE_Diet == "0.6% Fe", group[1], ifelse(dat$ATTRIBUTE_Diet == "Deficient", group[2], group[3]))
 dat$Study_Day <- as.numeric(as.character(dat$ATTRIBUTE_Study_Day))

 bilirubin<- ggplot(dat, aes(x = Study_Day, y = as.numeric(Bilirubin), colour = Diet)) +
   geom_point(size = 2) +
   scale_color_manual(values = c('royalblue1','red3','green4')) + 
   scale_y_continuous(limits = c(-2.2, 4.8))+
   labs(x = "Study Day", y = "Bilirubin") +
   theme_minimal() +
   scale_x_continuous(breaks = seq(0, 96, by = 8))+
   stat_summary(fun = "mean", geom = "line", size = 1, aes(group = Diet))+
   geom_text(aes(label = ATTIBUTE_mouse_no), vjust = -1, hjust = 1.2, check_overlap = TRUE)
   #theme(axis.line = element_line(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 
   ggsave("bilirubin.svg", plot = stercobilin, width = 8, height = 6, units = "in")
   
#BILIVERDIN
   dat <- Merged_Hits_Heme_cut2[order(Merged_Hits_Heme_cut2$ATTRIBUTE_Diet),] 
   group <- c("0.6% Fe","Deficient", "Normal" )
   dat$Diet <- ifelse(dat$ATTRIBUTE_Diet == "0.6% Fe", group[1], ifelse(dat$ATTRIBUTE_Diet == "Deficient", group[2], group[3]))
   dat$Study_Day <- as.numeric(as.character(dat$ATTRIBUTE_Study_Day))
   
   biliverdin<-ggplot(dat, aes(x = Study_Day, y = as.numeric(Biliverdin), colour = Diet)) +
     geom_point(size = 2) +
     scale_color_manual(values = c('royalblue1','red3','green4')) + 
     scale_y_continuous(limits = c(-2.2, 5.1))+
     labs(x = "Study Day", y = "Biliverdin") +
     theme_minimal() +
     scale_x_continuous(breaks = seq(0, 96, by = 8))+
     stat_summary(fun = "mean", geom = "line", size = 1, aes(group = Diet))+
     geom_text(aes(label = ATTIBUTE_mouse_no), vjust = -1, hjust = 1.2, check_overlap = TRUE)
   
   ggsave("biliverdin.svg", plot = stercobilin, width = 8, height = 6, units = "in")
    
#UROBILIN
   dat <- Merged_Hits_Heme_cut2[order(Merged_Hits_Heme_cut2$ATTRIBUTE_Diet),] 
   group <- c("0.6% Fe","Deficient", "Normal" )
   dat$Diet <- ifelse(dat$ATTRIBUTE_Diet == "0.6% Fe", group[1], ifelse(dat$ATTRIBUTE_Diet == "Deficient", group[2], group[3]))
   dat$Study_Day <- as.numeric(as.character(dat$ATTRIBUTE_Study_Day))
   
   urobilin<- ggplot(dat, aes(x = Study_Day, y = as.numeric(Urobilin), colour = Diet)) +
     geom_point(size = 2) +
     scale_color_manual(values = c('royalblue1','red3','green4')) + 
     scale_y_continuous(limits = c(-1.7, 7.8))+
     labs(x = "Study Day", y = "Urobilin") +
     theme_minimal() +
     scale_x_continuous(breaks = seq(0, 96, by = 8))+
     stat_summary(fun = "mean", geom = "line", size = 1, aes(group = Diet))+
     geom_text(aes(label = ATTIBUTE_mouse_no), vjust = -1, hjust = 1.2, check_overlap = TRUE)
    
     ggsave("urobilin.svg", plot = stercobilin, width = 8, height = 6, units = "in")
   
   
#STERCOBILIN
   dat <- Merged_Hits_Heme_cut2[order(Merged_Hits_Heme_cut2$ATTRIBUTE_Diet),] 
   group <- c("0.6% Fe","Deficient", "Normal" )
   dat$Diet <- ifelse(dat$ATTRIBUTE_Diet == "0.6% Fe", group[1], ifelse(dat$ATTRIBUTE_Diet == "Deficient", group[2], group[3]))
   dat$Study_Day <- as.numeric(as.character(dat$ATTRIBUTE_Study_Day))
   
   stercobilin<- ggplot(dat, aes(x = Study_Day, y = as.numeric(Stercobilin), colour = Diet)) +
     geom_point(size = 2) +
     scale_color_manual(values = c('royalblue1','red3','green4')) + 
     scale_y_continuous(limits = c(-1.9, 7.3))+
     labs(x = "Study Day", y = "Stercobilin") +
     theme_minimal() +
     scale_x_continuous(breaks = seq(0, 96, by = 8))+
     stat_summary(fun = "mean", geom = "line", size = 1, aes(group = Diet))+
     geom_text(aes(label = ATTIBUTE_mouse_no), vjust = -1, hjust = 1.2, check_overlap = TRUE)
  
      ggsave("stercobilin.svg", plot = stercobilin, width = 8, height = 6, units = "in")
      
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
################OLD#######################
   #Sort by diet
   Merged_Hits_Heme_all <- Merged_Hits_Heme_cut2[order(Merged_Hits_Heme_cut3$ATTRIBUTE_Diet),] 
   
   # Create plot with subsetted colors
   plot(x = Merged_Hits_Heme_all$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_all$Bilirubin, col = colors, pch = 16)
   
   #plot
   dat <-Merged_Hits_Heme_all
   plot(x = Merged_Hits_Heme_all$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_all$Bilirubin,   col = c("Normal"="green4","Deficient"="red3","0.6% Fe"="royalblue1"), pch = 16)
   plot(x = Merged_Hits_Heme_all$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_all$Biliverdin, col = c("Normal"="green4","Deficient"="red3","0.6% Fe"="royalblue1"),pch = 16)
   plot(x = Merged_Hits_Heme_all$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_all$Urobilin, col = c("Normal"="green4","Deficient"="red3","0.6% Fe"="royalblue1"), pch = 16)      
   plot(x = Merged_Hits_Heme_all$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_all$Stercobilin, col = c("Normal"="green4","Deficient"="red3","0.6% Fe"="royalblue1"), pch = 16)
   
#DEFICIENT
   #Sort by diet
   Merged_Hits_Heme_cut3 <- Merged_Hits_Heme_cut2[order(Merged_Hits_Heme_cut2$ATTRIBUTE_Diet),] 
   #Filter by diet
   Merged_Hits_Heme_D <- Merged_Hits_Heme_cut3 %>%  filter(Merged_Hits_Heme_cut3$ATTRIBUTE_Diet == "Deficient") 
   
   #plot
   dat <-Merged_Hits_Heme_D 
   plot(x = Merged_Hits_Heme_D$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_D$Bilirubin,  col = 'red3', pch = 16)
   plot(x = Merged_Hits_Heme_D$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_D$Biliverdin,  col = 'red3',pch = 16)
   plot(x = Merged_Hits_Heme_D$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_D$Urobilin,  col = 'red3', pch = 16)      
   plot(x = Merged_Hits_Heme_D$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_D$Stercobilin,  col = 'red3', pch = 16)
   
#0.6% Fe
   #Filter by diet
   Merged_Hits_Heme_E <- Merged_Hits_Heme_cut3 %>%  filter(Merged_Hits_Heme_cut3$ATTRIBUTE_Diet == "0.6% Fe") #Filter by diet
   
   #plot
   dat <-Merged_Hits_Heme_E 
   plot(x = Merged_Hits_Heme_E$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_E$Bilirubin,  col = 'royalblue1', pch = 16)
   plot(x = Merged_Hits_Heme_E$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_E$Biliverdin,  col = 'royalblue1', pch = 16)
   plot(x = Merged_Hits_Heme_E$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_E$Urobilin,  col ='royalblue1', pch = 16)      
   plot(x = Merged_Hits_Heme_E$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_E$Stercobilin,  col = 'royalblue1', pch = 16)
   
#NORMAL
   #Filter by diet
   Merged_Hits_Heme_N <- Merged_Hits_Heme_cut3 %>%  filter(Merged_Hits_Heme_cut3$ATTRIBUTE_Diet == "Normal") #Filter by diet
   
   #plot
   dat <-Merged_Hits_Heme_N 
   plot(x = Merged_Hits_Heme_N$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_N$Bilirubin,  col = 'green4', pch = 16)
   plot(x = Merged_Hits_Heme_N$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_N$Biliverdin,  col = 'green4', pch = 16)
   plot(x = Merged_Hits_Heme_N$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_N$Urobilin,  col = 'green4', pch = 16)      
   plot(x = Merged_Hits_Heme_N$ATTRIBUTE_Study_Day , y =  Merged_Hits_Heme_N$Stercobilin,  col = 'green4', pch = 16)
   
   