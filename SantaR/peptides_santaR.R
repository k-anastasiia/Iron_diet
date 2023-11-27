#install.packages("santaR")
# Install devtools
if(!require("devtools")) install.packages("devtools")
devtools::install_github("adwolfer/santaR", ref="master")
#p-value between groups based on distance between groupMeanCurves

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
library(RColorBrewer)
library(paletteer)
library(tidyverse)
library(gridExtra)

getwd() #to see the working directory 
#Setting a working directory in Google Colab:
setwd("/home/anastasiia/Desktop/mouse_studies/Heatmaps/Peptides/peptide_abundance")


scaled_table <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/peptide_abundance/CLR_Scaled.csv")
metadata<- read_csv("metadata_full (ak).csv")
Peptides <- read_csv("~/Desktop/mouse_studies/Heatmaps/Peptides/peptide_abundance/Peptides_no_blank.csv")
Normalised <- read_csv("Normalised_Quant_table.csv")



#NORMALISED TABLE
# Split name column into several names
colnames(Normalised)[1] ="row_ID"
Normalised[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(Normalised$'row_ID', '_', 3)

#Merge Tables
Library_Hits_Refined_AK <- read_delim("Peptides_no_blank.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Library_Hits_Refined_AK <- Library_Hits_Refined_AK [-26:-27,]
Library_Hits_Refined_AK <- Library_Hits_Refined_AK [-1,]
Merged_Hits_peptides <- merge(Normalised, Library_Hits_Refined_AK, by.x = "row_ID", by.y = "Row_ID",  
                              all.x = FALSE, all.y = FALSE)
Merged_Hits_peptides_cut <- Merged_Hits_peptides [,-193:-200]
Merged_Hits_peptides_cut <- Merged_Hits_peptides_cut [,-189:-191]

Merged_Hits_peptides_cut_f <- Merged_Hits_peptides_cut  %>% column_to_rownames ('Compound_Name') %>% t() %>% as.data.frame %>% rownames_to_column
Merged_Hits_peptides_cut_f <- Merged_Hits_peptides_cut_f[-1,]
metadata <- metadata[-188,-34:-37]
metadata <- metadata[,-2:-30]
metadata_ND <- metadata[order(metadata$ATTRIBUTE_Diet), ]
metadata_ND <- metadata_ND[-1:-62,]
data_merge_santaR_peptides <- merge(Merged_Hits_peptides_cut_f, metadata_ND, by.x = "rowname", by.y = "...1",  
                           all.x = FALSE, all.y = FALSE)

#all diets
inputData_peptides     <- data.frame(data_merge_santaR_peptides[,2:22])
ind_peptides           <- data_merge_santaR_peptides[["ATTIBUTE_mouse_no"]]
time_peptides          <- as.numeric(data_merge_santaR_peptides[["ATTRIBUTE_Study_Day"]])
group_peptides         <- data_merge_santaR_peptides[["ATTRIBUTE_Diet"]]
SANTAObj_peptides  <- santaR_auto_fit(inputData_peptides, ind_peptides, time_peptides, group_peptides, df=5)

data_merge_santaR_peptides$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_peptides$ATTRIBUTE_Study_Day)

santaR_plot(SANTAObj_peptides[[19]])



#
data_merge_santaR_peptides <- merge(Merged_Hits_peptides_cut_f, metadata, by.x = "rowname", by.y = "...1",  
                                    all.x = FALSE, all.y = FALSE)

#all diets
inputData_peptides     <- data.frame(data_merge_santaR_peptides[,2:22])
ind_peptides           <- data_merge_santaR_peptides[["ATTIBUTE_mouse_no"]]
time_peptides          <- as.numeric(data_merge_santaR_peptides[["ATTRIBUTE_Study_Day"]])
group_peptides         <- data_merge_santaR_peptides[["ATTRIBUTE_Diet"]]
SANTAObj_peptides  <- santaR_auto_fit(inputData_peptides, ind_peptides, time_peptides, group_peptides, df=5)

data_merge_santaR_peptides$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_peptides$ATTRIBUTE_Study_Day)


Ser_Lys  <-santaR_plot(SANTAObj_peptides[[6]], colorVect=c("royalblue1", "red3","green4"),  title = "Ser-Lys",xlab = "Days")
Ser_Lys  <- Ser_Lys    + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ser_Lys  <- Ser_Lys   + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank(),
                              axis.title.y = element_blank(),  
                              axis.text.y = element_blank()) 
Ser_Lys <-Ser_Lys +  theme(legend.position = "none")

Ile_Gly_Ile  <-santaR_plot(SANTAObj_peptides[[13]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Gly-Ile ",xlab = "Days")
Ile_Gly_Ile   <- Ile_Gly_Ile     + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Gly_Ile   <- Ile_Gly_Ile   + theme(panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank(),
                                       axis.title.y = element_blank(),  
                                       axis.text.y = element_blank())   
Ile_Gly_Ile <-Ile_Gly_Ile+ theme(legend.position = "none")


Leu_EK <-santaR_plot(SANTAObj_peptides[[19]], colorVect=c("royalblue1","red3","green4"),  title = "Leu-EK",xlab = "Days")
Leu_EK  <- Leu_EK  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_EK <- Leu_EK + theme(panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         axis.title.y = element_blank(),  
                         axis.text.y = element_blank(),
                  )   
Leu_EK<-Leu_EK+  theme(legend.position = "none")


Ile_Arg  <-santaR_plot(SANTAObj_peptides[[1]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Arg",xlab = "Days")
Ile_Arg <- Ile_Arg  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Arg<- Ile_Arg + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank(),
                          axis.title.y = element_blank(),  
                          axis.text.y = element_blank(),
                          axis.title.x = element_blank(), 
                          axis.text.x = element_blank())  
Ile_Arg
Ile_Arg<-Ile_Arg+  theme(legend.position = "none")





Ile_Leu_Lys2  <-santaR_plot(SANTAObj_peptides[[2]], colorVect=c("royalblue1","red3","green4"),  title = "Ile/Leu-Lys*",xlab = "Days")
Ile_Leu_Lys2  <- Ile_Leu_Lys2   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Leu_Lys2 <- Ile_Leu_Lys2  + theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    axis.title.y = element_blank(),  
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank(), 
                                    axis.text.x = element_blank())  
Ile_Leu_Lys2 
Ile_Leu_Lys2 <-Ile_Leu_Lys2 +  theme(legend.position = "none")



Ile_Leu_Lys  <-santaR_plot(SANTAObj_peptides[[4]], colorVect=c("royalblue1","red3","green4"),  title = "Ile/Leu-Lys",xlab = "Days")
Ile_Leu_Lys  <- Ile_Leu_Lys   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Leu_Lys <- Ile_Leu_Lys  + theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                                    axis.title.y = element_blank(),  
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank(), 
                                    axis.text.x = element_blank())   
Ile_Leu_Lys 

Ile_Leu_Lys<-Ile_Leu_Lys+  theme(legend.position = "none")


Ile_Lys  <-santaR_plot(SANTAObj_peptides[[5]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Lys",xlab = "Days")
Ile_Lys  <- Ile_Lys   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Lys <- Ile_Lys  + theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                            axis.title.y = element_blank(),  
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(), 
                            axis.text.x = element_blank())  
Ile_Lys

Ile_Lys <-Ile_Lys +  theme(legend.position = "none")


Ile_Arg2  <-santaR_plot(SANTAObj_peptides[[7]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Arg*",xlab = "Days")
Ile_Arg2  <- Ile_Arg2   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Arg2  <- Ile_Arg2   + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank(),
                              axis.title.y = element_blank(),  
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(), 
                              axis.text.x = element_blank())  
Ile_Arg2 
Ile_Arg2<-Ile_Arg2+  theme(legend.position = "none")




Leu_Leu_Arg  <-santaR_plot(SANTAObj_peptides[[8]], colorVect=c("royalblue1","red3","green4"),  title = "Leu-Leu-Arg ",xlab = "Days")
Leu_Leu_Arg   <- Leu_Leu_Arg   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_Leu_Arg  <- Leu_Leu_Arg  + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank(),
                              axis.title.y = element_blank(),  
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(), 
                              axis.text.x = element_blank())   
Leu_Leu_Arg 
Leu_Leu_Arg <-Leu_Leu_Arg +  theme(legend.position = "none")


Ile_Ile_Arg  <-santaR_plot(SANTAObj_peptides[[9]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Ile-Arg",xlab = "Days")
Ile_Ile_Arg  <- Ile_Ile_Arg + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Ile_Arg  <- Ile_Ile_Arg + theme(panel.grid.major.x = element_blank(),
                                     panel.grid.minor.x = element_blank(),
                                     panel.grid.major.y = element_blank(),
                                     panel.grid.minor.y = element_blank(),
                                    axis.title.y = element_blank(),  
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank(), 
                                    axis.text.x = element_blank())  
Ile_Ile_Arg
Ile_Ile_Arg<-Ile_Ile_Arg+  theme(legend.position = "none")



Ile_Tyr  <-santaR_plot(SANTAObj_peptides[[10]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Tyr ",xlab = "Days")
Ile_Tyr <- Ile_Tyr + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Tyr  <- Ile_Tyr + theme(panel.grid.major.x = element_blank(),
                                    panel.grid.minor.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    panel.grid.minor.y = element_blank(),
                            axis.title.y = element_blank(),  
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(), 
                            axis.text.x = element_blank())  
Ile_Tyr 
Ile_Tyr <-Ile_Tyr +  theme(legend.position = "none")



Leu_Glu_Leu  <-santaR_plot(SANTAObj_peptides[[11]], colorVect=c("royalblue1","red3","green4"),  title = "Leu-Glu-Leu ",xlab = "Days")
Leu_Glu_Leu <- Leu_Glu_Leu+ scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_Glu_Leu  <-Leu_Glu_Leu + theme(panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            axis.title.y = element_blank(),  
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(), 
                            axis.text.x = element_blank())   
Leu_Glu_Leu
Leu_Glu_Leu<-Leu_Glu_Leu+  theme(legend.position = "none")




Phe_Leu  <-santaR_plot(SANTAObj_peptides[[12]], colorVect=c("royalblue1","red3","green4"),  title = "Phe-Leu ",xlab = "Days")
Phe_Leu<- Phe_Leu+ scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Phe_Leu  <-Phe_Leu + theme(panel.grid.major.x = element_blank(),
                                   panel.grid.minor.x = element_blank(),
                                   panel.grid.major.y = element_blank(),
                                   panel.grid.minor.y = element_blank(),
                           axis.title.y = element_blank(),  
                           axis.text.y = element_blank(),
                           axis.title.x = element_blank(), 
                           axis.text.x = element_blank())   
Phe_Leu
Phe_Leu<-Phe_Leu+  theme(legend.position = "none")



Phe_Leu2 <-santaR_plot(SANTAObj_peptides[[14]], colorVect=c("royalblue1","red3","green4"),  title = "Phe-Leu*",xlab = "Days")
Phe_Leu2  <- Phe_Leu2  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Phe_Leu2  <- Phe_Leu2 + theme(panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank(),
                              axis.title.y = element_blank(),  
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(), 
                              axis.text.x = element_blank()) 
Phe_Leu2
Phe_Leu2<-Phe_Leu2+  theme(legend.position = "none")


Ile_Pro_Ile  <-santaR_plot(SANTAObj_peptides[[15]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Pro-Ile ",xlab = "Days")
Ile_Pro_Ile   <- Ile_Pro_Ile     + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Pro_Ile   <- Ile_Pro_Ile   + theme(panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank(),
                                       axis.title.y = element_blank(),  
                                       axis.text.y = element_blank(),
                                       axis.title.x = element_blank(), 
                                       axis.text.x = element_blank()) 
Ile_Pro_Ile 
Ile_Pro_Ile <-Ile_Pro_Ile +  theme(legend.position = "none")



Leu_Trp <-santaR_plot(SANTAObj_peptides[[16]], colorVect=c("royalblue1","red3","green4"),  title = "Leu-Trp ",xlab = "Days")
Leu_Trp  <- Leu_Trp   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_Trp  <- Leu_Trp  + theme(panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank(),
                             axis.title.y = element_blank(),  
                             axis.text.y = element_blank(),
                             axis.title.x = element_blank(), 
                             axis.text.x = element_blank())   
Leu_Trp
Leu_Trp<-Leu_Trp+  theme(legend.position = "none")



Leu_Trp2 <-santaR_plot(SANTAObj_peptides[[17]], colorVect=c("royalblue1","red3","green4"),  title = "Leu-Trp*",xlab = "Days")
Leu_Trp2  <- Leu_Trp2   + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_Trp2  <- Leu_Trp2  + theme(panel.grid.major.x = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             axis.title.y = element_blank(),  
                             axis.text.y = element_blank(),
                             axis.title.x = element_blank(), 
                             axis.text.x = element_blank()) 
Leu_Trp2
Leu_Trp2<-Leu_Trp2+  theme(legend.position = "none")


Phe_Phe <-santaR_plot(SANTAObj_peptides[[18]], colorVect=c("royalblue1","red3","green4"),  title = "Phe-Phe",xlab = "Days")
Phe_Phe  <- Phe_Phe  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Phe_Phe  <- Phe_Phe + theme(panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank(),
                            axis.title.y = element_blank(),  
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(), 
                            axis.text.x = element_blank())   
Phe_Phe
Phe_Phe<-Phe_Phe+  theme(legend.position = "none")


Ile_Lys <-santaR_plot(SANTAObj_peptides[[20]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Lys",xlab = "Days")
Ile_Lys  <- Ile_Lys  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Lys<-Ile_Lys+ theme(panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank(),
                        axis.title.y = element_blank(),  
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(), 
                        axis.text.x = element_blank())  
Ile_Lys
Ile_Lys<-Ile_Lys+  theme(legend.position = "none")


Leu_Glu_Lys <-santaR_plot(SANTAObj_peptides[[21]], colorVect=c("royalblue1","red3","green4"),  title = "Leu-Glu-Lys",xlab = "Days")
Leu_Glu_Lys <- Leu_Glu_Lys  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_Glu_Lys <-Leu_Glu_Lys + theme(panel.grid.major.x = element_blank(),  
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        axis.title.y = element_blank(),  
                        axis.text.y = element_blank(),
                        axis.title.x = element_blank(), 
                        axis.text.x = element_blank()) 
Leu_Glu_Lys
Leu_Glu_Lys<-Leu_Glu_Lys+ theme(legend.position = "none")


grid.arrange(Ser_Lys, Ile_Gly_Ile ,Leu_EK, Ile_Arg,Ile_Leu_Lys2 ,Ile_Leu_Lys ,Ile_Lys,Ile_Arg,Leu_Leu_Arg,Ile_Ile_Arg,
             Ile_Tyr ,Leu_Glu_Leu,Phe_Leu, Phe_Leu2,Ile_Pro_Ile ,Leu_Trp, Leu_Trp2,Phe_Phe,Ile_Lys,Leu_Glu_Lys, ncol = 5)

grid.arrange(Ser_Lys, Ile_Gly_Ile ,Leu_EK,ncol = 3)




###############
Leu_EK <-santaR_plot(SANTAObj_peptides[[19]], colorVect=c("red3","green4"),  title = "Leu-EK",xlab = "Days")
Leu_EK  <- Leu_EK  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_EK <- Leu_EK + theme(panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank()) 
Leu_EK
SANTAObj_peptides[[19]]$general$pval.dist

Ile_Gly_Ile  <-santaR_plot(SANTAObj_peptides[[13]], colorVect=c("red3","green4"),  title = "Ile-Gly-Ile ",xlab = "Days")
Ile_Gly_Ile   <- Ile_Gly_Ile     + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Gly_Ile   <- Ile_Gly_Ile   + theme(panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank()) 
Ile_Gly_Ile 
SANTAObj_peptides[[13]]$general$pval.dist

Ser_Lys  <-santaR_plot(SANTAObj_peptides[[6]], colorVect=c("red3","green4"),  title = "Ser-Lys",xlab = "Days")
Ser_Lys  <- Ser_Lys    + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ser_Lys  <- Ser_Lys   + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank()) 
Ser_Lys 
SANTAObj_peptides[[6]]$general$pval.dist










p1 <-santaR_plot(SANTAObj[[1]])
p2 <-santaR_plot(SANTAObj[[2]])
p3 <-santaR_plot(SANTAObj[[3]])
p4 <-santaR_plot(SANTAObj[[4]])
p5 <-santaR_plot(SANTAObj[[5]])
p6 <-santaR_plot(SANTAObj[[6]])
p7 <-santaR_plot(SANTAObj[[7]])
p8 <-santaR_plot(SANTAObj[[8]])
p9 <-santaR_plot(SANTAObj[[9]])
p10 <-santaR_plot(SANTAObj[[10]])
p11 <-santaR_plot(SANTAObj[[11]])
p12 <-santaR_plot(SANTAObj[[12]])
p13 <-santaR_plot(SANTAObj[[13]])
p14 <-santaR_plot(SANTAObj[[14]])
p15 <-santaR_plot(SANTAObj[[15]])
p16 <-santaR_plot(SANTAObj[[16]])
p17 <-santaR_plot(SANTAObj[[17]])
p18 <-santaR_plot(SANTAObj[[18]])
p19 <-santaR_plot(SANTAObj[[19]])
p20 <-santaR_plot(SANTAObj[[20]])
p21 <-santaR_plot(SANTAObj[[21]])
p6<-p19 + scale_x_continuous(breaks = data_merge_santaR$ATTRIBUTE_Study_Day)
p6<-p8 + scale_x_continuous(breaks = data_merge_santaR$ATTRIBUTE_Study_Day)

data_merge_santaR$ATTRIBUTE_Study_Day<-as.numeric (data_merge_santaR$ATTRIBUTE_Study_Day)
grid.arrange(p1, p2, p3, p4)
grid.arrange(p5, p6, p7, p8)
grid.arrange(p9, p10, p11, p12)
grid.arrange(p13, p14, p15, p16)
grid.arrange(p17, p18, p19, p20, p21)
#grid.arrange(p1, p2, ncol=2) # force them side by side
p17<-p17 + scale_x_continuous(breaks = data_merge_santaR$ATTRIBUTE_Study_Day)
p13<-p13 + scale_x_continuous(breaks = data_merge_santaR$ATTRIBUTE_Study_Day)
p1<-p1 + scale_x_continuous(breaks = data_merge_santaR$ATTRIBUTE_Study_Day)
grid.arrange(p17, p13, p1)


santaR_plot(SANTAObj[[1]], title = "Ile-Arg", legend = TRUE, showIndPoint = TRUE,
            showIndCurve = TRUE, showGroupMeanCurve = TRUE,
            showTotalMeanCurve = FALSE, showConfBand = TRUE, colorVect = NA,
            sampling = 250, xlab = "x", ylab = "y", shortInd = FALSE)


# Plot significant metabolites
SantaR_met <- ggarrange(Succinate_timeFluc_Group, Sarcosine_timeFluc_Group, Inosine_timeFluc_Sex,
                        legend = "right", align = "hv", ncol = 1, nrow = 3)














#########EXTRA
#exclude 0.6% Fe
data_merge_santaR<-data_merge_santaR[order(data_merge_santaR$ATTRIBUTE_Diet), ]
data_merge_santa_ND<-data_merge_santaR[-1:-62,]

inputData     <- data.frame(data_merge_santa_ND[,2:22])
ind           <- data_merge_santa_ND[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santa_ND[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santa_ND[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

#exclude mouse 6
data_6<-data_merge_santaR[order(data_merge_santaR$ATTIBUTE_mouse_no), ]
data_6<-data_6[-140:-148,]


inputData     <- data.frame(data_6[,2:22])
ind           <- data_6[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_6[["ATTRIBUTE_Study_Day"]])
group         <- data_6[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)

#CLR Scaled table
Library_Hits_Refined_AK <- read_delim("Peptides_no_blank.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Library_Hits_Refined_AK <- Library_Hits_Refined_AK [-26:-27,]
Library_Hits_Refined_AK <- Library_Hits_Refined_AK [-1,]
View(Library_Hits_Refined_AK)

Merged_Hits_peptides <- merge(scaled_table_t, Library_Hits_Refined_AK, by.x = "row_ID", by.y = "Row_ID",  
                              all.x = FALSE, all.y = FALSE)
Merged_Hits_peptides_cut <- Merged_Hits_peptides [,-193:-200]
Merged_Hits_peptides_cut <- Merged_Hits_peptides_cut [,-189:-191]

Merged_Hits_peptides_cut_f <- Merged_Hits_peptides_cut  %>% column_to_rownames ('Compound_Name') %>% t() %>% as.data.frame %>% rownames_to_column
Merged_Hits_peptides_cut_f <- Merged_Hits_peptides_cut_f[-1,]
metadata <- metadata[,-34:-37]
metadata <- metadata[,-2:-30]
data_merge_santaR <- merge(Merged_Hits_peptides_cut_f, metadata, by.x = "rowname", by.y = "...1",  
                              all.x = FALSE, all.y = FALSE)

#change the column name
colnames(scaled_table)[1] <- "row_ID"
scaled_table_t <- scaled_table %>% column_to_rownames ('row_ID') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")

# Split name column into several names
colnames(scaled_table_t)[1] ="row_ID"
scaled_table_t[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(scaled_table_t$'row_ID', '_', 3)
View(scaled_table_t)
 
data_merge_santaR<-data_merge_santaR[order(data_merge_santaR$ATTRIBUTE_Diet), ]
data_merge_santa_ND<-data_merge_santaR[-1:-50,]



inputData     <- data.frame(data_merge_santa_ND[,2:22])
ind           <- data_merge_santa_ND[["ATTIBUTE_mouse_no"]]
time          <- as.numeric(data_merge_santa_ND[["ATTRIBUTE_Study_Day"]])
group         <- data_merge_santa_ND[["ATTRIBUTE_Diet"]]
SANTAObj  <- santaR_auto_fit(inputData, ind, time, group, df=5)


santaR_plot(SANTAObj[[1]], title = "", legend = TRUE, showIndPoint = TRUE,
            showIndCurve = TRUE, showGroupMeanCurve = TRUE,
            showTotalMeanCurve = FALSE, showConfBand = TRUE, colorVect = NA,
            sampling = 250, xlab = "x", ylab = "y", shortInd = FALSE)


santaR_auto_fit(inputData, ind, time, group, df=5, ncores = 0,
                CBand = TRUE, pval.dist = TRUE, pval.fit = FALSE, nBoot = 1000,
                alpha = 0.05, nPerm = 1000, nStep = 5000, alphaPval = 0.05,
                forceParIndTimeMat = FALSE)

santaR_plot(SANTAObj, title = "", legend = TRUE, showIndPoint = TRUE,
            showIndCurve = TRUE, showGroupMeanCurve = TRUE,
            showTotalMeanCurve = FALSE, showConfBand = TRUE, colorVect = NA,
            sampling = 250, xlab = "x", ylab = "y", shortInd = FALSE)
#santaR_start_GUI(browser = TRUE)

SRs <- santaR_auto_fit(inputData = data_merge_santa_ND[,2:5], ind = data_merge_santa_ND[["ATTIBUTE_mouse_no"]],
                       time = as.numeric(data_merge_santa_ND[["ATTRIBUTE_Study_Day"]]), group = data_merge_santa_ND[["ATTRIBUTE_Diet"]],
                       df = 5)# Generate plots

santaR_plot(SRs, title = "", legend = TRUE, showIndPoint = TRUE,
            showIndCurve = TRUE, showGroupMeanCurve = TRUE,
            showTotalMeanCurve = FALSE, showConfBand = TRUE, colorVect = NA,
            sampling = 250, xlab = "x", ylab = "y", shortInd = FALSE)
############EXTRA END




# Simone's Code
SRg <- santaR_auto_fit(inputData = metlogINFO[,met_names], ind = metlogINFO[["Subject"]],
                       time = as.numeric(metlogINFO[["Timepoint"]]), group = metlogINFO[["Group"]],
                       df = 5)
SRs <- santaR_auto_fit(inputData = metlogINFO[,met_names], ind = metlogINFO[["Subject"]],
                       time = as.numeric(metlogINFO[["Timepoint"]]), group = metlogINFO[["Sex"]],
                       df = 5)# Generate plots
for (n in 1:2) {
  for (k in 1:length(met_names)) {
    if (n == 1) {
      x <- SRg
      ob <- "Group"
      pal <- PalGroup
    } else {
      x <- SRs
      ob <- "Sex"
      pal <- PalSex
    }    p <- santaR_plot(x[[met_names[k]]], showIndPoint = FALSE, showIndCurve = FALSE,
                          xlab = NULL, ylab = NULL, colorVect = pal,
                          title = paste(met_names[k],"-", ob, "Fluctuation - p value =",
                                        round(x[[met_names[k]]]$general$pval.dist, digits = 3))) +
      scale_x_continuous(breaks = 1:length(Sampling), labels = Sampling) +
      theme(plot.title = element_text(size = 8),
            axis.text = element_text(size = 5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())  assign(paste(met_names[k], "_timeFluc_", ob, sep = ""), p)
  }
}# Plot significant metabolites
SantaR_met <- ggarrange(Succinate_timeFluc_Group, Sarcosine_timeFluc_Group, Inosine_timeFluc_Sex,
                        legend = "right", align = "hv", ncol = 1, nrow = 3)
      
