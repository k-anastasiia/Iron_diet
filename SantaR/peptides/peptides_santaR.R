if(!require("devtools")) install.packages("devtools")
devtools::install_github("adwolfer/santaR", ref="master")
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

setwd('~/Desktop/GitHub/SantaR/peptides')
getwd()

metadata<- read_csv("md_new.csv")
metadata <- metadata[-1:-10,-1] 
metadata <- metadata[,-34:-37]
metadata <- metadata[,-2:-30]
Peptides <- read_csv("Peptides_no_blank.csv")
Normalised <- read_csv("Normalised_Quant_table.csv")
colnames(Normalised)[1] ="row_ID"
Normalised[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(Normalised$'row_ID', '_', 3)

Library_Hits_Refined_AK <- read_delim("Peptides_no_blank.csv", 
                                      delim = ",", escape_double = FALSE, 
                                      trim_ws = TRUE)
Merged_Hits_peptides <- merge(Normalised, Library_Hits_Refined_AK, by.x = "row_ID", by.y = "Row_ID",  
                              all.x = FALSE, all.y = FALSE)
Merged_Hits_peptides_cut <- Merged_Hits_peptides [,-194:-201]
Merged_Hits_peptides_cut <- Merged_Hits_peptides_cut [,-189:-192]

Merged_Hits_peptides_cut_f <- Merged_Hits_peptides_cut  %>% column_to_rownames ('Compound_Name') %>% t() %>% as.data.frame %>% rownames_to_column
Merged_Hits_peptides_cut_f <- Merged_Hits_peptides_cut_f[-1,]
data_merge_santaR_peptides <- merge(Merged_Hits_peptides_cut_f, metadata, by.x = "rowname", by.y = "filename",  
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



Ile_Leu_Lys2  <-santaR_plot(SANTAObj_peptides[[2]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Leu-Lys*",xlab = "Days")
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



Ile_Leu_Lys  <-santaR_plot(SANTAObj_peptides[[4]], colorVect=c("royalblue1","red3","green4"),  title = "Ile-Leu-Lys",xlab = "Days")
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

#Main figure
metadata_ND <- metadata[order(metadata$ATTRIBUTE_Diet), ]
metadata_ND <- metadata_ND[-1:-62,]
data_merge_santaR_peptides_ND <- merge(Merged_Hits_peptides_cut_f, metadata_ND, by.x = "rowname", by.y = "filename",  
                                    all.x = FALSE, all.y = FALSE)
inputData_peptides_ND     <- data.frame(data_merge_santaR_peptides_ND[,2:22])
ind_peptides_ND           <- data_merge_santaR_peptides_ND[["ATTIBUTE_mouse_no"]]
time_peptides_ND          <- as.numeric(data_merge_santaR_peptides_ND[["ATTRIBUTE_Study_Day"]])
group_peptides_ND         <- data_merge_santaR_peptides_ND[["ATTRIBUTE_Diet"]]
SANTAObj_peptides_ND  <- santaR_auto_fit(inputData_peptides_ND, ind_peptides_ND, time_peptides_ND, group_peptides_ND, df=5)


data_merge_santaR_peptides_ND$ATTRIBUTE_Study_Day<-as.numeric(data_merge_santaR_peptides_ND$ATTRIBUTE_Study_Day)

Leu_EK <-santaR_plot(SANTAObj_peptides_ND[[19]], colorVect=c("red3","green4"),  title = "Leu-EK",xlab = "Days")
Leu_EK  <- Leu_EK  + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Leu_EK <- Leu_EK + theme(panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank()) 
Leu_EK
SANTAObj_peptides_ND[[19]]$general$pval.dist

Ile_Gly_Ile  <-santaR_plot(SANTAObj_peptides_ND[[13]], colorVect=c("red3","green4"),  title = "Ile-Gly-Ile ",xlab = "Days")
Ile_Gly_Ile   <- Ile_Gly_Ile     + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ile_Gly_Ile   <- Ile_Gly_Ile   + theme(panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       panel.grid.major.y = element_blank(),
                                       panel.grid.minor.y = element_blank()) 
Ile_Gly_Ile 
SANTAObj_peptides_ND[[13]]$general$pval.dist

Ser_Lys  <-santaR_plot(SANTAObj_peptides_ND[[6]], colorVect=c("red3","green4"),  title = "Ser-Lys",xlab = "Days")
Ser_Lys  <- Ser_Lys    + scale_x_continuous(breaks = data_merge_santaR_peptides$ATTRIBUTE_Study_Day)
Ser_Lys  <- Ser_Lys   + theme(panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.y = element_blank()) 
Ser_Lys 
SANTAObj_peptides_ND[[6]]$general$pval.dist


grid.arrange(Ser_Lys, Ile_Gly_Ile ,Leu_EK,ncol = 3)



