#https://github.com/simonezuffa/Casein_Manuscript/blob/main/Casein_Analysis.Rmd#L1031
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
library(mixOmics)

setwd('/home/anastasiia/Desktop/mouse_studies/Spearman_cor/')
getwd() 
EC <- read_csv("~/Desktop/mouse_studies/Spearman_cor/EC_counts_mouse_n.csv")
EC<- EC[order(EC$mouse_number), ]
Normalised_day96 <- read_csv("Normalised_day96.csv")
Normalised <- Normalised_day96 %>% column_to_rownames ('ID') %>% t() %>% as.data.frame %>% rownames_to_column("mouse_number")
Diet <- c(rep("Normal", 5), rep("Deficient", 5), rep("Overload", 5))
EC <- cbind(EC, Diet)
Normalised <- cbind(Normalised, Diet)
EC$'mouse_number' <- gsub("", "EC-", EC$`mouse_number`)
#merge<- merge(EC, Normalised, by = "mouse_number")
#data1 <- merge %>% column_to_rownames ('mouse_number') %>% t() %>% as.data.frame %>% rownames_to_column("SampleID")
#new_row <- c("Diet", rep("Normal", 5), rep("Deficient", 5), rep("Overload", 5))
#data1 <- rbind(data1, new_row)



#Normalised <- log10(Normalised)  # Apply log10 to the entire data frame
#EC[, c("Diet", "mouse_number")] <- NULL  # Remove the specified columns
#EC <- log10(EC)  # Apply log10 to the entire data frame
# Log transformation
#EC <- log10(select(EC, -c("Diet", "mouse_number")))
#FrontalCortex <- log10(select(A1A2IntegralsFC_filter, -c("FC_Sample ID", "FC_Group")))



Normalised1<-Normalised[,-5732]
Normalised1<-Normalised1[,-1]

EC1<-EC[,-1025]
EC1<-EC1[,-1]


# PLS between the two datasets. Extract correlation structure between components for DIABLO design matrix 
A1A2pls <- pls(Normalised1, EC1, ncomp = 2, scale = TRUE)
plotIndiv(A1A2pls, comp = 1:2, rep.space= 'XY-variate', group = Normalised$Diet, ind.names = FALSE,
          legend = TRUE, title = "PLS comp 1 - 2, XY-space", pch = 20, style = "lattice", centroid = TRUE, ellipse = TRUE)
cor(A1A2pls$variates$X, A1A2pls$variates$Y) %>% diag() # correlation close to 0.9

# Prepare datasets
data <- list(Normalised1 = as.matrix(Normalised1), 
             `EC1` = as.matrix(EC1))
Y <- Normalised$Diet
design <- matrix(0.8, ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
diag(design) = 0



# Preliminary DIABLO 
DIABLO <- block.splsda(X = data, Y = Y, ncomp = 4, design = design)
perf_DIABLO <- perf(DIABLO, validation = 'loo')
plot(perf_DIABLO)  # 3 components and centroid distance

# Tune number of variables to retain per component
test_keepX <- list(Urine = c(4:9, seq(10, 18, 2), seq(20,30,5)),
                   `Frontal Cortex` = c(4:9, seq(10, 18, 2), seq(20,40,5)))

tune_DIABLO <- tune.block.splsda(X = data, Y = Y, ncomp = 3, 
                                 test.keepX = test_keepX, design = design,
                                 validation = 'loo', dist = "centroids.dist")

# Final DIABLO model
DIABLO <- block.splsda(X = data, Y = Y, ncomp = 3, keepX = tune_DIABLO$choice.keepX, design = design)
plotDiablo(DIABLO)



# Plot 
plotIndiv(DIABLO, ind.names = FALSE, legend = TRUE, title = 'DIABLO', col.per.group = c("#3A383F", "#85BEDC", "#A6B0BB"),
          style = "lattice", ellipse = TRUE, centroid = TRUE, pch = 20, X.label = "Component 1", Y.label = "Component 2")

# Circos Plot - Evaluate correlation between variables in the two omics datasets
circosPlot(DIABLO, cutoff = 0.7, line = FALSE, size.labels = 0.5, comp = 1,
           color.blocks = c("#006a4e", "grey51"), showIntraLinks = FALSE, legend = TRUE,
           size.variables = 0.5, size.legend = 0.5)







































# extract training data
data = list(mRNA = data1[2:1024,], 
            miRNA = data1[1024:6753,])

# check dimension
lapply(data, dim)

# outcome
Y = data1$Diet
summary(Y)


###################################################
### code chunk number 3: DIABLO-analysis.Rnw:175-179
###################################################
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0
design 


###################################################
### code chunk number 4: DIABLO-analysis.Rnw:186-198
###################################################
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                           design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
t1 = proc.time()
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
t2 = proc.time()
running_time = t2 - t1; running_time

#perf.diablo  # lists the different outputs
plot(perf.diablo) 



###################################################
### code chunk number 5: DIABLO-analysis.Rnw:204-206
###################################################
perf.diablo$choice.ncomp$WeightedVote
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]


###################################################
### code chunk number 6: DIABLO-analysis.Rnw:216-233 (eval = FALSE)
###################################################
## set.seed(123) # for reproducibility
## test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
## 
## t1 = proc.time()
## tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
##                               test.keepX = test.keepX, design = design,
##                               validation = 'Mfold', folds = 10, nrepeat = 1,
##                               cpus = 2, dist = "centroids.dist")
## t2 = proc.time()
## running_time = t2 - t1; running_time
## 
## list.keepX = tune.TCGA$choice.keepX
## list.keepX
## 
## #save(tune.TCGA,list.keepX, running_time, file = 'RData/result-TCGA-diablo_design0.1.RData')


###################################################
### code chunk number 7: DIABLO-analysis.Rnw:238-253 (eval = FALSE)
###################################################
## #set.seed(123) # for reproducibility, only when the `cpus' argument is not used
## test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
## 
## t1 = proc.time()
## tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
##                               test.keepX = test.keepX, design = design, 
##                               validation = 'Mfold', folds = 10, nrepeat = 1, 
##                               dist = "centroids.dist", cpus = 2)
## t2 = proc.time()
## running_time = t2 - t1; running_time
## 
## list.keepX = tune.TCGA$choice.keepX
## list.keepX


###################################################
### code chunk number 8: DIABLO-analysis.Rnw:257-259
###################################################
load('RData/result-TCGA-diablo_design0.1.RData')
running_time


###################################################
### code chunk number 9: DIABLO-analysis.Rnw:264-266
###################################################
list.keepX = tune.TCGA$choice.keepX
tune.TCGA$choice.keepX


###################################################
### code chunk number 10: DIABLO-analysis.Rnw:272-275
###################################################
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)
#sgccda.res   # lists the different functions of interest related to that object


###################################################
### code chunk number 11: DIABLO-analysis.Rnw:280-281
###################################################
sgccda.res$design


###################################################
### code chunk number 12: DIABLO-analysis.Rnw:286-288
###################################################
# mRNA variables selected on component 1
selectVar(sgccda.res, block = 'mRNA', comp = 1)


###################################################
### code chunk number 13: DIABLO-analysis.Rnw:296-297
###################################################
plotDiablo(sgccda.res, ncomp = 1)


###################################################
### code chunk number 14: DIABLO-analysis.Rnw:305-306
###################################################
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE,  style="ggplot2")


###################################################
### code chunk number 15: DIABLO-analysis.Rnw:311-312
###################################################
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')


###################################################
### code chunk number 16: DIABLO-analysis.Rnw:321-323
###################################################
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2,2,2), col = c('darkorchid', 'brown1', 'lightgreen'))


###################################################
### code chunk number 17: DIABLO-analysis.Rnw:330-333
###################################################
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2,2,2), col = c('darkorchid', 'brown1', 'lightgreen'),
        comp.select = 1)


###################################################
### code chunk number 18: DIABLO-analysis.Rnw:340-343
###################################################
circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)


###################################################
### code chunk number 19: DIABLO-analysis.Rnw:352-354
###################################################
network(sgccda.res, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)


###################################################
### code chunk number 20: DIABLO-analysis.Rnw:359-364 (eval = FALSE)
###################################################
## # not run
## library(igraph)
## my.network = network(sgccda.res, blocks = c(1,2,3),
##         color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
## write.graph(my.network$gR, file = "myNetwork.gml", format = "gml")


###################################################
### code chunk number 21: DIABLO-analysis.Rnw:370-371
###################################################
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')


###################################################
### code chunk number 22: DIABLO-analysis.Rnw:376-378
###################################################
cimDiablo(sgccda.res, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8,20), legend.position = "right")


###################################################
### code chunk number 23: DIABLO-analysis.Rnw:387-402
###################################################
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
t1 = proc.time()
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10, 
                   dist = 'centroids.dist')
t2 = proc.time()
running_time = t2 - t1; running_time


#perf.diablo  # lists the different outputs

# Performance with Majority vote
perf.diablo$MajorityVote.error.rate

# Performance with Weighted vote
perf.diablo$WeightedVote.error.rate


###################################################
### code chunk number 24: DIABLO-analysis.Rnw:413-414
###################################################
auc.diablo = auroc(sgccda.res, roc.block = "miRNA", roc.comp = 2)


###################################################
### code chunk number 25: DIABLO-analysis.Rnw:421-428
###################################################
# prepare test set data: here one block (proteins) is missing
data.test.TCGA = list(mRNA = breast.TCGA$data.test$mrna, 
                      miRNA = breast.TCGA$data.test$mirna)

predict.diablo = predict(sgccda.res, newdata = data.test.TCGA)
# the warning message will inform us that one block is missing
#predict.diablo # list the different outputs


###################################################
### code chunk number 26: DIABLO-analysis.Rnw:433-437
###################################################
confusion.mat = get.confusion_matrix(truth = breast.TCGA$data.test$subtype, 
                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat
get.BER(confusion.mat)


###################################################
### code chunk number 27: DIABLO-analysis.Rnw:442-443
###################################################
sessionInfo()


###################################################
### code chunk number 28: DIABLO-analysis.Rnw:447-449
###################################################
# the end of that chapter. Command not to be run
Stangle('DIABLO-analysis.Rnw', encoding = 'utf8')


