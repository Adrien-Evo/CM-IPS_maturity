library(yaml)
require(scales)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cluster)
library(ComplexHeatmap)

rootfolder = "/ceph-recherche/shares/u1087/afoucal/CM-IPS_maturity/"
radio = read_yaml(file.path(rootfolder,"radiofile.yml"))

#' This script is used to generate clustering from s the intersected chromHMM states for ENCODE samples and our merged states. It only looks at the TSS start point.
#' Since the variables are categorical, it uses the gower distance 

tss_intersect.dcast = read.table(file.path(rootfolder,radio$intersect.tss),h=F,sep = "\t")


#test = tss_intersect[1:1000,]
tss_intersect.dcast = dcast(tss_intersect.dcast,V1~V2)
rownames(tss_intersect.dcast) <- tss_intersect.dcast$V1
tss_intersect.dcast = tss_intersect.dcast[,2:length(tss_intersect.dcast[1,])]
tss_intersect.dcast[sapply(tss_intersect.dcast, is.character)] <- lapply(tss_intersect.dcast[sapply(tss_intersect.dcast, is.character)], 
                                                                         as.factor)

head(tss_intersect.dcast)
####Testing byt reducing the number of ENCODE states (HEt, ZNF, Tx and Quies into regulation that we cna't reach with our histones marks)


tss_intersect.dcast[tss_intersect.dcast == "Quies"] <- "Het"
tss_intersect.dcast[tss_intersect.dcast == "Tx"] <- "Het"
tss_intersect.dcast[tss_intersect.dcast == "ZNF_Rpts"] <- "Het"

## Removing some ENCODE samples
tss_intersect.dcast <- tss_intersect.dcast[,-c(1,2)]
head(tss_intersect.dcast)


  
  ####################################
  ## Clustering gene
  ####################################
  gower.dist.gene = daisy(matrix, metric = c("gower"))
  write.table(as.matrix(gower.dist.gene), radio$cluster.full.gower, row.names=FALSE, col.names=FALSE,sep="\t")


divisive.clust.gene <- diana(as.matrix(gower.dist.gene), 
                               diss = TRUE, keep.diss = FALSE)
  