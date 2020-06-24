library(yaml)
require(scales)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cluster)
library(ComplexHeatmap)
library(klaR)

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
####Testing by reducing the number of ENCODE states (HEt, ZNF, Tx and Quies into regulation that we cna't reach with our histones marks)


tss_intersect.dcast[tss_intersect.dcast == "Quies"] <- "Het"
tss_intersect.dcast[tss_intersect.dcast == "Tx"] <- "Het"
tss_intersect.dcast[tss_intersect.dcast == "ZNF_Rpts"] <- "Het"

## Removing some ENCODE samples
tss_intersect.dcast <- tss_intersect.dcast[,-c(1,2)]
head(tss_intersect.dcast)


  
####################################
## Clustering gene with gower distance
####################################
gower.dist.gene = daisy(tss_intersect.dcast, metric = c("gower"))


divisive.clust.gene <- diana(as.matrix(gower.dist.gene),diss = TRUE, keep.diss = FALSE)
  
write.table(as.matrix(divisive.clust.gene), file.path(rootfolder,radio$cluster.gower.diana), row.names=FALSE, col.names=FALSE,sep="\t")

#########################################
## Kmeans with Klar
#####################################

# http://www.cs.ust.hk/~qyang/Teaching/537/Papers/huang98extensions.pdf
# https://cran.r-project.org/web/packages/klaR/klaR.pdf
# https://dabblingwithdata.wordpress.com/2016/10/10/clustering-categorical-data-with-r/

#Using differents clustering
# 3, 5 , 10, 15, 20 

kmodes_3_all <- kmodes(tss_intersect.dcast[,2:6], 3, iter.max = 3, weighted = FALSE )

kmodes_5_all <- kmodes(tss_intersect.dcast[,2:6], 5, iter.max = 5, weighted = FALSE )

kmodes_10_all <- kmodes(tss_intersect.dcast[,2:6], 10, iter.max = 10, weighted = FALSE )

kmodes_15_all <- kmodes(tss_intersect.dcast[,2:6], 15, iter.max = 10, weighted = FALSE )

kmodes_20_all <- kmodes(tss_intersect.dcast[,2:6], 20, iter.max = 10, weighted = FALSE )

tss_intersect.dcast$kmeans3 = kmodes_3_all$cluster
tss_intersect.dcast$kmeans5 = kmodes_5_all$cluster
tss_intersect.dcast$kmeans10 = kmodes_10_all$cluster
tss_intersect.dcast$kmeans15 = kmodes_15_all$cluster
tss_intersect.dcast$kmeans20 = kmodes_20_all$cluster

tss_intersect.dcast <- tss_intersect.dcast[,c(1,2,4,5,6,3,7,8,9,10,11)]

write.table(tss_intersect.dcast, file.path(rootfolder,radio$kmeans.all), row.names=TRUE, col.names=TRUE,sep="\t")
