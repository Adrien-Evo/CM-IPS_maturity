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

tss_intersect = read.table(file.path(rootfolder,radio$atac.chrom),h=F,sep = "\t", stringsAsFactor = FALSE)


#test = tss_intersect[1:1000,]
atac_chromhmm = dcast(tss_intersect,V1~V2)
rownames(atac_chromhmm) <- atac_chromhmm$V1
atac_chromhmm = atac_chromhmm[,2:length(atac_chromhmm[1,])]
atac_chromhmm[sapply(atac_chromhmm, is.character)] <- lapply(atac_chromhmm[sapply(atac_chromhmm, is.character)], 
                                                                         as.factor)

head(atac_chromhmm)
####Testing by reducing the number of ENCODE states (HEt, ZNF, Tx and Quies into regulation that we cna't reach with our histones marks)

###I don't think that's necessary for clustering purposes
# 
atac_chromhmm[atac_chromhmm == "Quies"] <- "Het"
atac_chromhmm[atac_chromhmm == "Tx"] <- "Het"
atac_chromhmm[atac_chromhmm == "ZNF_Rpts"] <- "Het"
head(atac_chromhmm)

## Removing some ENCODE samples
atac_chromhmm <- atac_chromhmm[,-c(1,2)]
head(atac_chromhmm)
atac_chromhmm <- atac_chromhmm[,c(3,9,1,4,2,5,7,8,10,6)]
head(atac_chromhmm)

colnames(atac_chromhmm) <- c("ATAC.CMiPS","ATAC.iPS","ATAC.Brug","ATAC.Control","CM-iPS","HRT.ATR.R","HRT.VENT.L","HRT.VENT.R","IPS","HRT.FET")
head(atac_chromhmm)

kmeans = read.table(file.path(rootfolder,radio$kmeans.all), h=T, sep ="\t", stringsAsFactors=FALSE,check.names=FALSE)


indexmatch <- match(rownames(atac_chromhmm),rownames(kmeans))
###This methods works
##This is how I verified that the rownames do match propper. COuld maybe haved used merged using all = 0 but it takes some time.
#test = cbind(atac_chromhmm[!is.na(indexmatch),c(1,2)],kmeans[,], rownames(kmeans))
# all(rownames(test)== test[,15])

df = cbind(atac_chromhmm[!is.na(indexmatch),c(1,2)],kmeans[,])

write.table(df, file.path(rootfolder,radio$atac.chrom.kmeans), row.names=TRUE, col.names=TRUE,sep="\t")
