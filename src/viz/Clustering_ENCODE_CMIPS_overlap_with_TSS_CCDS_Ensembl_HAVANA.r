
library(yaml)
library(tidyr)
require(scales)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cluster)
library(ComplexHeatmap)

radio = read_yaml("/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/radiofile.yml")
#######################################################
####  Color scheme
#######################################################
colour15 = read.table(radio$color.mapping.encode.15,sep="\t")

colours18 = read.table(radio$color.mapping.encode.18,sep="\t")

colourcmips= read.table(radio$color.mapping.cmips,sep="\t")

all_colors = rbind(colour15,colours18,colourcmips)

getrgb <-function(x){
  rgb(matrix(as.vector(unlist(strsplit(as.character(x),","))),ncol=3),maxColorValue=255)
}

chromHMM_color_scheme = sapply(all_colors$V2,getrgb)
names(chromHMM_color_scheme) <- all_colors$V1
chromHMM_color_scheme


tss_intersect.dcast = read.table(radio$intersect.tss,h=F,sep = "\t")


#test = tss_intersect[1:1000,]
tss_intersect.dcast = dcast(tss_intersect.dcast,V1~V2)
rownames(tss_intersect.dcast) <- tss_intersect.dcast$V1
tss_intersect.dcast = tss_intersect.dcast[,2:length(tss_intersect.dcast[1,])]
tss_intersect.dcast[sapply(tss_intersect.dcast, is.character)] <- lapply(tss_intersect.dcast[sapply(tss_intersect.dcast, is.character)], 
                                                                         as.factor)

head(tss_intersect.dcast)
####Testing byt reducing the number of ENCODE states (HEt, ZNF, Tx and Quies into regulation)


tss_intersect.dcast[tss_intersect.dcast == "Quies"] <- "Het"
tss_intersect.dcast[tss_intersect.dcast == "Tx"] <- "Het"
tss_intersect.dcast[tss_intersect.dcast == "ZNF_Rpts"] <- "Het"

##REmoving some ENCODE samples
tss_intersect.dcast <- tss_intersect.dcast[,-2]
head(tss_intersect.dcast)

plotheatmap <- function(matrix,nb_clust_gene,filename){
  
  ####################################
  ## Clustering gene
  ####################################
  gower.dist.gene = daisy(matrix, metric = c("gower"))
  divisive.clust.gene <- diana(as.matrix(gower.dist.gene), 
                               diss = FALSE, keep.diss = FALSE)
  
  #plot(divisive.clust.gene, main = "Divisive")
  
  row_dend = as.dendrogram(divisive.clust.gene)
  row_dend.cut = cutree(divisive.clust.gene,k=nb_clust_gene)
  
  ####################################
  ## Clustering samples
  ####################################
  
  matrix.t=t(matrix)
  matrix.t=as.data.frame(matrix.t)
  
  gower.dist.tissue = daisy(matrix.t, metric = c("gower"))
  
  divisive.clust.tissue <- diana(as.matrix(gower.dist.tissue), 
                                 diss = TRUE, keep.diss = TRUE)
  col_dend = as.dendrogram(divisive.clust.tissue)
  #######################################################
  ####  Clustering
  #######################################################
  split = data.frame(row_dend.cut)
  
  pdf(paste0(filename,".pdf"),width=8, height=8)
  ComplexHeatmap::Heatmap(matrix,
                          col = chromHMM_color_scheme,
                          row_split = split,
                          cluster_columns=col_dend,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
  dev.off()
}

##Trying it with only Active promoter sites
tssA = tss_intersect.dcast[tss_intersect.dcast$cmips == "Active_promoter",]
plotheatmap(tssA,10,"tss_heatmap_k10_diana")

biv =  tss_intersect.dcast[tss_intersect.dcast$cmips == "Poised_promoter",]
plotheatmap(biv,10,"bivtss_heatmap_k10_diana")

ENCODE_samps = c("HRT.VENT.L","HRT.VENT.R","HRT.ATR.R","HRT.FET","BRN.ADLT","BRN.FET","IPS")  
ENCODE.col = "#e41a1c"
CM_IPS.col = "#377eb8"

fetal.col = "#c51b7d"
adult.col = "#4d9221"
ips.col = "#f1b6da"

right.atrium.col = "#66c2a5"
left.ventricule.col = "#8da0cb"
right.ventricule.col ="#a6d854"
brain.col = "#e5c494"
heart.col = "#e41a1c"
#Colors from 8-class Set2 R colorbrewer
source.cols = c(rep(ENCODE.col,length(ENCODE_samps)),rep(CM_IPS.col,length(our_samps)))
source.cols
stage.cols = c(adult.col,adult.col,adult.col,
               fetal.col,adult.col,
               fetal.col,
               ips.col,
               rep(CM_IPS.col,length(our_samps)))
stage.cols
tissue.cols = c(left.ventricule.col,right.ventricule.col,right.atrium.col,
                fetal.col,
                brain.col,brain.col,
                ips.col,
                CM_IPS.col,CM_IPS.col,CM_IPS.col,CM_IPS.col)
tissue.cols
rlab=t(cbind(source.cols,stage.cols,tissue.cols))

rownames(rlab)=c("Data origin","Cell type", "Tissue type")
# Heatmap(test.dcast, name = "mat",col = chromHMM_color_scheme,  cluster_rows=row_dend)
