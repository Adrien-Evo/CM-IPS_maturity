library(yaml)
require(scales)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cluster)
library(ComplexHeatmap)



rootfolder = "/ceph-recherche/shares/u1087/afoucal/CM-IPS_maturity/"
radio = read_yaml(file.path(rootfolder,"radiofile.yml"))

#######################################################
####  Color scheme
#######################################################
colour15 = read.table(file.path(rootfolder,radio$color.mapping.encode.15),sep="\t")

colours18 = read.table(file.path(rootfolder,radio$color.mapping.encode.18),sep="\t")

colourcmips= read.table(file.path(rootfolder,radio$color.mapping.cmips),sep="\t")

all_colors = rbind(colour15,colours18,colourcmips)

getrgb <-function(x){
  rgb(matrix(as.vector(unlist(strsplit(as.character(x),","))),ncol=3),maxColorValue=255)
}

chromHMM_color_scheme = sapply(all_colors$V2,getrgb)
names(chromHMM_color_scheme) <- all_colors$V1
chromHMM_color_scheme
####################################
## Clustering samples
####################################
divisive.clust.tissue = readRDS(radio$cluster.tissue) 
col_dend = as.dendrogram(divisive.clust.tissue)


#######################################
## kmeans
#######################################

kmeans = read.table(file.path(rootfolder,radio$kmeans.all), h=T, sep ="\t", stringsAsFactors=FALSE)

mat = as.matrix(kmeans[,1:6])

  pdf(file.path(rootfolder,"plots","clustering","simple_mat.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat,
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
dev.off()


#######################################################
####  Active TSS plotting
#######################################################
tssA = tss_intersect.dcast[tss_intersect.dcast$cmips == "Active_promoter",]
matrix = tssA

divisive.clust.gene = readRDS(radio$cluster.active.tss)
row_dend = as.dendrogram(divisive.clust.gene)

nb_clust_gene =5
row_dend.cut = cutree(divisive.clust.gene,k=nb_clust_gene)
  

  ComplexHeatmap::Heatmap(matrix,
                          col = chromHMM_color_scheme,
                          cluster_columns=col_dend,
                          row_order = divisive.clust.gene$order,
                          show_row_names = FALSE,
                          jitter = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
  
  pdf(paste0("tss_cluster_5",".pdf"),width=8, height=8)
  
  ComplexHeatmap::Heatmap(matrix,
                          col = chromHMM_color_scheme,
                          cluster_columns=col_dend,
                          row_order = divisive.clust.gene$order,
                          show_row_names = FALSE,
                          jitter = TRUE,
                          border=TRUE,
                          split = data.frame(row_dend.cut),
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          use_raster = TRUE, raster_device = "png",raster_quality =1)

dev.off()

#######################################################
####  Biv TSS plotting
#######################################################
biv =  tss_intersect.dcast[tss_intersect.dcast$cmips == "Poised_promoter",]
matrix = biv

divisive.clust.gene = readRDS(radio$cluster.biv.tss)
row_dend = as.dendrogram(divisive.clust.gene)

nb_clust_gene = 6
row_dend.cut = cutree(divisive.clust.gene,k=nb_clust_gene)

head(divisive.clust.gene)

pdf(paste0("test",".pdf"),width=8, height=8)
ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        use_raster = TRUE, raster_device = "png",raster_quality =1)

ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        cluster_columns=col_dend,
                        #row_order = divisive.clust.gene$order,
                        show_row_dend = TRUE,
                        cluster_rows = row_dend,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        use_raster = TRUE, raster_device = "png",raster_quality =1)



ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        cluster_columns=col_dend,
                        row_order = divisive.clust.gene$order,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        use_raster = TRUE, raster_device = "png",raster_quality =1)

##This is the good plot
pdf(paste0("biv_cluster_6",".pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        cluster_columns=col_dend,
                        row_order = divisive.clust.gene$order,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        border=TRUE,
                        split = data.frame(row_dend.cut),
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        use_raster = TRUE, raster_device = "png",raster_quality =1)
dev.off()

#######################################################
####  Weak TSS plotting
#######################################################

weak =  tss_intersect.dcast[tss_intersect.dcast$cmips == "Weak_promoter",]
matrix = weak

divisive.clust.gene = readRDS(radio$cluster.weak.tss)
nb_clust_gene = 10

row_dend = as.dendrogram(divisive.clust.gene)
row_dend.cut = cutree(as.hclust(divisive.clust.gene),k=nb_clust_gene)



pdf(paste0("weak_clustering",".pdf"),width=8, height=8)  

ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        cluster_columns=col_dend,
                        row_order = order(row_dend.cut),
                        row_split = 5,
                        #show_row_dend = FALSE,
                        cluster_rows = row_dend,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        border = TRUE,
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        use_raster = TRUE, raster_device = "png",raster_quality =1)

dev.off()

#######################################################
####  Full plotting
#######################################################
matrix = tss_intersect.dcast
divisive.clust.gene = readRDS(radio$cluster.full.hclust)

nb_clust_gene = 10
row_dend.cut = cutree(divisive.clust.gene,k=nb_clust_gene)

ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        cluster_columns=col_dend,
                        #row_order = divisive.clust.gene$order,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        use_raster = TRUE, raster_device = "png",raster_quality =1)

ComplexHeatmap::Heatmap(matrix,
                        col = chromHMM_color_scheme,
                        cluster_columns=col_dend,
                        row_order = divisive.clust.gene$order,
                        show_row_names = FALSE,
                        jitter = TRUE,
                        border=TRUE,
                        split = data.frame(row_dend.cut),
                        heatmap_legend_param=list(title = "ChromHMM annot"),
                        cluster_row_slices = FALSE,
                        use_raster = TRUE, raster_device = "png",raster_quality = 1)

#######################################################
####  Encode colored matrix for niceness
#######################################################

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
