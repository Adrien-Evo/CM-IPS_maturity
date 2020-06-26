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
##Keeping only annotation that we have
colourcmips = colourcmips[c(1,2,3,4,5,10,11),]

all_colors = rbind(colour15,colours18,colourcmips)

getrgb <-function(x){
  rgb(matrix(as.vector(unlist(strsplit(as.character(x),","))),ncol=3),maxColorValue=255)
}

chromHMM_color_scheme = sapply(all_colors$V2,getrgb)
sapply(colourcmips$V2,getrgb)

names(chromHMM_color_scheme) <- all_colors$V1
chromHMM_color_scheme

#######################################################
####  Encode colored matrix for niceness
#######################################################

ENCODE_samps = c("HRT.ATR.R","HRT.VENT.L","HRT.VENT.R","IPS","HRT.FET")  
ENCODE.col = "#377eb8"
CM_IPS.col = "#e41a1c"

fetal.col = "#c51b7d"
adult.col = "#4d9221"
ips.col = "#f1b6da"

right.atrium.col = "#66c2a5"
left.ventricule.col = "#8da0cb"
right.ventricule.col ="#a6d854"
brain.col = "#e5c494"
heart.col = "#e41a1c"
#Colors from 8-class Set2 R colorbrewer
############# Source
source.cols = c(CM_IPS.col,rep(ENCODE.col,length(ENCODE_samps)))
names(source.cols ) <- c("CM-iPS",rep("ENCODE",length(ENCODE_samps)))
source.cols

###Maturity 
stage.cols = c(CM_IPS.col,adult.col,adult.col,adult.col,ips.col,fetal.col)
names(stage.cols) <- c("CM-iPS","Mature","Mature","Mature","iPS","Fetal")

tissue.cols = c(CM_IPS.col,right.atrium.col,left.ventricule.col,right.ventricule.col, ips.col,fetal.col)

tissue.cols

rlab=t(cbind(source.cols,stage.cols,tissue.cols))

rownames(rlab)=c("Data origin","Cell type", "Tissue type")

matrix.annotation = HeatmapAnnotation(Source = names(source.cols ),Maturity = names(stage.cols ), col = list(Source = source.cols,Maturity = stage.cols ))

###############################################################
### Levels order
###############################################################
merged_HMM_order = c("Active_promoter","Weak_promoter","Poised_promoter","Strong_enhancer","Poised_enhancer","Polycomb_repressed","Heterochrom_lowsignal")
encode_HMM_order = c("TssA","TssAFlnk","TssBiv","EnhG","EnhA","EnhWk","EnhBiv","BivFlnk","ReprPC","Het")

####################################
## Clustering samples
####################################
#divisive.gower.diana 
lgd = Legend(labels=colourcmips$V1, legend_gp = colourcmips)
#ENCODE legend
all_colors = rbind(colour15,colours18)
all_colors = unique
ENCODE_leg = 
#######################################
## KMEANS data
#######################################

kmeans = read.table(file.path(rootfolder,radio$kmeans.all), h=T, sep ="\t", stringsAsFactors=FALSE,check.names=FALSE)


mat = as.matrix(kmeans[,1:6])


pdf(file.path(rootfolder,"plots","clustering","all.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat,
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          heatmap_legend_param=lgd,
                          bottom_annotation = matrix.annotation,
                          column_names_side = "top",
                          column_names_rot = 45,
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
dev.off()



### Levels ordering

mat = as.matrix(kmeans[,1:6])
# reorder column into CMIPS HRT.VENT.L HRT.VENT.R HRT.ATR.R IPS HRT.FET

sort_order = order(match(mat[,1],merged_HMM_order),
                            match(mat[,2],encode_HMM_order),
                            match(mat[,3],encode_HMM_order),
                            match(mat[,4],encode_HMM_order),
                            match(mat[,5],encode_HMM_order),
                            match(mat[,6],encode_HMM_order))
mat <- mat[sort_order,]

pdf(file.path(rootfolder,"plots","clustering","all_sorted_levels.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat,
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
dev.off()



#################################################
##### KMEANS 3
#################################################

mat = as.matrix(kmeans[order(kmeans$kmeans3),1:6])

pdf(file.path(rootfolder,"plots","clustering","kmeans3.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat,
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          row_split= kmeans$kmeans3[order(kmeans$kmeans3)],
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
dev.off()



### Levels ordering
mat = as.matrix(kmeans[,1:6])

sort_order = order(kmeans$kmeans3,
                            match(mat[,1],merged_HMM_order),
                            match(mat[,2],encode_HMM_order),
                            match(mat[,3],encode_HMM_order),
                            match(mat[,4],encode_HMM_order),
                            match(mat[,5],encode_HMM_order),
                            match(mat[,6],encode_HMM_order))
#mat <- mat[sort_order,]

pdf(file.path(rootfolder,"plots","clustering","kmeans3_sorted_levels.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat[sort_order,],
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          row_split= factor(kmeans$kmeans3[sort_order],levels = c(1,2,3)),
                          cluster_row_slices = FALSE,
                          use_raster = FALSE, raster_device = "png",raster_quality =1)
dev.off()

#################################################
##### KMEANS 10
#################################################

mat = as.matrix(kmeans[order(kmeans$kmeans10),1:6])

pdf(file.path(rootfolder,"plots","clustering","kmeans10.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat,
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          row_split= kmeans$kmeans10[order(kmeans$kmeans10)],
                          use_raster = TRUE, raster_device = "png",raster_quality =1)
dev.off()


### Levels ordering
mat = as.matrix(kmeans[,1:6])

sort_order = order(kmeans$kmeans10,
                            match(mat[,1],merged_HMM_order),
                            match(mat[,2],encode_HMM_order),
                            match(mat[,3],encode_HMM_order),
                            match(mat[,4],encode_HMM_order),
                            match(mat[,5],encode_HMM_order),
                            match(mat[,6],encode_HMM_order))

pdf(file.path(rootfolder,"plots","clustering","kmeans10_sorted_levels.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat[sort_order,],
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          row_split= factor(kmeans$kmeans10[sort_order],levels = 1:10),
                          cluster_row_slices = FALSE,
                          use_raster = FALSE, raster_device = "png",raster_quality =1)
dev.off()

#################################################
##### KMEANS 15
#################################################



### Levels ordering
mat = as.matrix(kmeans[,1:6])

sort_order = order(kmeans$kmeans15,
                            match(mat[,1],merged_HMM_order),
                            match(mat[,2],encode_HMM_order),
                            match(mat[,3],encode_HMM_order),
                            match(mat[,4],encode_HMM_order),
                            match(mat[,5],encode_HMM_order),
                            match(mat[,6],encode_HMM_order))

pdf(file.path(rootfolder,"plots","clustering","kmeans15_sorted_levels.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat[sort_order,],
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          row_split= factor(kmeans$kmeans15[sort_order],levels = 1:15),
                          cluster_row_slices = FALSE,
                          use_raster = FALSE, raster_device = "png",raster_quality =1)
dev.off()



#################################################
##### KMEANS 20
#################################################



### Levels ordering
mat = as.matrix(kmeans[,1:6])

sort_order = order(kmeans$kmeans20,
                            match(mat[,1],merged_HMM_order),
                            match(mat[,2],encode_HMM_order),
                            match(mat[,3],encode_HMM_order),
                            match(mat[,4],encode_HMM_order),
                            match(mat[,5],encode_HMM_order),
                            match(mat[,6],encode_HMM_order))

pdf(file.path(rootfolder,"plots","clustering","kmeans20_sorted_levels.pdf"),width=8, height=8)

ComplexHeatmap::Heatmap(mat[sort_order,],
                          col = chromHMM_color_scheme,
                          jitter = TRUE,
                          show_row_names = FALSE,
                          border = TRUE,
                          heatmap_legend_param=list(title = "ChromHMM annot"),
                          row_split= factor(kmeans$kmeans20[sort_order],levels = 1:20),
                          cluster_row_slices = FALSE,
                          use_raster = FALSE, raster_device = "png",raster_quality =1)
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

