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
colour15 = read.table(file.path(rootfolder,radio$color.mapping.encode.15),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

colours18 = read.table(file.path(rootfolder,radio$color.mapping.encode.18),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

colourcmips= read.table(file.path(rootfolder,radio$color.mapping.cmips),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
##Keeping only annotation that we have
colourcmips = colourcmips[c(1,2,3,4,5,10,11),]

all_colors = rbind(colour15,colours18,colourcmips)

getrgb <-function(x){
rgb(matrix(as.vector(unlist(strsplit(as.character(x),","))),ncol=3),maxColorValue=255)
}



#Color scheme
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

#######################################
## ENCODE legend
#######################################


ENCODE_colours = colours18[c(1,2,3,4,5,6),]
ENCODE_colours = rbind(ENCODE_colours,colour15[c(6,7,8,9,10,11,12),])
ENCODE_colours$V2 = sapply(ENCODE_colours$V2,getrgb)
#named_list_ENCODE = ENCODE_colours$V2
#names(named_list_ENCODE) <- ENCODE_colours$V1
lgd_encode = Legend(title = "ChromHMM ENCODE",labels=ENCODE_colours$V1, legend_gp = gpar(fill = ENCODE_colours$V2),border = "black")

#######################################
## ChromHMM CM-iPS legend
#######################################

colourcmips$V2 = sapply(colourcmips$V2,getrgb)
lgd_cmips = Legend(title = "ChromHMM CM-iPS",labels=colourcmips$V1, legend_gp = gpar(fill = colourcmips$V2),border = "black")


lgd = packLegend(lgd_encode,lgd_cmips)
#######################################
## KMEANS data
#######################################

kmeans = read.table(file.path(rootfolder,radio$kmeans.all), h=T, sep ="\t", stringsAsFactors=FALSE,check.names=FALSE)


mat = as.matrix(kmeans[,1:6])


pdf(file.path(rootfolder,"plots","clustering","all.pdf"),width=8, height=8)

dd <- ComplexHeatmap::Heatmap(mat,
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png",raster_quality =1)
draw(dd, annotation_legend_list = lgd)
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

dd <- ComplexHeatmap::Heatmap(mat,
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
heatmap_legend_param=list(title = "ChromHMM annot"),
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png", raster_quality =1)
draw(dd, annotation_legend_list = lgd)
dev.off()

#################################################
##### KMEANS 3
#################################################

mat = as.matrix(kmeans[order(kmeans$kmeans3),1:6])

pdf(file.path(rootfolder,"plots","clustering","kmeans3.pdf"),width=8, height=8)

dd <- ComplexHeatmap::Heatmap(mat,
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= kmeans$kmeans3[order(kmeans$kmeans3)],
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png",raster_quality =1)
draw(dd, annotation_legend_list = lgd)
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

pdf(file.path(rootfolder,"plots","clustering","kmeans3_sorted_levels.pdf"),width=8, height=8)

dd <- ComplexHeatmap::Heatmap(mat[sort_order,],
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= factor(kmeans$kmeans3[sort_order],levels = c(1,2,3)),
cluster_row_slices = FALSE,
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png",raster_quality =1)
draw(dd, annotation_legend_list = lgd)
dev.off()

#################################################
##### KMEANS 10
#################################################
mat = as.matrix(kmeans[order(kmeans$kmeans10),1:6])

pdf(file.path(rootfolder,"plots","clustering","kmeans10.pdf"),width=8, height=8)

dd <- ComplexHeatmap::Heatmap(mat,
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= kmeans$kmeans10[order(kmeans$kmeans10)],
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
gap = unit(2, "mm"),
row_title_rot = 0,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png",raster_quality =1)
draw(dd, annotation_legend_list = lgd)
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

dd <- ComplexHeatmap::Heatmap(mat[sort_order,],
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= factor(kmeans$kmeans10[sort_order],levels = 1:10),
cluster_row_slices = FALSE,
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
gap = unit(2, "mm"),
row_title_rot = 0,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png",raster_quality =1)
draw(dd, annotation_legend_list = lgd)
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

png(file.path(rootfolder,"plots","clustering","kmeans15_sorted_levels.png"),width = 640, height = 640, units = "px")

dd <- ComplexHeatmap::Heatmap(mat[sort_order,],
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= factor(kmeans$kmeans15[sort_order],levels = 1:15),
cluster_row_slices = FALSE,
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
gap = unit(2, "mm"),
row_title_rot = 0,
use_raster = TRUE, raster_device = "png",raster_quality =2)
draw(dd, annotation_legend_list = lgd)
dev.off()


pdf(file.path(rootfolder,"plots","clustering","kmeans15_sorted_levels.pdf"))

dd <- ComplexHeatmap::Heatmap(mat[sort_order,],
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= factor(kmeans$kmeans15[sort_order],levels = 1:15),
cluster_row_slices = FALSE,
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
gap = unit(2, "mm"),
row_title_rot = 0,
use_raster = FALSE, raster_device = "png",raster_quality =2)
draw(dd, annotation_legend_list = lgd)
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

dd <- ComplexHeatmap::Heatmap(mat[sort_order,],
col = chromHMM_color_scheme,
jitter = TRUE,
show_row_names = FALSE,
border = TRUE,
heatmap_legend_param=list(title = "ChromHMM annot"),
row_split= factor(kmeans$kmeans20[sort_order],levels = 1:20),
cluster_row_slices = FALSE,
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,gap = unit(2, "mm"),
row_title_rot = 0,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png",raster_quality =1)
draw(dd, annotation_legend_list = lgd)
dev.off()

#######################################################
####  Plotting for each states in CMIPS
#######################################################
for(state in colourcmips$V1){
lgd_encode = Legend(title = "ChromHMM ENCODE",labels=ENCODE_colours$V1, legend_gp = gpar(fill = ENCODE_colours$V2),border = "black")

lgd_cmips = Legend(title = "ChromHMM CM-iPS",labels=state, legend_gp = gpar(fill = colourcmips$V2[which(colourcmips$V1 == state)]),border = "black")


lgd = packLegend(lgd_encode,lgd_cmips)

select = kmeans[which(kmeans[,1] == state),]
mat = as.matrix(select[,1:6])



sort_order = order(select$singular_kmeans,
match(mat[,1],merged_HMM_order),
match(mat[,2],encode_HMM_order),
match(mat[,3],encode_HMM_order),
match(mat[,4],encode_HMM_order),
match(mat[,5],encode_HMM_order),
match(mat[,6],encode_HMM_order))


png(file.path(rootfolder,"plots","clustering",paste0(state,"_kmeans10.png")),width = 640, height = 640, units = "px")

dd <-ComplexHeatmap::Heatmap(mat[sort_order,],col = chromHMM_color_scheme,jitter = TRUE,show_row_names = FALSE,border = TRUE,heatmap_legend_param=list(title = "ChromHMM annot"),row_split= factor(select$singular_kmeans[sort_order],levels = 1:10),
  cluster_row_slices = FALSE,bottom_annotation = matrix.annotation,column_names_side = "top",column_names_rot = 45,show_heatmap_legend = FALSE,gap = unit(2, "mm"),row_title_rot = 0,
  use_raster = TRUE, raster_device = "png",raster_quality = 2)
draw(dd, annotation_legend_list = lgd)
dev.off()
}



#######################################
## ChromHMM CM-iPS legend
#######################################

lgd_cmips = Legend(title = "ChromHMM CM-iPS",labels=state, legend_gp = gpar(fill = colourcmips$V2[which(colourcmips$V1 == state)]),border = "black")


lgd = packLegend(lgd_encode,lgd_cmips)