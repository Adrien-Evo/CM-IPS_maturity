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
colourcmips$V2 <- sapply(colourcmips$V2,getrgb)

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

lgd_cmips = Legend(title = "ChromHMM CM-iPS",labels=colourcmips$V1, legend_gp = gpar(fill = colourcmips$V2),border = "black")


lgd = packLegend(lgd_encode,lgd_cmips)
#######################################
## KMEANS data
#######################################

atac = read.table(file.path(rootfolder,radio$atac.chrom.kmeans), h=T, sep ="\t", stringsAsFactors=FALSE,check.names=FALSE)

###############################
#########Looking at peaks in one group and note the other. No clustering.

##Sizing
sum(!is.na(atac$ATAC.CMiPS))                                                                                                                                                                              

sum(!is.na(atac$ATAC.iPS))

##Sizing
sum(!is.na(atac$ATAC.Brug))                                                                                                                                                                                    

sum(!is.na(atac$ATAC.Control))




mat = as.matrix(atac[!is.na(atac$ATAC.CMiPS),5:10])
# reorder column into CMIPS HRT.VENT.L HRT.VENT.R HRT.ATR.R IPS HRT.FET

sort_order = order(match(mat[,1],merged_HMM_order),
match(mat[,2],encode_HMM_order),
match(mat[,3],encode_HMM_order),
match(mat[,4],encode_HMM_order),
match(mat[,5],encode_HMM_order),
match(mat[,6],encode_HMM_order))
mat <- mat[sort_order,]

pdf(file.path(rootfolder,"plots","Atac-seq-integration/clustering","all_sorted_CM.pdf"),width=8, height=8)

dd <- ComplexHeatmap::Heatmap(mat,
col = chromHMM_color_scheme,
show_row_names = FALSE,
heatmap_legend_param=list(title = "ChromHMM annot"),
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png", raster_quality =1)
draw(dd, annotation_legend_list = lgd)
dev.off()


mat = as.matrix(atac[!is.na(atac$ATAC.iPS),5:10])
# reorder column into CMIPS HRT.VENT.L HRT.VENT.R HRT.ATR.R IPS HRT.FET

sort_order = order(match(mat[,1],merged_HMM_order),
match(mat[,2],encode_HMM_order),
match(mat[,3],encode_HMM_order),
match(mat[,4],encode_HMM_order),
match(mat[,5],encode_HMM_order),
match(mat[,6],encode_HMM_order))
mat <- mat[sort_order,]

pdf(file.path(rootfolder,"plots","Atac-seq-integration/clustering","all_sorted_iPS.pdf"),width=8, height=8)

dd <- ComplexHeatmap::Heatmap(mat,
col = chromHMM_color_scheme,
show_row_names = FALSE,
heatmap_legend_param=list(title = "ChromHMM annot"),
bottom_annotation = matrix.annotation,
column_names_side = "top",
column_names_rot = 45,
show_heatmap_legend = FALSE,
use_raster = TRUE, raster_device = "png", raster_quality =1)
draw(dd, annotation_legend_list = lgd)
dev.off()

###Testons singular kmeans sur poised promoter and ips. It's super similar to singular keamsn without selecting atac-seq peaks
for(state in colourcmips$V1){
lgd_encode = Legend(title = "ChromHMM ENCODE",labels=ENCODE_colours$V1, legend_gp = gpar(fill = ENCODE_colours$V2),border = "black")

lgd_cmips = Legend(title = "ChromHMM CM-iPS",labels=state, legend_gp = gpar(fill = colourcmips$V2[which(colourcmips$V1 == state)]),border = "black")


lgd = packLegend(lgd_encode,lgd_cmips)

select = atac[which(atac[,5] == state),]
select = select[!is.na(select$ATAC.iPS),]
mat = as.matrix(select[,5:10])



sort_order = order(select$singular_kmeans,
match(mat[,1],merged_HMM_order),
match(mat[,2],encode_HMM_order),
match(mat[,3],encode_HMM_order),
match(mat[,4],encode_HMM_order),
match(mat[,5],encode_HMM_order),
match(mat[,6],encode_HMM_order))


png(file.path(rootfolder,"plots","Atac-seq-integration","clustering",paste0(state,"_iPSpeaks_singular_kmeans10.png")),width = 640, height = 640, units = "px")

dd <-ComplexHeatmap::Heatmap(mat[sort_order,],col = chromHMM_color_scheme,show_row_names = FALSE,border = TRUE,heatmap_legend_param=list(title = "ChromHMM annot"),row_split= factor(select$singular_kmeans[sort_order],levels = 1:10),
  cluster_row_slices = FALSE,bottom_annotation = matrix.annotation,column_names_side = "top",column_names_rot = 45,show_heatmap_legend = FALSE,gap = unit(2, "mm"),row_title_rot = 0,
  use_raster = TRUE, raster_device = "png",raster_quality = 2)
draw(dd, annotation_legend_list = lgd)
dev.off()
write.table(rownames(mat),file.path(rootfolder,"data","processed","ATAC","intersect-sets", paste0(state,"_iPSpeaks.txt")),col.names = F,row.names = F,quote=F)

}

###CMCMCMCMCMCMCMCMCMCMCMCMCMC
###Testons singular kmeans sur poised promoter and CM. It's super similar to singular keamsn without selecting atac-seq peaks
for(state in colourcmips$V1){
lgd_encode = Legend(title = "ChromHMM ENCODE",labels=ENCODE_colours$V1, legend_gp = gpar(fill = ENCODE_colours$V2),border = "black")

lgd_cmips = Legend(title = "ChromHMM CM-iPS",labels=state, legend_gp = gpar(fill = colourcmips$V2[which(colourcmips$V1 == state)]),border = "black")


lgd = packLegend(lgd_encode,lgd_cmips)

select = atac[which(atac[,5] == state),]
select = select[!is.na(select$ATAC.CMiPS),]
mat = as.matrix(select[,5:10])



sort_order = order(select$singular_kmeans,
match(mat[,1],merged_HMM_order),
match(mat[,2],encode_HMM_order),
match(mat[,3],encode_HMM_order),
match(mat[,4],encode_HMM_order),
match(mat[,5],encode_HMM_order),
match(mat[,6],encode_HMM_order))


png(file.path(rootfolder,"plots","Atac-seq-integration","clustering",paste0(state,"_CMpeaks_singular_kmeans10.png")),width = 640, height = 640, units = "px")

dd <-ComplexHeatmap::Heatmap(mat[sort_order,],col = chromHMM_color_scheme,show_row_names = FALSE,border = TRUE,heatmap_legend_param=list(title = "ChromHMM annot"),row_split= factor(select$singular_kmeans[sort_order],levels = 1:10),
  cluster_row_slices = FALSE,bottom_annotation = matrix.annotation,column_names_side = "top",column_names_rot = 45,show_heatmap_legend = FALSE,gap = unit(2, "mm"),row_title_rot = 0,
  use_raster = TRUE, raster_device = "png",raster_quality = 2)
draw(dd, annotation_legend_list = lgd)
dev.off()

write.table(rownames(mat),file.path(rootfolder,"data","processed","ATAC","intersect-sets", paste0(state,"_CMpeaks.txt")),col.names = F,row.names = F,quote=F)
}


