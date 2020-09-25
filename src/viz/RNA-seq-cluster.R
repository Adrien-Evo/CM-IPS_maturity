library(biomaRt)
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


#######################################
## RNA-SEQ cluster data from Bastien
########################################

rna = read.table(file.path(rootfolder,radio$rna.cluster), h=T, sep ="\t", stringsAsFactors=FALSE,check.names=FALSE)

mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
genes=getBM(attributes = c("hgnc_symbol", "entrezgene_id","ensembl_gene_id"),filters = "hgnc_symbol", values = rna$GeneID[!duplicated(rna$GeneID)], bmHeader = TRUE, mart = mart)


for (cluster in names(table(rna$Cluster))){
    ids = rna$GeneID[which(rna$Cluster == cluster)]
    print(ids)
    genes=getBM(attributes = c("hgnc_symbol", "entrezgene_id","ensembl_gene_id","ensembl_transcript_id"),filters = "hgnc_symbol", values = ids[!duplicated(ids)], bmHeader = TRUE, mart = mart)
    select = atac[rownames(atac) %in% genes[,4],]
    mat = as.matrix(select[,5:10])
    # reorder column into CMIPS HRT.VENT.L HRT.VENT.R HRT.ATR.R IPS HRT.FET

    sort_order = order(match(mat[,1],merged_HMM_order),
    match(mat[,2],encode_HMM_order),
    match(mat[,3],encode_HMM_order),
    match(mat[,4],encode_HMM_order),
    match(mat[,5],encode_HMM_order),
    match(mat[,6],encode_HMM_order))
    mat <- mat[sort_order,]




    pdf(file.path(rootfolder,"plots","RNA-seq-integration",paste0("cluster_",cluster,"_all_sorted.pdf")),width=8, height=8)

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

    png(file.path(rootfolder,"plots","RNA-seq-integration",paste0("cluster_",cluster,"_all_sorted.png")),width = 640, height = 640, units = "px")
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


}

for(cluster in 1:9){
    go = read.table(paste0(file.path(rootfolder,radio$rna.go),cluster,".txt"), h=T, sep ="\t", stringsAsFactors=FALSE,check.names=FALSE)

    for (i in 1:length(go$geneID)){
        ids = strsplit(go$geneID[i],"/")[[1]]
        print(ids)
        genes=getBM(attributes = c("hgnc_symbol", "entrezgene_id","ensembl_gene_id","ensembl_transcript_id"),filters = "hgnc_symbol", values = ids[!duplicated(ids)], bmHeader = TRUE, mart = mart)
        select = atac[rownames(atac) %in% genes[,4],]
        mat = as.matrix(select[,5:10])
        goname=gsub(" ","-",go$Description)

        sort_order = order(match(mat[,1],merged_HMM_order),
        match(mat[,2],encode_HMM_order),
        match(mat[,3],encode_HMM_order),
        match(mat[,4],encode_HMM_order),
        match(mat[,5],encode_HMM_order),
        match(mat[,6],encode_HMM_order))
        mat <- mat[sort_order,]

        pdf(file.path(rootfolder,"plots","RNA-seq-integration","GO",paste0("cluster_",cluster,"_",go$ID[i],"_",goname,".pdf")),width=8, height=8)

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

        png(file.path(rootfolder,"plots","RNA-seq-integration","GO",paste0("cluster_",cluster,"_",go$ID[i],"_",goname,".png")),width = 640, height = 640, units = "px")
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

    }
}

