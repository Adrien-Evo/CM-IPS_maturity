library(yaml)
library(ggthemes)
library(tidyr)
require(scales)
#library(GenomicRanges)
#library(bedr)
library(reshape2)
library(ggplot2)
library(dplyr)

radio = read_yaml("/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/radiofile.yml")

outDir=paste0(radio$viz,"/exploration-of-annotated-states_tss_distance")
dir.create(outDir)
setwd(outDir)


##Discretize distance to TSS with the cut function


barplot_tss_encode <- function(encode,nb_encode_states,plotname){
  
  encode$tss.dis = cut(abs(encode$Distance.to.TSS), c(0,500,1000,2000,3000,4000,5000,10000,Inf),include.lowest=TRUE,labels= c("0-0.5","0.5-1","1-2","2-3","3-4","4-5","5-10",">10"))
  
  #Simple overlap of TSS:
  overlap=(encode$End- encode$Start)/2 > abs(encode$Distance.to.TSS)
  encode_overlap = encode[overlap,]
  encode_overlap$tss.dis <-rep("overlap",length(encode_overlap$tss.dis))
  
  overlap_melted = melt(data = encode_overlap, id.vars = "tss.dis", measure.vars = c("name"))
  overlap_toplot = overlap_melted %>% group_by(tss.dis) %>% count(value)
  head(overlap_toplot)
  
  ##Color scheme
  if(nb_encode_states == 15){
    color = read.table(radio$color.mapping.encode.15,sep="\t")
  }else if(nb_encode_states == 18){
    color = read.table(radio$color.mapping.encode.18,sep="\t")
  }
  
  getrgb <-function(x){
    rgb(matrix(as.vector(unlist(strsplit(as.character(x),","))),ncol=3),maxColorValue=255)
  }
  
  chromHMM_color_scheme = sapply(color$V2,getrgb)
  names(chromHMM_color_scheme) <- color$V1
  chromHMM_color_scheme
  
  
  melted_intersect=melt(data = encode, id.vars = "tss.dis", measure.vars = c("name"))
  head(melted_intersect)
  toplot = melted_intersect %>% group_by(tss.dis) %>% count(value)
  head(toplot)
  toplot=rbind(toplot,overlap_toplot)
  toplot$tss.dis <- factor(toplot$tss.dis,levels= c("overlap","0-0.5","0.5-1","1-2","2-3","3-4","4-5","5-10",">10"))
  ####Rearrange levels for better plot organization
  
  if(nb_encode_states == 15){
    toplot$value <- factor(toplot$value,levels = rev(c("TssA","TssAFlnk","TssBiv","BivFlnk","EnhG","EnhA","EnhBiv",
                                                       "Tx","ReprPC","ZNF_Rpts","Het",
                                                       "Quies")))
  }else if(nb_encode_states == 18){
    toplot$value <- factor(toplot$value,levels = rev(c("TssA","TssAFlnk","TssBiv","EnhG","EnhA",
                                                       "EnhWk","EnhBiv","Tx","ReprPC","ZNF_Rpts","Het",
                                                       "Quies")))
  }
  
  theme_set(theme_classic())
  g <- ggplot(toplot) + 
    geom_bar(stat="identity", aes(fill=value, x=tss.dis, y=n),colour="black", width = 0.95,position = position_fill()) + 
    coord_flip()  +   labs(fill="ENCODE ChromHMM annotation") + theme_tufte()+   scale_fill_manual(values = chromHMM_color_scheme)+
    xlab("Distance to TSS (kb)") + ylab("Fraction") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(paste0(plotname))
  g
  
   ggsave(g,filename = paste0("ENCODE_Tss_fraction_", plotname,".png"),width = 8, height = 6 )
    
  gg <- ggplot(toplot) +
    geom_bar(stat="identity", aes(fill=value, x=tss.dis, y=n),colour="black", width = 0.95) + 
    labs(fill="ENCODE ChromHMM annotation") + theme_tufte()+   scale_fill_manual(values = chromHMM_color_scheme)+ xlab("Distance to TSS (kb)") +
    ylab("Count") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(paste0(plotname)) + scale_y_continuous(labels = scales::comma)
  gg
  ggsave(gg,filename = paste0("ENCODE_Tss_count_", plotname,".png"),width = 8, height = 6 )
}

###For new tss distances
vent.R = read.table(radio$HRT.VENT.R.a,sep="\t", header=T, quote="", fill = TRUE)
barplot_tss_encode(vent.R,18,"VENT.R")

HRT.FET = read.table(radio$HRT.FET.a,sep="\t", header=T, quote="", fill = TRUE)
barplot_tss_encode(HRT.FET,15,"FET.HRT")

ENCODE.IPS = read.table(radio$IPS.a,sep="\t", header=T, quote="", fill = TRUE)
barplot_tss_encode(ENCODE.IPS,18,"ENCODE.IPS")


encode = vent.R
