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

outDir=paste0(radio$viz,"/barplot_of_ENCODE_states_for_each_CMIPS_chromHMM_states")
dir.create(outDir)
setwd(outDir)

barplot_intersect <- function(intersect_cmips_encode,nb_encode_states,plotname){
  
    
  colnames(intersect_cmips_encode) <-c("chr1","start1","end1","cmips","chr2","start2","end2","encode")

  
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
  
  melted_intersect=melt(data = intersect_cmips_encode, id.vars = "cmips", measure.vars = c("encode"))
  head(melted_intersect)
  toplot = melted_intersect %>% group_by(cmips) %>% count(value)
  head(toplot)
  ####Rearrange levels for better plot organization
  toplot$cmips <- factor(toplot$cmips,levels = rev(c("Active_promoter","Weak_promoter","Poised_promoter","Strong_enhancer","Poised_enhancer",
                                                                      "Polycomb_repressed","Heterochrom_lowsignal")))
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
    geom_bar(stat="identity", aes(fill=value, x=cmips, y=n),colour="black", width = 0.95,position = position_fill()) + 
    coord_flip()  +   labs(fill="ENCODE ChromHMM annotation") + theme_tufte()+   scale_fill_manual(values = chromHMM_color_scheme)+ xlab("CM-IPS") + ylab("Fraction") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(paste0("CM-IPS vs ENCODE ", plotname))
  g
  
  ggsave(g,filename = paste0(plotname,"_fraction.png"),width = 8, height = 6 )
  
  gg <- ggplot(toplot) + 
    geom_bar(stat="identity", aes(fill=value, x=cmips, y=n),colour="black", width = 0.95) + 
    labs(fill="ENCODE ChromHMM annotation") + theme_tufte()+   scale_fill_manual(values = chromHMM_color_scheme)+ xlab("CM-IPS") + ylab("Count") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(paste0("CM-IPS vs ENCODE ", plotname))+ scale_y_continuous(labels = scales::comma)
  gg
  
  ggsave(gg,filename = paste0(plotname,"_count.png"),width = 8, height = 6 )
}

intersect_cmips_encode=read.table(radio$intersect.all.CMIPS.FET,h=F,sep = "\t",stringsAsFactors=F)

barplot_intersect(intersect_cmips_encode,15,"FETAL")

intersect_cmips_encode=read.table(radio$intersect.all.CMIPS.VENT.R,h=F,sep = "\t",stringsAsFactors=F)

barplot_intersect(intersect_cmips_encode,18,"VENT.R")
  
intersect_cmips_encode=read.table(radio$intersect.all.CMIPS.FET.tss_windows_1kb,h=F,sep = "\t",stringsAsFactors=F)

barplot_intersect(intersect_cmips_encode,15,"FETAL_1kb_TSS")

intersect_cmips_encode=read.table(radio$intersect.all.CMIPS.VENT.R.tss_windows_1kb,h=F,sep = "\t",stringsAsFactors=F)

barplot_intersect(intersect_cmips_encode,18,"VENT.R_1kb_TSS")

####Tss CCDS

intersect_cmips_encode=read.table(radio$intersect.all.CMIPS.FET.tss.ccds,h=F,sep = "\t",stringsAsFactors=F)

barplot_intersect(intersect_cmips_encode,15,"FETAL_TSS_CCDS")

intersect_cmips_encode=read.table(radio$intersect.all.CMIPS.VENT.R.tss.ccds,h=F,sep = "\t",stringsAsFactors=F)

barplot_intersect(intersect_cmips_encode,18,"VENT.R_TSS_CCDS")


### TSS CCDS just looking at he base level
tss_intersect = read.table(radio$intersect.tss,h=F,sep = "\t", stringsAsFactors = FALSE)

tss_intersect.dcast = dcast(tss_intersect,V1~V2)


tss_intersect.dcast[sapply(tss_intersect.dcast, is.character)] <- lapply(tss_intersect.dcast[sapply(tss_intersect.dcast, is.character)], 
                                                                         as.factor)

head(tss_intersect.dcast)

test = tss_intersect.dcast[,c(1,1,1,4,1,1,1,8)]
barplot_intersect(test,18,"VENT.R_TSS_CCDS_1bp")

test = tss_intersect.dcast[,c(1,1,1,4,1,1,1,6)]
barplot_intersect(test,18,"HRT.FET_TSS_CCDS_1bp")

test = tss_intersect.dcast[,c(1,1,1,4,1,1,1,9)]
barplot_intersect(test,18,"IPS_TSS_CCDS_1bp")

sum(is.na(test$cmips))
