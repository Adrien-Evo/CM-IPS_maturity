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



barplot_tss_cmips <- function(cmips,plotname){
##Discretize distance to TSS with the cut function

cmips$tss.dis = cut(abs(cmips$Distance.to.TSS), c(0,500,1000,2000,3000,4000,5000,10000,Inf),include.lowest=TRUE,labels= c("0-0.5","0.5-1","1-2","2-3","3-4","4-5","5-10",">10"))

#Simple overlap of TSS:
overlap=(cmips$End- cmips$Start)/2 > abs(cmips$Distance.to.TSS)
cmips_overlap = cmips[overlap,]
cmips_overlap$tss.dis <-rep("overlap",length(cmips_overlap$tss.dis))

overlap_melted = melt(data = cmips_overlap, id.vars = "tss.dis", measure.vars = c("name"))
overlap_toplot = overlap_melted %>% group_by(tss.dis) %>% count(value)
head(overlap_toplot)

##Color scheme
color = read.table(radio$color.mapping.cmips,sep="\t")

getrgb <-function(x){
  rgb(matrix(as.vector(unlist(strsplit(as.character(x),","))),ncol=3),maxColorValue=255)
}

chromHMM_color_scheme = sapply(color$V2,getrgb)
names(chromHMM_color_scheme) <- color$V1
chromHMM_color_scheme

# 
melted_intersect=melt(data = cmips, id.vars = "tss.dis", measure.vars = c("name"))
head(melted_intersect)
toplot = melted_intersect %>% group_by(tss.dis) %>% count(value)
head(toplot)
toplot=rbind(toplot,overlap_toplot)
####Rearrange levels for better plot organization
toplot$value <- factor(toplot$value,levels = rev(c("Active_promoter","Weak_promoter","Poised_promoter","Strong_enhancer","Poised_enhancer",
                                                                    "Polycomb_repressed","Heterochrom_lowsignal","Repetitive_CNV")))
toplot$tss.dis <- factor(toplot$tss.dis,levels= c("overlap","0-0.5","0.5-1","1-2","2-3","3-4","4-5","5-10",">10"))


theme_set(theme_classic())
g <- ggplot(toplot) + 
  geom_bar(stat="identity", aes(fill=value, x=tss.dis, y=n),colour="black", width = 0.95,position = position_fill()) + 
  coord_flip()  +   labs(fill="CMIPS ChromHMM annotation") + theme_tufte()+   scale_fill_manual(values = chromHMM_color_scheme)+ xlab("Distance to TSS (kb)") + 
  ylab("Fraction") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(paste0(plotname))
g

ggsave(g,filename = paste0("CMIPS_Tss_fraction_", plotname,".png"),width = 8, height = 6 )

gg <- ggplot(toplot) + 
  geom_bar(stat="identity", aes(fill=value, x=tss.dis, y=n),colour="black", width = 0.95) + 
  labs(fill="CMIPS ChromHMM annotation") + theme_tufte()+   scale_fill_manual(values = chromHMM_color_scheme)+ xlab("Distance to TSS (kb)") +
  ylab("Count") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ggtitle(paste0(plotname)) + scale_y_continuous(labels = scales::comma)
gg
ggsave(gg,filename = paste0("CMIPS_Tss_count_", plotname,".png"),width = 8, height = 6 )

}

##Read intersect file
cmips = read.table(radio$chromhmm_homer,sep="\t", header=T, quote="")
head(cmips)
###For new tss distances
cmips_ccds = read.table(radio$chromhmm_homer_tss_ccds,sep="\t", header=T, quote="")


barplot_tss_cmips(cmips,"gtf_ensembl")

barplot_tss_cmips(cmips_ccds,"CCDS")

  