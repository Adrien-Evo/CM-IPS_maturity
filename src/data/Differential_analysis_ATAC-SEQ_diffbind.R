###Works with the DiffBind environnement

library(DiffBind)
library(ggplot2)
library(dplyr)
library(yaml)


rootfolder = "/ceph-recherche/shares/u1087/afoucal/CM-IPS_maturity/"
radio = yaml.load_file(file.path(rootfolder,"radiofile.yml"))



sample=read.csv(file.path(rootfolder,radio$atac),h=T,sep=",")

## Keeping only QC positiv samples
goodSamples=sample[which(sample$QC==1),]


###Doing the CM vs IPS analysis


###QC selection
DBsample=dba(sampleSheet=goodSamples)





##Plot folder
setwd(file.path(rootfolder,radio$viz,"DiffBind"))

#######
####### OCCUPANCY ANALYSIS
########

png("Clustering_Occupancy.png", width= 840, height = 840)
plot(DBsample)
dev.off()



###Overlap
cm_overlap <- dba.overlap(DBsample,DBsample$masks$CM ,mode=DBA_OLAP_RATE)
png("Overlap_CM.png", width= 840, height = 840)
plot(cm_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for CM ')
dev.off()


i_overlap <- dba.overlap(DBsample,DBsample$masks$i ,mode=DBA_OLAP_RATE)

png("Overlap_iPS.png", width= 840, height = 840)
plot(i_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for iPS ')
dev.off()


# ===> consensus with half looks good

DBsample_simple <- dba.peakset(DBsample, consensus=c(DBA_TISSUE), minOverlap=0.5)


DBsample_brugada <- dba.peakset(DBsample, consensus=c(DBA_TISSUE, DBA_FACTOR), minOverlap=0.5)



png("Venn_CM_i_consensus_with_brugada_peakset.png", width= 840, height = 840)

dba.plotVenn(DBsample_brugada,DBsample_brugada$masks$Consensus)

dev.off()

png("Venn_CM_i_consensus peakset.png", width= 840, height = 840)

dba.plotVenn(DBsample_simple,DBsample_simple$masks$Consensus)

dev.off()





#################### DIfferential
################## Affinity analysis
##################################################################
DBcount <- dba.count(DBsample)




png("PCA_tissue.png", width= 840, height = 840)

dba.plotPCA(DBcount,attributes=DBA_TISSUE,title="CM vs IPS")
dev.off()


png("PCA_brugada.png", width= 840, height = 840)

dba.plotPCA(DBcount,attributes=DBA_FACTOR,title="Control vs Brugada")

dev.off()


png("PCA_batch.png", width= 840, height = 840)

dba.plotPCA(DBcount,attributes=DBA_TREATMENT,title="Batch")

dev.off()


png("Clustering_Affinity.png", width= 840, height = 840)
plot(DBcount)
dev.off()


contrast <- dba.contrast(DBcount, categories=DBA_TISSUE)

CMvsIPS<- dba.analyze(contrast)

png("MA_plot_CM_IPS.png", width= 840, height = 840)

dba.plotMA(CMvsIPS)
dev.off()

png("volcano_plot_CM_IPS.png", width= 840, height = 840)
dba.plotVolcano(CMvsIPS)
dev.off()

CMvsIPS.report <- dba.report(CMvsIPS)



chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")

df <- data.frame(seqnames=seqnames(CMvsIPS.report),
  starts=start(CMvsIPS.report)-1,
  ends=end(CMvsIPS.report),
  names=c(rep(".", length(CMvsIPS.report))),
  scores=c(rep(".", length(CMvsIPS.report))),
  strands=strand(CMvsIPS.report),
  conc = CMvsIPS.report$Conc,
                 conc_CM=CMvsIPS.report$Conc_CM,
                 conc_i=CMvsIPS.report$Conc_i,
                 fold = CMvsIPS.report$Fold,
                 pval = CMvsIPS.report$"p-value",
                 FDR = CMvsIPS.report$FDR
           )

df <- df[order(df[,1],df[,2]),]
df[,4] <- paste0("peak",1:length(brugg_block_report))
write.table(df,paste0("brugg_DBA_24",i,".bed"),quote=F,row.names=F,col.names=F,sep="\t")