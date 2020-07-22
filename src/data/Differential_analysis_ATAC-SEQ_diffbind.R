###Works with the DiffBind environnement

library(DiffBind)
library(ggplot2)
library(dplyr)
library(yaml)


rootfolder = "/ceph-recherche/shares/u1087/afoucal/CM-IPS_maturity/"
radio = yaml.load_file(file.path(rootfolder,"radiofile.yml"))



sample=read.csv(file.path(rootfolder,radio$atac.sample.sheet),h=T,sep=",", stringsAsFactor = FALSE)

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

saveRDS(DBcount, file = file.path(rootfolder,radio$dba.count))

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

png("volcano_plot_CM_IPS_pval.png", width= 840, height = 840)
dba.plotVolcano(CMvsIPS, bUsePval=TRUE, th = 0.01, fold=1)
dev.off()


png("volcano_plot_CM_IPS_FDR.png", width= 840, height = 840)
dba.plotVolcano(CMvsIPS, bUsePval=FALSE, th = 0.01, fold=1)
dev.off()


CMvsIPS.report <- dba.report(CMvsIPS, th = 1,bCalled=TRUE)
##Reordering df

CMvsIPS.report <- sort(CMvsIPS.report)

df <- data.frame(seqnames=seqnames(CMvsIPS.report),
  starts=start(CMvsIPS.report)-1,
  ends=end(CMvsIPS.report),
  names=c(paste0("peak",1:length(CMvsIPS.report))),
  scores=c(rep(".", length(CMvsIPS.report))),
  strands=strand(CMvsIPS.report),
  conc = CMvsIPS.report$Conc,
                 conc_CM=CMvsIPS.report$Conc_CM,
                 conc_i=CMvsIPS.report$Conc_i,
                 fold = CMvsIPS.report$Fold,
                 pval = CMvsIPS.report$"p-value",
                 FDR = CMvsIPS.report$FDR
           )


write.table(df,file.path(rootfolder,"data","processed","DiffBind","ATAC","CM_vs_IPS_full.csv"),quote=F,row.names=F,col.names=F,sep="\t")

###Selecting based on asb(fold change) > 1 and FDR < 1%
df_iPS <- df[which(df$fold <= -1),]
df_CM <- df[which(df$fold >= 1),]

df_iPS <- df_iPS[which(df_iPS$FDR <= 0.01),]
df_CM <- df_CM[which(df_CM$FDR <= 0.01),]

df_iPS <- df_iPS[,c(1,2,3,4)]
df_CM <- df_CM[,c(1,2,3,4)]

write.table(df_iPS,file.path(rootfolder,"data","processed","DiffBind","ATAC","IPS_diff_ATAC_fdr-0.01_fold-1.csv"),quote=F,row.names=F,col.names=F,sep="\t")

write.table(df_CM,file.path(rootfolder,"data","processed","DiffBind","ATAC","CM_diff_ATAC_fdr-0.01_fold-1.csv"),quote=F,row.names=F,col.names=F,sep="\t")










################################Brugada vs COntrol

brugg_overlap <- dba.overlap(DBsample,DBsample$masks$CM & DBsample$masks$brugada ,mode=DBA_OLAP_RATE)
png("Overlap_Brugada_CM.png", width= 840, height = 840)
plot(brugg_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for CM BRUGADA')
dev.off()


control_overlap <- dba.overlap(DBsample,DBsample$masks$CM & DBsample$masks$control ,mode=DBA_OLAP_RATE)

png("Overlap_control_CM.png", width= 840, height = 840)
plot(control_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for iPS ')
dev.off()



contrast <- dba.contrast(DBcount, categories=DBA_FACTOR, block=DBA_TISSUE)

BruggvsControl<- dba.analyze(contrast,method= DBA_DESEQ2)



png("MA_plot_Brugg_Control_Block_tissu.png", width= 840, height = 840)

dba.plotMA(BruggvsControl,method=DBA_DESEQ2_BLOCK)
dev.off()

png("volcano_plot_BRUGADA_vs_CONTROL_Block_tissu_PVAL.png", width= 840, height = 840)
dba.plotVolcano(BruggvsControl, bUsePval=TRUE, th = 0.05, fold=1,method=DBA_DESEQ2_BLOCK)
dev.off()

png("volcano_plot_BRUGADA_vs_CONTROL_Block_tissu_FDR.png", width= 840, height = 840)
dba.plotVolcano(BruggvsControl, bUsePval=FALSE, th = 0.1, fold=1,method=DBA_DESEQ2_BLOCK)
dev.off()

BruggvsControl.report <- dba.report(BruggvsControl,th = 1,bCalled=TRUE, method = c(DBA_DESEQ2_BLOCK) )
##Reordering df

BruggvsControl.report <- sort(BruggvsControl.report)

df <- data.frame(seqnames=seqnames(BruggvsControl.report),
  starts=start(BruggvsControl.report)-1,
  ends=end(BruggvsControl.report),
  names=c(paste0("peak",1:length(BruggvsControl.report))),
  scores=c(rep(".", length(BruggvsControl.report))),
  strands=strand(BruggvsControl.report),
  conc = BruggvsControl.report$Conc,
                 conc_brugada=BruggvsControl.report$Conc_brugada,
                 conc_control=BruggvsControl.report$Conc_control,
                 fold = BruggvsControl.report$Fold,
                 pval = BruggvsControl.report$"p-value",
                 FDR = BruggvsControl.report$FDR
           )


write.table(df,file.path(rootfolder,"data","processed","ATAC","DiffBind","BRUGADA_vs_CONTROL_full.csv"),quote=F,row.names=F,col.names=F,sep="\t")

###Selecting based on asb(fold change) > 1 and FDR < 10 %
df_CONTROL <- df[which(df$fold <= -1),]
df_CONTROL <- df_CONTROL[which(df_CONTROL$FDR <= 0.1),]
df_CONTROL <- df_CONTROL[,c(1,2,3,4)]

df_BRUG <- df[which(df$fold >= 1),]
df_BRUG <- df_BRUG[which(df_BRUG$FDR <= 0.1),]
df_BRUG <- df_BRUG[,c(1,2,3,4)]

write.table(df_CONTROL,file.path(rootfolder,"data","processed","ATAC","DiffBind","CONTROL_diff_ATAC_block_fdr-0.1_fold-1.csv"),quote=F,row.names=F,col.names=F,sep="\t")

write.table(df_BRUG,file.path(rootfolder,"data","processed","ATAC","DiffBind","BRUG_diff_ATAC_block_fdr-0.1_fold-1.csv"),quote=F,row.names=F,col.names=F,sep="\t")

