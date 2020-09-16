###Works with the DiffBind environnement

library(DiffBind)
library(ggplot2)
library(dplyr)
library(yaml)


rootfolder = "/ceph-recherche/shares/u1087/afoucal/CM-IPS_maturity/"
radio = yaml.load_file(file.path(rootfolder,"radiofile.yml"))



sample=read.csv(file.path(rootfolder,radio$atac.sample.sheet),h=T,sep=",", stringsAsFactor = FALSE)

sample$Condition <- paste0(sample$Tissue,"_",sample$Factor)
## Keeping only QC positiv samples
goodSamples=sample[which(sample$QC==1),]


###Doing the CM vs IPS analysis


###QC selection
DBsample=dba(sampleSheet=goodSamples)





##Plot folder
setwd(file.path(rootfolder,radio$viz,"DiffBind"))

#######
####### OCCUPANCY ANALYSIS
#######

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


DBsample_brugada <- dba.peakset(DBsample, consensus=c(DBA_CONDITION), minOverlap=0.33)



png("Venn_CM_i_consensus_with_brugada_peakset_33.png", width= 840, height = 840)

dba.plotVenn(DBsample_brugada,DBsample_brugada$masks$Consensus)

dev.off()

png("Venn_CM_i_consensus peakset.png", width= 840, height = 840)

dba.plotVenn(DBsample_simple,DBsample_simple$masks$Consensus)

dev.off()

#Venn_Data
venn = dba.overlap(DBsample_brugada, DBsample_brugada$masks$Consensus)

write_report_occupancy <- function(report,reportname){

df <- data.frame(seqnames=seqnames(report),
    starts=start(report)-1,
    ends=end(report),
    names=c(paste0(reportname,"_peak",names(report))),
    scores=c(rep("0", length(report))),
    strands=c(rep(".", length(report)))
            )
  write.table(df,file.path(rootfolder,"data","processed","ATAC","DiffBind",paste0(reportname,"_Occupancy.csv")),quote=F,row.names=F,col.names=F,sep="\t")
}


write_report_occupancy(venn$onlyA,"CM_Brugada")
write_report_occupancy(venn$onlyB,"CM_Control")
write_report_occupancy(venn$onlyC,"iPS_Control")
write_report_occupancy(venn$onlyD,"iPS_Brugada")
write_report_occupancy(venn$BandC,"Control")
write_report_occupancy(venn$AandD,"Brugada")


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
###########Function to write the report. Here fold = fold_cond1 -fold_cond2
write_report <- function(report,reportname,cond1,cond2,FDR,fold){


  df <- data.frame(chr=seqnames(report),
    starts=start(report)-1,
    ends=end(report),
    names=c(paste0("peak",1:length(report))),
    scores=c(rep("0", length(report))),
    strands=c(rep(".", length(report))),
    conc = report$Conc,
    conc_cond1=mcols(report)[,2],
    conc_cond2=mcols(report)[,3],
    fold =report$Fold,
    pval = report$"p-value",
    FDR = report$FDR
            )
colnames(df)[8] <- paste0(cond1,"_conc")
colnames(df)[9] <- paste0(cond2,"_conc")

  write.table(df,file.path(rootfolder,"data","processed","ATAC","DiffBind",paste0(reportname,"_full.csv")),quote=F,col.names=T, row.names = F,sep="\t")

  temp = df$names
  temp[which(df$fold <= 0)] <- paste0(cond2,"_")
  temp[which(df$fold > 0)] <- paste0(cond1,"_")

  df$names=paste0(temp,df$names)
  ###Selecting based on asb(fold change) > 1 and FDR < 1%
  df_select <- df[which(abs(df$fold) >= fold),]

  df_select <- df_select[which(df_select$FDR <= FDR),]


  write.table(df_select[,c(1,2,3,4)],file.path(rootfolder,"data","processed","ATAC","DiffBind",paste0(reportname,"_fdr-",FDR,"_fold-",fold,".csv")),quote=F,row.names=F,col.names=F,sep="\t")

}

write_report(CMvsIPS.report,"CM_vs_IPS","CM","iPS",0.01,1)


################################Brugada vs COntrol

brugg_overlap <- dba.overlap(DBsample,DBsample$masks$CM & DBsample$masks$brugada ,mode=DBA_OLAP_RATE)
png("Overlap_Brugada_CM.png", width= 840, height = 840)
plot(brugg_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for CM BRUGADA')
dev.off()


control_overlap <- dba.overlap(DBsample,DBsample$masks$CM & DBsample$masks$control ,mode=DBA_OLAP_RATE)

png("Overlap_control_CM.png", width= 840, height = 840)
plot(control_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for CM Control ')
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

write_report(BruggvsControl.report,"BRUGADA_vs_CONTROL","BRUGADA","CONTROL",0.1,1)


######################
###################### Differential analyses Brugada vs COntrol using only CM
######################

goodSamples = sample[which(sample$Tissue=="CM"),]


## Keeping only QC positiv samples
goodSamples=goodSamples[which(goodSamples$QC==1),]



###QC selection
DBsample=dba(sampleSheet=goodSamples)
DBcount <- dba.count(DBsample)


contrast <- dba.contrast(DBcount, categories=DBA_FACTOR)

BruggvsControl<- dba.analyze(contrast,method= DBA_DESEQ2)



png("MA_plot_Brugg_Control_CM.png", width= 840, height = 840)

dba.plotMA(BruggvsControl,method=DBA_DESEQ2)
dev.off()

png("volcano_plot_BRUGADA_vs_CONTROL_CM_PVAL.png", width= 840, height = 840)
dba.plotVolcano(BruggvsControl, bUsePval=TRUE, th = 0.05, fold=1,method=DBA_DESEQ2)
dev.off()

png("volcano_plot_BRUGADA_vs_CONTROL_CM_FDR.png", width= 840, height = 840)
dba.plotVolcano(BruggvsControl, bUsePval=FALSE, th = 0.1, fold=1,method=DBA_DESEQ2)
dev.off()


BruggvsControl.report <- dba.report(BruggvsControl,th = 1,bCalled=TRUE, method = c(DBA_DESEQ2) )
##Reordering df

BruggvsControl.report <- sort(BruggvsControl.report)

write_report(BruggvsControl.report,"BRUGADA_vs_CONTROL_onlyCM","BRUGADA","CONTROL",0.1,1)



