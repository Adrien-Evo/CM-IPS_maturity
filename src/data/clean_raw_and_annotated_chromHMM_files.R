library(yaml)
library(ggthemes)
library(tidyr)
require(scales)
#library(GenomicRanges)
library(bedr)
library(reshape2)

radio = read_yaml("/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/radiofile.yml")


chrhmmHomer = read.table(radio$chromhmm_merged_annotated,
                         sep="\t", header=T, quote="", stringsAsFactors = FALSE)

cleanup_chromhmm_homer <- function(chrhmmHomer){
  names(chrhmmHomer)[1] <- "name"
  head(chrhmmHomer)
  
  # recode the annotation to remove mention of the transcript
  chrhmmHomer <-chrhmmHomer[!is.na(chrhmmHomer$Annotation),]
  chrhmmHomer$Annotation[grep("promoter", chrhmmHomer$Annotation)] <- "Promoter-TSS"
  chrhmmHomer$Annotation[grep("TTS", chrhmmHomer$Annotation)] <- "TTS"
  chrhmmHomer$Annotation[grep("intron", chrhmmHomer$Annotation)] <- "Intron"
  chrhmmHomer$Annotation[grep("exon", chrhmmHomer$Annotation)] <- "Exon"
  chrhmmHomer$name <- gsub("-[0-9]*","",chrhmmHomer$name, perl = T)
  chrhmmHomer <- chrhmmHomer[-grep(pattern="GL*",chrhmmHomer$Chr),]
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "Y"),]
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "X"),]
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "MT"),]
  
  
  chrhmmHomer$Annotation <- factor(chrhmmHomer$Annotation, levels =rev(c('Exon',
                                                                         'Intron',
                                                                         'Promoter-TSS',
                                                                         'TTS','Intergenic')))
  chrhmmHomer$name <- factor(chrhmmHomer$name, levels = c("E1", "E2", "E3", "E4", "E5",
                                                          "E6", "E7", "E8", "E9", "E10", "E11"))
  head(chrhmmHomer)
  table(chrhmmHomer$name)
  
  chrhmmHomer$state <- chrhmmHomer$name
  chrhmmHomer$name <- gsub("E","",chrhmmHomer$name, perl = T)
  
  chrhmmHomer$name <- dplyr::recode(as.numeric(chrhmmHomer$name),
                                         "Heterochrom_lowsignal","Polycomb_repressed",
                                         "Heterochrom_lowsignal","Poised_promoter",
                                         "Weak_promoter","Active_promoter","Strong_enhancer",
                                         "Weak_promoter","Poised_enhancer","Strong_enhancer",
                                         "Heterochrom_lowsignal")
  head(chrhmmHomer)
  table(chrhmmHomer$name)
  return(chrhmmHomer)
}
chrhmmHomer_towrite <- cleanup_chromhmm_homer(chrhmmHomer)
write.table(chrhmmHomer_towrite,radio$chromhmm_homer,quote=F,sep = "\t", col.names = T, row.names = F)

chrhmmHomer_towrite_clean_bed <- chrhmmHomer_towrite[,c(2,3,4,1)]
tempname = tempfile(pattern = "file", fileext = ".bed")

write.table(chrhmmHomer_towrite_clean_bed,tempname,quote=F,sep = "\t", col.names = F, row.names = F)
system2('sort', args = c('-k1,1', '-k2,2n', tempname),
        stdout = radio$chromhmm.c, stderr = 'stderr.txt')


chrhmmHomer = read.table(radio$chromhmm_merged_annotated_tss_ccds,
                         sep="\t", header=T, quote="", stringsAsFactors = FALSE, fill = TRUE)
chrhmmHomer_towrite <- cleanup_chromhmm_homer(chrhmmHomer)
write.table(chrhmmHomer_towrite,radio$chromhmm_homer_tss_ccds,quote=F,sep = "\t", col.names = T, row.names = F)
################################
# For an ENCODE sample annotated by HOMER
################################

recode18 <- function(chrhmmHomer){
  names(chrhmmHomer)[1] <- "name"
  chrhmmHomer <-chrhmmHomer[!is.na(chrhmmHomer$Annotation),]
  chrhmmHomer$Annotation[grep("promoter", chrhmmHomer$Annotation)] <- "Promoter-TSS"
  chrhmmHomer$Annotation[grep("TTS", chrhmmHomer$Annotation)] <- "TTS"
  chrhmmHomer$Annotation[grep("intron", chrhmmHomer$Annotation)] <- "Intron"
  chrhmmHomer$Annotation[grep("exon", chrhmmHomer$Annotation)] <- "Exon"
  chrhmmHomer$name <- gsub("-[0-9]*","",chrhmmHomer$name, perl = T)
  chrhmmHomer$name <- as.numeric(chrhmmHomer$name)
  chrhmmHomer$name <- dplyr::recode(chrhmmHomer$name,"TssA","TssAFlnk","TssAFlnk","TssAFlnk","Tx","Tx","EnhG","EnhG","EnhA",
                                    "EnhA","EnhWk","ZNF_Rpts","Het","TssBiv","EnhBiv","ReprPC",
                                    "ReprPC","Quies")
  chrhmmHomer$name <- factor(chrhmmHomer$name,levels = rev(c("TssA","TssAFlnk","Tx","EnhG","EnhA",
                                                             "EnhWk","ZNF_Rpts","Het","TssBiv","EnhBiv",
                                                             "ReprPC","Quies")))
  chrhmmHomer$Annotation <- factor(chrhmmHomer$Annotation, levels =rev(c('Exon',
                                                                         'Intron',
                                                                         'Promoter-TSS',
                                                                         'TTS','Intergenic')))
  if(length(which(chrhmmHomer$Chr == "M")) != 0){
    chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "M"),]
  }
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "X"),]
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "Y"),]
  return(chrhmmHomer)
}

recode15 <- function(chrhmmHomer){
  names(chrhmmHomer)[1] <- "name"
  chrhmmHomer <-chrhmmHomer[!is.na(chrhmmHomer$Annotation),]
  chrhmmHomer$Annotation[grep("promoter", chrhmmHomer$Annotation)] <- "Promoter-TSS"
  chrhmmHomer$Annotation[grep("TTS", chrhmmHomer$Annotation)] <- "TTS"
  chrhmmHomer$Annotation[grep("intron", chrhmmHomer$Annotation)] <- "Intron"
  chrhmmHomer$Annotation[grep("exon", chrhmmHomer$Annotation)] <- "Exon"
  chrhmmHomer$name <- gsub("-[0-9]*","",chrhmmHomer$name, perl = T)
  chrhmmHomer$name <- as.numeric(chrhmmHomer$name)
  chrhmmHomer$name <- dplyr::recode(chrhmmHomer$name,"TssA","TssAFlnk","Tx","Tx","Tx","EnhG","EnhA",
                                    "ZNF_Rpts","Het","TssBiv","BivFlnk","EnhBiv","ReprPC",
                                    "ReprPC","Quies")
  chrhmmHomer$name <- factor(chrhmmHomer$name,levels = rev(c("TssA","TssAFlnk","Tx","EnhG","EnhA",
                                                             "ZNF_Rpts","Het","TssBiv","BivFlnk","EnhBiv",
                                                             "ReprPC","Quies")))
  chrhmmHomer$Annotation <- factor(chrhmmHomer$Annotation, levels =rev(c('Exon',
                                                                         'Intron',
                                                                         'Promoter-TSS',
                                                                         'TTS','Intergenic')))
  if(length(which(chrhmmHomer$Chr == "M")) != 0){
    chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "M"),]
  }
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "X"),]
  chrhmmHomer = chrhmmHomer[-which(chrhmmHomer$Chr == "Y"),]
  return(chrhmmHomer)
}


###Right ventricule from ENCODE

HRT.VENT.R=read.table(radio$ENCODE.HRT.VENT.R.a,
                      sep="\t", header=T, quote="", stringsAsFactors = FALSE, fill=TRUE)
HRT.VENT.R.a <- recode18(HRT.VENT.R)
# head(HRT.VENT.R.a)
write.table(HRT.VENT.R.a,radio$HRT.VENT.R.a,quote=F,sep = "\t", col.names = T, row.names = F)

###FETAL heart from ENCODE

HRT.FET.a=read.table(radio$ENCODE.HRT.FET.a,
                   sep="\t", header=T, quote="", stringsAsFactors = FALSE, fill=TRUE)
HRT.FET.a <- recode15(HRT.FET.a)
head(HRT.FET.a)
write.table(HRT.FET.a,radio$HRT.FET.a, quote=F,sep = "\t", col.names = T, row.names = F)


###IPS heart from ENCODE

ENCODE.IPS = read.table(radio$ENCODE.IPS.a,sep="\t", header=T, quote="", stringsAsFactors = FALSE, fill = TRUE)

ENCODE.IPS <- recode18(ENCODE.IPS)
head(ENCODE.IPS)
write.table(ENCODE.IPS,radio$IPS.a, quote=F,sep = "\t", col.names = T, row.names = F)


# Read in ENCODE 18 states data
# String as factor =F very important here
# Don't know why but the ENCODE chromHMM have the header of the bed as the last line of the file
HRT.VENT.L = read.table(radio$ENCODE.HRT.VENT.L,stringsAsFactors=F)
HRT.VENT.R = read.table(radio$ENCODE.HRT.VENT.R,stringsAsFactors=F)
HRT.ATR.R = read.table(radio$ENCODE.HRT.ATR.R,stringsAsFactors=F)
HRT.FET = read.table(radio$ENCODE.HRT.FET,stringsAsFactors=F)

#IPS
IPS = read.table(radio$ENCODE.IPS,stringsAsFactors=F)
IPS= IPS[-length(IPS$V1),]


#Fetal male brain
BRN.FET = read.table(radio$ENCODE.BRN.FET,stringsAsFactors=F)
BRN.FET = BRN.FET[-length(BRN.FET$V1),]

#Male brain
BRN.ADLT = read.table(radio$ENCODE.BRN.ADLT,stringsAsFactors=F)
BRN.ADLT= BRN.ADLT[-length(BRN.ADLT$V1),]

# Recode values with human readable categories
recode18 <-function(chromHMM18){
  chromHMM18$V2 = as.integer(chromHMM18$V2)
  chromHMM18$V3 = as.integer(chromHMM18$V3)
  chromHMM18$V4 <- as.numeric(chromHMM18$V4)
  chromHMM18$V4 = dplyr::recode(chromHMM18$V4,"TssA","TssAFlnk","TssAFlnk","TssAFlnk","Tx","Tx","EnhG","EnhG","EnhA",
                                "EnhA","EnhWk","ZNF_Rpts","Het","TssBiv","EnhBiv","ReprPC",
                                "ReprPC","Quies")
  if(length(which(chromHMM18$V1 == "chrM")) != 0){
    chromHMM18 = chromHMM18[-which(chromHMM18$V1 == "M"),]
  }
  chromHMM18 = chromHMM18[-which(chromHMM18$V1 == "X"),]
  chromHMM18 = chromHMM18[-which(chromHMM18$V1 == "Y"),]
  chromHMM18$V1 <- gsub("chr","",chromHMM18$V1, perl = T)
  
  return(chromHMM18)
}

# Recode values with human readable categories
recode15 <-function(chromHMM15){
  chromHMM15$V2 = as.integer(chromHMM15$V2)
  chromHMM15$V3 = as.integer(chromHMM15$V3)
  chromHMM15$V4 <- as.numeric(chromHMM15$V4)
  chromHMM15$V4 = dplyr::recode(chromHMM15$V4,"TssA","TssAFlnk","Tx","Tx","Tx","EnhG","EnhA",
                                "ZNF_Rpts","Het","TssBiv","BivFlnk","EnhBiv","ReprPC",
                                "ReprPC","Quies")
  if(length(which(chromHMM15$V1 == "chrM")) != 0){
    
    chromHMM15 = chromHMM15[-which(chromHMM15$V1 == "M"),]
  }
  chromHMM15 = chromHMM15[-which(chromHMM15$V1 == "X"),]
  chromHMM15 = chromHMM15[-which(chromHMM15$V1 == "Y"),]
  chromHMM15$V1 <- gsub("chr","",chromHMM15$V1, perl = T)
  
  return(chromHMM15)
}

HRT.VENT.L = recode18(HRT.VENT.L)
head(HRT.VENT.L)
write.table(HRT.VENT.L,radio$HRT.VENT.L.c, quote=F,sep = "\t",col.names = F, row.names = F)

HRT.VENT.R = recode18(HRT.VENT.R)
write.table(HRT.VENT.R,radio$HRT.VENT.R.c, quote=F,sep = "\t", col.names = F,row.names = F)

HRT.ATR.R = recode18(HRT.ATR.R)
write.table(HRT.ATR.R,radio$HRT.ATR.R.c, quote=F,sep = "\t", col.names = F,row.names = F)

HRT.FET = recode15(HRT.FET)
write.table(HRT.FET,radio$HRT.FET.c, quote=F,sep = "\t",col.names = F, row.names = F)

IPS = recode18(IPS)
write.table(IPS,radio$IPS.c, quote=F,sep = "\t",col.names = F, row.names = F)

BRN.ADLT = recode18(BRN.ADLT)
write.table(BRN.ADLT,radio$BRN.ADLT.c, quote=F,sep = "\t", col.names = F,row.names = F)

BRN.FET = recode15(BRN.FET)
write.table(BRN.FET,radio$BRN.FET.c, quote=F,sep = "\t", col.names = F,row.names = F)


