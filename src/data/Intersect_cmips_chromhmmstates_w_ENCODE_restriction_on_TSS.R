library(yaml)
library(ggthemes)
library(tidyr)
require(scales)
#library(GenomicRanges)
#library(bedr)
library(reshape2)

radio = read_yaml("/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/radiofile.yml")



####fonction
get_inter<-function(sample,encode,namesample,nameencode,outdircompa){
  dir.create(outdircompa,recursive=TRUE)

  mat = matrix(NA,length(table(sample$name)),length(levels(encode$V4)))
  
  rownames(mat) <- names(table(sample$name))
  colnames(mat) <- levels(encode$V4)
  
  samp1 = sample[,c(2,3,4)]
  samp1.valid =  is.valid.region(samp1, verbose = F,check.chr = FALSE)
  tempname = tempfile(pattern = "file", fileext = ".bed")
  write.table(sample[samp1.valid,c(2,3,4,1)], file = tempname, sep = "\t", col.names = F, quote = F, row.names = F)
  system(paste0("sort -k1,1 -k2,2n ", tempname, " > ", outdir,"/",namesample,"_sorted.bed"))
  
  for( i in rownames(mat)){
    
    #Get the region with the specific annotation
    samp1 = sample[which(sample$name == i),c(2,3,4)]
    
    samp1.valid =  is.valid.region(samp1, verbose = F,check.chr = FALSE)
    tempname = tempfile(pattern = "file", fileext = ".bed")
    write.table(samp1[samp1.valid,], file = tempname, sep = "\t", col.names = F, quote = F, row.names = F)
    system(paste0("sort -k1,1 -k2,2n ", tempname, " > ", outdircompa,"/",i,"_sorted.bed"))
  }
  samp1 = encode[,c(1,2,3)]
  
  samp1.valid =  is.valid.region(samp1, verbose = F,check.chr = FALSE)
  head(samp1.valid)
  tempname = tempfile(pattern = "file", fileext = ".bed")
  write.table(encode[samp1.valid,], file = tempname, sep = "\t", col.names = F, quote = F, row.names = F)
  system(paste0("sort -k1,1 -k2,2n ", tempname, " > ",outdircompa,"/", nameencode,"_sorted.bed"))
  
  for( i in colnames(mat)){
    
    samp1 = encode[which(encode$V4 == i),c(1,2,3)]
    
    samp1.valid =  is.valid.region(samp1, verbose = F,check.chr = FALSE)
    head(samp1.valid)
    tempname = tempfile(pattern = "file", fileext = ".bed")
    write.table(samp1[samp1.valid,], file = tempname, sep = "\t", col.names = F, quote = F, row.names = F)
    system(paste0("sort -k1,1 -k2,2n ", tempname, " > ",outdircompa,"/", i,"_sorted.bed"))
  }
  
  
  #FULL intersect
  system(paste0("bedtools intersect -a ",outdircompa,"/", namesample,"_sorted.bed", " -b ",outdircompa,"/", nameencode,"_sorted.bed"," -wa -wb > ",outdircompa,"/",namesample,"_",nameencode,"_full_intersect.bed" ))     
  
  
  for( j in rownames(mat)){
    #Get the sample
    system(paste0("bedtools intersect -a ",outdircompa,"/", j,"_sorted.bed", " -b ",radio$HRT.VENT.R.c," -wa -wb > ",outdircompa,"/",j,"_inter_all_encode.bed" ))     
    
    for (i in colnames(mat)){
      
      #print(paste0("bedtools jaccard -a ", all_samps[i],"_sorted.bed", " -b ", all_samps[j],"_sorted.bed ", "> jaccard.temp" ))
      # Not in A
      system(paste0("bedtools intersect -a ",outdircompa,"/", j,"_sorted.bed", " -b ",outdircompa,"/", i,"_sorted.bed -v > ",outdircompa,"/","In_",j,"_but_not_in_",i,".bed" ))     
      system(paste0("bedtools intersect -a ",outdircompa,"/", i,"_sorted.bed", " -b ",outdircompa,"/", j,"_sorted.bed -v > ",outdircompa,"/","In_",i,"_but_not_in_",j,".bed" ))     
      system(paste0("bedtools intersect -a ",outdircompa,"/", j,"_sorted.bed", " -b ",outdircompa,"/", i,"_sorted.bed -wa -u > ",outdircompa,"/",j,"_inter_",i,".bed" ))
      
      
      
    }
  }         
  return(NULL)
}



cmips = read.table(radio$chromhmm_homer_tss_ccds,sep="\t", header=T, quote="")
head(cmips)

encodevent = read.table(radio$HRT.VENT.R.c, quote="")
encodefet = read.table(radio$HRT.FET.c, quote="")

##For ventricul Right
outdir = "/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/processed/ventricule-right-comparison"
get_inter(cmips,encodevent,"CM-IPS","VENT.R",outdir)
#For fetal

outdir = "/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/processed/fetal-heart-comparison"
get_inter(cmips,encodefet,"CM-IPS","FET.HRT",outdir)

####################only 1kb windows around tss from emsembl_havana CCDS
cmips_tss_1kb = cmips[abs(cmips$Distance.to.TSS)<501,]
outdir = "/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/processed/ventricule-right-comparison-only_TSS_1kb_window/"
get_inter(cmips_tss_1kb,encodevent,"CM-IPS","VENT.R",outdir)

outdir = "/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/processed/fetal-heart-comparison-only_TSS_1kb_window/"
get_inter(cmips_tss_1kb,encodefet,"CM-IPS","FET.HRT",outdir)

##Only states that overlap a TSS from ensembl_havana CCDS

overlap=(cmips$End- cmips$Start)/2 > abs(cmips$Distance.to.TSS)
cmips_overlap = cmips[overlap,]

outdir = "/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/processed/ventricule-right-comparison-TSS-overlap"
get_inter(cmips_overlap,encodevent,"CM-IPS","VENT.R",outdir)

outdir = "/home/itx1433/Documents/PROJ/CM-IPS_maturity/data/processed/fetal-heart-comparison-TSS-overlap"
get_inter(cmips_overlap,encodefet,"CM-IPS","FET.HRT",outdir)


