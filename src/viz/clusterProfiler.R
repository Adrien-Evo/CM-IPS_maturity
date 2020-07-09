
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(DOSE))
suppressMessages(library(pathview))
suppressMessages(library(enrichplot))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library(KEGG.db))
suppressMessages(library(magrittr))
suppressMessages(library(msigdbr))
suppressMessages(library(rWikiPathways))
library(tidyr)
library(yaml)
library(biomaRt)




rootfolder = "/ceph-recherche/shares/u1087/afoucal/CM-IPS_maturity/"
radio = read_yaml(file.path(rootfolder,"radiofile.yml"))

## Read the input bed file
kmeans = read.table(file.path(rootfolder,radio$kmeans.all), h=T, sep ="\t", stringsAsFactors=FALSE)

gmt = file.path(rootfolder,radio$wikipathway.gmt)


summarize_gsea <- function(gsea_result, database){

if(length(gsea_result) == 0){
return(c(database,"GeneSetEnrichment",deparse(substitute(gsea_result)),0,0.05,"BH"))
} else {
pval = gsea_result@params$pvalueCutoff
p_adjust_method = gsea_result@params$pAdjustMethod

nb_enriched_term = sum(gsea_result@result$p.adjust < pval)

return(c(database,"GeneSetEnrichment",deparse(substitute(gsea_result)),nb_enriched_term,pval,p_adjust_method))
}
}

summarize_ora <- function(ora_result, database){
if(length(ora_result) == 0){
return(c(database,"OverRepresentation",deparse(substitute(ora_result)),0,0.05,"BH"))
} else {
pval = ora_result@pvalueCutoff
p_adjust_method = ora_result@pAdjustMethod

nb_enriched_term = sum(ora_result@result$p.adjust < pval)

return(c(database,"OverRepresentation",deparse(substitute(ora_result)),nb_enriched_term,pval,p_adjust_method))
}
}

###################################################
######## Select a cluster
###################################################

colourcmips= read.table(file.path(rootfolder,radio$color.mapping.cmips),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
##Keeping only annotation that we have
colourcmips = colourcmips[c(1,2,3,4,5,10,11),]

for(state in colourcmips$V1){
  #state= "Active_promoter"
plotfolder = file.path(rootfolder,paste0("plots/clustering/clusterProfiler/",state))

## Here cluster 3 from kmeans 10
select = kmeans[which(kmeans[,1] == state),]
dir.create(plotfolder)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

for(k in 1:10){


total_k = 10
############
############
ENST = rownames(select[which(select$singular_kmeans == k),])

#ENST = rownames(kmeans[which(kmeans$kmeans15 == k),])

### Looking at cluster 1 in details 

#This will need the geneId slot for ENTREZ gene ID

###############################
#####BiomaRt conversion
###############################
results <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"), filters = "ensembl_transcript_id",
           values = ENST, mart = mart)


#ids <- bitr(ENST, fromType="ENSEMBLTRANS", toType=c("ENTREZID","GENENAME","SYMBOL"), OrgDb="org.Hs.eg.db")

#gene = ids$ENTREZID
gene= results$entrezgene_id
fold_flag = FALSE

enrich_function(k,plotfolder,total_k,gene)
}
}

enrich_function <- function(k, plotfolder,total_k, gene){
# WIKIPATHWAY

# SUMMARIZED
SUMMARIZED = vector()

cat("WikiPathways\n")

#Downloading and reading the latest human gmt archive using the rWikiPathways package
wp2gene <- read.gmt(gmt)

wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

#Over Representation Analysis with only gene name
ora_wp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

ora_wp <- setReadable(ora_wp, org.Hs.eg.db, keyType = "ENTREZID")

# --- Summarizing outputs -----

SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_wp,"WikiPathways") )



# CELL MARKER

cell_markers <- read.table('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt',h=T,sep = "\t", stringsAsFactors = FALSE,quote="\"") %>%
  dplyr::select(tissueType, cancerType, cellName, geneID) %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep="-") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
  summary(unlist(lapply(cell_markers$cellMarker,nchar)))

ora_cm <- enricher(gene,
  TERM2GENE=cell_markers,
  minGSSize=1,
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_cm,"cellMarker") )

# MSIGdb

# Molecular signature from http://software.broadinstitute.org/gsea/msigdb

### Here there is a lot of categories to be used. Here is used the curated gene sets C2
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)


ora_ms <- enricher(gene,
  TERM2GENE=m_t2g,
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05)
# --- Summarizing outputs -----
SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_ms,"MSigDb_C2") )


# GO

ora_CC_go <- enrichGO(gene = gene,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable     = TRUE)
SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_CC_go,"GO_CellularComponent") )

ora_BP_go <- enrichGO(gene = gene,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable     = TRUE)
SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_BP_go,"GO_BiologicalProcess") )

ora_MF_go <- enrichGO(gene = gene,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable     = TRUE)
SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_MF_go,"GO_MolecularFonction") )


#KEGG 
search_kegg_organism('hsa', by='kegg_code')

# Here you need to say human because theres a function called organismMapper that will translate it (horribly) in the proper organism
ora_kegg <- enrichKEGG(gene = gene,organism     = 'hsa',pvalueCutoff = 0.05)

SUMMARIZED = rbind(SUMMARIZED, summarize_ora(ora_kegg,"KEGG") )



# SUMARRIZED

summarized.df = data.frame(DB = SUMMARIZED[,1],test=SUMMARIZED[,2],varName=SUMMARIZED[,3],enrichedTerms=SUMMARIZED[,4],pVal=SUMMARIZED[,5],pAdjustMethod=SUMMARIZED[,6],stringsAsFactors = FALSE)

write.table(summarized.df,file.path(rootfolder,"data/processed/clustering/clusterProfiler",paste0("kmeans_",total_k,"_clusterNB_",k,".clusterProfiler_summarized.txt")),quote = FALSE, sep = "\t",row.names = FALSE)

# VIZ

#  Visualization for Over Representation test

viz_ORA <- function(enrichR_output, database, GO = FALSE){

  # --------- Dotplot ----------
  dotplot_plot = dotplot(enrichR_output, showCategory=30) + ggtitle(database)
  try(ggsave(paste0("kmeans_",total_k,"_clusterNB_",k,"_",database,"_dotplot_ora.png"),dotplot_plot,width = 10, height = 6))
  
  # ----- Enrichment map or goplot ------
  if(GO == FALSE){
    emapplot_plot = emapplot(enrichR_output) + ggtitle(database)
    try(ggsave(paste0("kmeans_",total_k,"_clusterNB_",k,"_",database,"_enrichment_map_ora.png"),emapplot_plot,width = 16, height = 16))
    
  }else{
    goplot_plot = goplot(enrichR_output) + ggtitle(database)
    try(ggsave(paste0("kmeans_",total_k,"_clusterNB_",k,"_",database,"_enrichment_map_ora_GO.png"),goplot_plot,width = 16, height = 16))
  }

  # --------- Heatmap ---------- 
  heatplot_plot = heatplot(enrichR_output) + ggtitle(database)
  try(ggsave(paste0("kmeans_",total_k,"_clusterNB_",k,"_",database,"_heatmap_ora.png"),heatplot_plot,width = 24, height = 12))

  # ---- Network w/t genes ----- 

  cnetplot_plot = cnetplot(enrichR_output) + ggtitle(database)
  try(ggsave(paste0("kmeans_",total_k,"_clusterNB_",k,"_",database,"_cnetplot_ora.png"),cnetplot_plot,width = 12, height = 12))
}


# ~~~~~~~~ Move to plots directory ~~~~~~~~~~
setwd(plotfolder)
for( i in 1:length(summarized.df[,1])){

anaOut = get(summarized.df$varName[i])
anaType = summarized.df$varName[i]
GO = FALSE
print(summarized.df[i,])
cat("\n\n")
#  Checking if the analyses is using GO 
if(length(grep("go",anaType) == 1)){
GO = TRUE
}
if(as.numeric(summarized.df$enrichedTerms[i]) > 1){
if(length(grep("ora",anaType)) == 1){viz_ORA(enrichR_output = anaOut, database = summarized.df$DB[i], GO = GO)}
}
}

}

# hsa05414 <- pathview(gene.data  = geneList,
#                      pathway.id = "hsa05414",
#                      species    = "hsa",
#                      limit      = list(gene=max(abs(geneList)), cpd=1))

# hsa04310 <- pathview(gene.data  = geneList,
#                      pathway.id = "hsa04310",
#                      species    = "hsa",
#                      limit      = list(gene = 4, cpd = 4))


# data(gse16873.d)
# v.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
#                   species = "hsa", out.suffix = "gse16873")


