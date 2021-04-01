#####################task1 boxplot of E9.5 scRNA marker gene expression
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("/scratch/hz9fq/toronto_project/E95_scRNA/")
E95<-readRDS("E95.rds")

mtx<-E95@assays$RNA@data

#read in organ annotation
meta <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95_scRNA/matrics/GSE136689_MetaFile_DE_and_SM_Cells_E8.5_E9.0_E9.5.txt",header=T,row.names=1,stringsAsFactors=F,sep='\t')

meta.use <- meta[which(meta$Type=="Endoderm" & (meta$Stages=="E9.5Anterior" | meta$Stages=="E9.5Posterior")),]

organs=c("intestine","pancreas","stomach","lung","liver","pharynx","esophagus")
pal=c("#D51F26","#D3C035","#C06CAB","#89C75F","#90D5E4","#E6C2DC","#F47D2B")

meta.use$celltype <- meta.use$LineageAnnotations
meta.use$celltype[meta.use$celltype=="duodenum"] <- "intestine"
meta.use$celltype[meta.use$celltype=="dorsal pancreas"] <- "pancreas"
meta.use$celltype[meta.use$celltype=="anaterior stomach" | meta.use$celltype=="posterior stomach"] <- "stomach"
meta.use$celltype[meta.use$celltype=="respiratory-lung"] <- "lung"
meta.use$celltype[meta.use$celltype=="hepatoblast" | meta.use$celltype=="liver-hepatocyte"] <- "liver"
meta.use$celltype[meta.use$celltype=="pharynx"] <- "pharynx"
meta.use$celltype[meta.use$celltype=="esophagus-1" | meta.use$celltype=="esophagus-2"] <- "esophagus"

meta.use <- meta.use[which(meta.use$celltype%in%organs),]

geneList = c("Hhex","Creb3l3","Pdx1","Neurog3","Cdx2","Tff3","Osr1","Sox2","Nkx2-1","Nkx6-1","Ptf1a")
#geneList = c("Mafb", "Hnf4g", "Grhl3", "Hopx")

pdf("result_plots/E95scRNA_markerGene_boxplot.pdf")
for(K in geneList){
    plotdat <- data.frame(exp=mtx[K,rownames(meta.use)],celltype=meta.use$celltype)
    plotdat$celltype <- factor(plotdat$celltype,levels=organs)

    boxplot(exp~celltype,plotdat,outline=T,cex.axis=0.8,col=pal,main=K)
}
dev.off()
