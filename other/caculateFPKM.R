#get gene exons length
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("mm10_refgenes.gtf",format="gtf")

exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

exons_gene_lens.V <- c()
for(i in 1:length(exons_gene_lens)){
    exons_gene_lens.V<-c(exons_gene_lens.V,exons_gene_lens[[i]])
}

out.dat <- data.frame(gene=names(exons_gene_lens),len=exons_gene_lens.V)

write.table(out.dat,"mm10_refgenes_exon_gene_length.txt",col.names=F,row.names=F,quote=F,sep='\t')

#read in count matrix
library(Seurat)
library(dplyr)
library(Matrix)

setwd("/scratch/hz9fq/toronto_project/E95_scRNA")
E95 <- readRDS("E95.rds")

count.mtx <- E95@assays$RNA@counts
ct.v <- as.character(E95@meta.data$celltype)
names(ct.v) <- colnames(count.mtx)

libsize <- apply(count.mtx,2,sum)

lens.tab <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_refgenes_exon_gene_length.txt",header=F)

rownames(lens.tab) <- lens.tab$V1
gene_use <- rownames(count.mtx)[rownames(count.mtx)%in%rownames(lens.tab)]

count.mtx <- count.mtx[gene_use,]
lens.v <- lens.tab[rownames(count.mtx),2]

#caculate FPKM
FPKM.mtx <- count.mtx
for(i in 1:nrow(FPKM.mtx)){
    FPKM.mtx[i,] <- FPKM.mtx[i,]/lens.v[i]
}

for(i in 1:ncol(FPKM.mtx)){
    FPKM.mtx[,i] <- FPKM.mtx[,i]/libsize[i]
}

FPKM.mtx <- FPKM.mtx*1e9

write.table(FPKM.mtx,"E95_FPKM_matrix.txt",row.names=T,col.names=T,sep='\t',quote=F)

FPKM_value <- c()
for(i in names(table(ct.v))){
    cells <- names(ct.v)[which(ct.v==i)]
    print(length(cells))
    tmp.mtx <- FPKM.mtx[,cells]
    FPKM_value <- c(FPKM_value,apply(tmp.mtx,1,mean))
}

final_mtx <- matrix(FPKM_value,ncol=length(table(ct.v)))

colnames(final_mtx) <- names(table(ct.v))
rownames(final_mtx) <- rownames(FPKM.mtx)

write.table(final_mtx,"E95_FPKM_matrix_byOrgan.txt",row.names=T,col.names=T,sep='\t',quote=F)

organs = c("intestine","pancreas","stomach","lung")
ct_byOrgan = list(c("duodenum"),c("dorsal pancreas"),c("anaterior stomach","posterior stomach"),c("respiratory-lung"))

FPKM_value <- c()
for(i in 1:4){
    cells <- names(ct.v)[which(ct.v%in%ct_byOrgan[[i]])]
    print(length(cells))
    tmp.mtx <- FPKM.mtx[,cells]
    FPKM_value <- c(FPKM_value,apply(tmp.mtx,1,mean))
}

final_mtx <- matrix(FPKM_value,ncol=4)

colnames(final_mtx) <- organs
rownames(final_mtx) <- rownames(FPKM.mtx)

write.table(final_mtx,"E95_FPKM_matrix_4Organs.txt",row.names=T,col.names=T,sep='\t',quote=F)

plotV <- final_mtx[,4]
plotV <- plotV[which(plotV<=100)]
hist(plotV,breaks=seq(0,100,length.out=101),freq=F)
