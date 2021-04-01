#for gene, use signature genes(hg38) called from TCGA_sigGene_GSEA
#for ATAC-seq peaks, use full DESeq2 bed for bart
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/DESeq2_full_table")
res<-read.table("STAD_peaks_DESeq2.txt",header=T)
options(scipen = 200)

coord <- sapply(rownames(res),function(x){
    tmp1 <- strsplit(x,":")
    y <- tmp1[[1]][2]
    tmp2 <- strsplit(y,"-")
    return(c(tmp1[[1]][1],tmp2[[1]][1],tmp2[[1]][2]))
})

output <- as.data.frame(matrix(as.vector(unlist(coord)),ncol=3,byrow=T))
output$V2 <- as.numeric(output$V2)
output$V3 <- as.numeric(output$V3)
output <- cbind(output,".",res$log2FoldChange,".")
output[which(is.na(output[,5])==TRUE),5] <- 0
write.table(output,"STAD_peakswithFC_bed6.bed",row.names=F,col.names=F,quote=F,sep="\t")

#reverse
output <- as.data.frame(matrix(as.vector(unlist(coord)),ncol=3,byrow=T))
output$V2 <- as.numeric(output$V2)
output$V3 <- as.numeric(output$V3)
output <- cbind(output,".",-res$log2FoldChange,".")
output[which(is.na(output[,5])==TRUE),5] <- 0
write.table(output,"STAD_peakswithFC_reverse_bed6.bed",row.names=F,col.names=F,quote=F,sep="\t")
