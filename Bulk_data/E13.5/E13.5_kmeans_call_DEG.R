set.seed(1)
t0<-read.table("FPKM_collection.txt",header=T,stringsAsFactors=F)
t1<-t0[,7:14]
rownames(t1)<-t0$Gene.ID

t1 <- data.frame(apply(t1[,1:2],1,mean),apply(t1[,3:4],1,mean),apply(t1[,5:6],1,mean),apply(t1[,7:8],1,mean))
colnames(t1) <- c("E135_int_RNA","E135_lung_RNA","E135_pancreas_RNA","E135_stomach_RNA")

t1$var <- apply(t1,1,var)
t1$mean <- apply(t1,1,mean)
t1$vm <- t1$var/t1$mean
t1$vm[which(is.na(t1$vm))] <- 0
t1<-t1[order(t1$vm,decreasing=T),]

t1$max_FPKM <- apply(t1[,1:4],1,max)

#2500,14
#t1top<-t1[1:4000,1:4]
t1top <- t1[which(t1$max_FPKM>1),1:4]
t1top<-t1[1:8000,1:4]
t1top <- log2(t1top+1)

#quantile normalization
mtx <- as.matrix(t1top)
mtx <- preprocessCore::normalize.quantiles(mtx)
rownames(mtx) <- rownames(t1top)
colnames(mtx) <- colnames(t1top)
t1top <- as.data.frame(mtx)

kc <- kmeans(t1top,14,iter.max=20,nstart=25)
t1top$cluster<-kc$cluster
t1top <- t1top[order(t1top$cluster),]

#t1top.use <- t1top[which(t1top$cluster %in% c(1,2,6,8,13)),]
t1top.use<-t1top
plotdat <- as.matrix(t1top.use[,1:4])
print(table(t1top.use$cluster))
gaps<-cumsum(table(t1top.use$cluster))
gaps<-gaps[1:(length(gaps)-1)]

anno_row<-data.frame(tag=t1top.use$cluster)
rownames(anno_row)<-rownames(t1top.use)
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
anno_color<-list(tag=stallion_pal[1:length(table(t1top.use$cluster))])
names(anno_color[[1]])<-c(1:length(table(t1top.use$cluster)))

library(pheatmap)
pdf("heatmap.pdf")
pheatmap(plotdat,breaks=seq(min(plotdat),quantile(plotdat,0.95),length.out=101),cluster_cols=F,cluster_rows=F,show_rownames=F,gaps_row=gaps,
annotation_row=anno_row,annotation_colors=anno_color)
dev.off()

saveRDS(kc,"kmeansResults_top8000exp3.rds")

#write out DEG
#read in gene range file
setwd("../DEG_kmeans/")

promoter <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_TSS3kb.bed",header=F,stringsAsFactors=F)
rownames(promoter) <- promoter$V5
genebody <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_genebody50kb.bed",header=F,stringsAsFactors=F)
rownames(genebody) <- genebody$V5

promoter.out<-promoter[rownames(t1top)[which(t1top$cluster==3)],]
genebody.out<-genebody[rownames(t1top)[which(t1top$cluster==3)],]
print(nrow(promoter.out))
genelist<-data.frame(V1=rownames(promoter.out))

write.table(genelist,"intestine_DEG_genelist.txt",col.names=F,row.names=F,sep='\t',quote=F)
write.table(promoter.out,"intestine_DEG_promoter.bed",col.names=F,row.names=F,sep='\t',quote=F)
write.table(genebody.out,"intestine_DEG_genebody.bed",col.names=F,row.names=F,sep='\t',quote=F)

