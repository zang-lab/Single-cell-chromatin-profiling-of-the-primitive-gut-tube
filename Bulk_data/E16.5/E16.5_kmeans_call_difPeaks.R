##########################task3 E16.5 dif peaks
filenames=c("E16_FS_ATAC.peak_1bin.bed","E16_HS_ATAC.peak_1bin.bed","E16_SI_ATAC.peak_1bin.bed","E16_CL_ATAC.peak_1bin.bed")
signal=c()
for(fn in filenames){
    tmp.t <- read.table(fn,header=F,stringsAsFactors=F)
    signal<-c(signal,tmp.t$V4)
}
signal <- as.data.frame(matrix(signal,ncol=4))

rownames(signal) <- paste0(tmp.t$V1,":",tmp.t$V2,"-",tmp.t$V3)
colnames(signal) <- c("E16_FS","E16_HS","E16_SI","E16_CL")

signal$var <- apply(signal,1,var)
signal$mean <- apply(signal,1,mean)
signal$vm <- signal$var/signal$mean
signal$vm[which(is.na(signal$vm))] <- 0
signal<-signal[order(signal$vm,decreasing=T),]

#quantile normalization
mtx <- as.matrix(signal[,1:4])
mtx <- preprocessCore::normalize.quantiles(mtx)
rownames(mtx) <- rownames(signal)
colnames(mtx) <- colnames(signal)[1:4]

signaltop_s1 <- as.data.frame(mtx)
signaltop_s1 <- log2(signaltop_s1+1)

kc <- kmeans(signaltop_s1,13,iter.max=200,nstart=25,algorithm="MacQueen")
signaltop_s2 <- signaltop_s1
signaltop_s2$cluster<-kc$cluster
signaltop_s2 <- signaltop_s2[order(signaltop_s2$cluster),]

plotdat <- as.matrix(signaltop_s2[,1:4])
print(table(signaltop_s2$cluster))
gaps<-cumsum(table(signaltop_s2$cluster))
gaps<-gaps[1:(length(gaps)-1)]

anno_row<-data.frame(tag=signaltop_s2$cluster)
rownames(anno_row)<-rownames(signaltop_s2)
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
anno_color<-list(tag=stallion_pal[1:length(table(signaltop_s2$cluster))])
names(anno_color[[1]])<-c(1:length(table(signaltop_s2$cluster)))

library(pheatmap)
pdf("heatmap_allPeaks13.pdf")
pheatmap(plotdat,breaks=seq(min(plotdat),quantile(plotdat,0.95),length.out=101),cluster_cols=F,cluster_rows=F,show_rownames=F,gaps_row=gaps,
annotation_row=anno_row,annotation_colors=anno_color)
dev.off()

saveRDS(kc,"kmeansResults_allPeaks13.rds")

#add fold change and filter
signal.qnorm <- as.data.frame(mtx)

signal.qnorm <- signal.qnorm[rownames(signaltop_s2),]
signal.qnorm$cluster <- signaltop_s2$cluster

fc_cutoff=2.5

#FS
clustersig<-signal.qnorm[which(signal.qnorm$cluster==4),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(2,3,4)])
    foreground=x[1]
    return(foreground/background)
})
difpeak_FS <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#HS
clustersig<-signal.qnorm[which(signal.qnorm$cluster==10),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,3,4)])
    foreground=x[2]
    return(foreground/background)
})
difpeak_HS <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#SI
clustersig<-signal.qnorm[which(signal.qnorm$cluster==2),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,2,4)])
    foreground=x[3]
    return(foreground/background)
})
difpeak_SI <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#CL
clustersig<-signal.qnorm[which(signal.qnorm$cluster==5),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,2,3)])
    foreground=x[4]
    return(foreground/background)
})
difpeak_CL <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#Stm
clustersig<-signal.qnorm[which(signal.qnorm$cluster==3 | signal.qnorm$cluster==7),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(3,4)])
    foreground=mean(x[c(1,2)])
    return(foreground/background)
})
difpeak_Stm <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#Int
clustersig<-signal.qnorm[which(signal.qnorm$cluster==12),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,2)])
    foreground=mean(x[c(3,4)])
    return(foreground/background)
})
difpeak_Int <- clustersig[which(clustersig$foldchange>=fc_cutoff),]


nrow(difpeak_FS)
nrow(difpeak_HS)
nrow(difpeak_SI)
nrow(difpeak_CL)
nrow(difpeak_Stm)
nrow(difpeak_Int)
difpeaks <- list(difpeak_FS,difpeak_HS,difpeak_SI,difpeak_CL,difpeak_Stm,difpeak_Int)

signaltop_s3<-signaltop_s2[c(rownames(difpeak_FS),rownames(difpeak_HS),rownames(difpeak_SI),rownames(difpeak_CL),rownames(difpeak_Stm),rownames(difpeak_Int)),]

plotdat <- as.matrix(signaltop_s3[,1:4])
gaps<-cumsum(lapply(difpeaks,nrow))
gaps<-gaps[1:(length(gaps)-1)]

anno_row<-data.frame(tag=signaltop_s3$cluster)
rownames(anno_row)<-rownames(signaltop_s3)
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
anno_color<-list(tag=stallion_pal[1:length(table(signaltop_s3$cluster))])
names(anno_color[[1]])<-c(names(table(signaltop_s3$cluster)))

library(pheatmap)
pdf("heatmap_difPeaks13.pdf")
pheatmap(plotdat,breaks=seq(min(plotdat),quantile(plotdat,0.95),length.out=101),cluster_cols=F,cluster_rows=F,show_rownames=F,gaps_row=gaps,
annotation_row=anno_row,annotation_colors=anno_color)
dev.off()


#write out dif peak
#read in peak file
t1 <- read.table("../narrowpeak/mergedPeak.bed",header=F,stringsAsFactors=F)
rownames(t1) <- paste0(t1$V1,":",t1$V2,"-",t1$V3)

t1.out<-t1[rownames(difpeak_FS),]
write.table(t1.out,"difPeaks_kmeans/E16_FS_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_HS),]
write.table(t1.out,"difPeaks_kmeans/E16_HS_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_SI),]
write.table(t1.out,"difPeaks_kmeans/E16_SI_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_CL),]
write.table(t1.out,"difPeaks_kmeans/E16_CL_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_Stm),]
write.table(t1.out,"difPeaks_kmeans/E16_Stm_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_Int),]
write.table(t1.out,"difPeaks_kmeans/E16_Int_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')
