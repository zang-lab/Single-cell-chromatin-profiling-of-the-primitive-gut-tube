########################kmeans call dif peaks
filenames=c("GSM394317_E13_int_ATAC_r1.peak_1bin.bed","E135_panc_atac.peak_1bin.bed","GSM394318_E13_Stomach_ATAC_r1.peak_1bin.bed","E135_lung_atac.peak_1bin.bed")
signal=c()
for(fn in filenames){
    tmp.t <- read.table(fn,header=F,stringsAsFactors=F)
    signal<-c(signal,tmp.t$V4)
}
signal <- as.data.frame(matrix(signal,ncol=4))

rownames(signal) <- paste0(tmp.t$V1,":",tmp.t$V2,"-",tmp.t$V3)
colnames(signal) <- c("E135_intestine_ATAC","E135_pancreas_ATAC","E135_stomach_ATAC","E135_lung_ATAC")

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

signaltop <- as.data.frame(mtx)
signaltop <- log2(signaltop+1)

kc <- readRDS("kmeansResults_allPeaks13.rds")
#kc <- kmeans(signaltop,13,iter.max=20,nstart=50)
signaltop$cluster<-kc$cluster
signaltop <- signaltop[order(signaltop$cluster),]

plotdat <- as.matrix(signaltop[,1:4])
print(table(signaltop$cluster))
gaps<-cumsum(table(signaltop$cluster))
gaps<-gaps[1:(length(gaps)-1)]

anno_row<-data.frame(tag=signaltop$cluster)
rownames(anno_row)<-rownames(signaltop)
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
anno_color<-list(tag=stallion_pal[1:length(table(signaltop$cluster))])
names(anno_color[[1]])<-c(1:length(table(signaltop$cluster)))

library(pheatmap)
pdf("heatmap_allPeaks13.pdf")
pheatmap(plotdat,breaks=seq(min(plotdat),quantile(plotdat,0.95),length.out=101),cluster_cols=F,cluster_rows=F,show_rownames=F,gaps_row=gaps,
annotation_row=anno_row,annotation_colors=anno_color)
dev.off()

saveRDS(kc,"kmeansResults_allPeaks13.rds")

#add fold change and filter
signal.qnorm <- as.data.frame(mtx)

signal.qnorm <- signal.qnorm[rownames(signaltop),]
signal.qnorm$cluster <- signaltop$cluster

fc_cutoff=2.5
#stomach
clustersig<-signal.qnorm[which(signal.qnorm$cluster==9),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,2,4)])
    foreground=x[3]
    return(foreground/background)
})
difpeak_stom <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#pancreas
clustersig<-signal.qnorm[which(signal.qnorm$cluster==3),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,3,4)])
    foreground=x[2]
    return(foreground/background)
})
difpeak_panc <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#intestine
clustersig<-signal.qnorm[which(signal.qnorm$cluster==11),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(2,3,4)])
    foreground=x[1]
    return(foreground/background)
})
difpeak_inte <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

#lung
clustersig<-signal.qnorm[which(signal.qnorm$cluster==2),]
clustersig$foldchange <- apply(clustersig[,1:4],1,function(x){
    background=mean(x[c(1,2,3)])
    foreground=x[4]
    return(foreground/background)
})
difpeak_lung <- clustersig[which(clustersig$foldchange>=fc_cutoff),]

dim(difpeak_stom)
dim(difpeak_panc)
dim(difpeak_inte)
dim(difpeak_lung)

signaltop2<-signaltop[c(rownames(difpeak_inte),rownames(difpeak_panc),rownames(difpeak_stom),rownames(difpeak_lung)),]
signaltop2$cluster <- factor(signaltop2$cluster,level=c(9,3,11,2))

plotdat <- as.matrix(signaltop2[,1:4])
print(table(signaltop2$cluster))
gaps<-cumsum(c(nrow(difpeak_inte),nrow(difpeak_panc),nrow(difpeak_stom)))

anno_row<-data.frame(tag=signaltop2$cluster)
rownames(anno_row)<-rownames(signaltop2)
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
anno_color<-list(tag=stallion_pal[1:length(table(signaltop2$cluster))])
names(anno_color[[1]])<-c(names(table(signaltop2$cluster)))

library(pheatmap)
pdf("heatmap_difPeaks13.pdf")
pheatmap(plotdat,breaks=seq(min(plotdat),quantile(plotdat,0.95),length.out=101),cluster_cols=F,cluster_rows=F,show_rownames=F,gaps_row=gaps,
annotation_row=anno_row,annotation_colors=anno_color)
dev.off()


#write out dif peak
#read in peak file
t1 <- read.table("mergedPeak.bed",header=F,stringsAsFactors=F)
rownames(t1) <- paste0(t1$V1,":",t1$V2,"-",t1$V3)

t1.out<-t1[rownames(difpeak_inte),]
write.table(t1.out,"difPeaks_kmeans/intestine_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_panc),]
write.table(t1.out,"difPeaks_kmeans/pancreas_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_stom),]
write.table(t1.out,"difPeaks_kmeans/stomach_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

t1.out<-t1[rownames(difpeak_lung),]
write.table(t1.out,"difPeaks_kmeans/lung_difPeak_kmeans.bed",col.names=F,row.names=F,quote=F,sep='\t')

difpeaks <- rbind(difpeak_inte,difpeak_panc,difpeak_stom,difpeak_lung)
rangeDF <- tmp.t[,1:3]
rownames(rangeDF) <- paste0(rangeDF$V1,":",rangeDF$V2,"-",rangeDF$V3)

difpeaks_range <- rangeDF[rownames(difpeaks),]
difpeaks_range$cluster <- difpeaks$cluster
difpeaks_range[which(difpeaks_range$cluster==11),4] <- "intestine"
difpeaks_range[which(difpeaks_range$cluster==3),4] <- "pancreas"
difpeaks_range[which(difpeaks_range$cluster==9),4] <- "stomach"
difpeaks_range[which(difpeaks_range$cluster==2),4] <- "lung"

write.table(difpeaks_range,"difPeaks_kmeans/E135_organs_difPeaks.bed",col.names=F,row.names=F,sep="\t",quote=F)
