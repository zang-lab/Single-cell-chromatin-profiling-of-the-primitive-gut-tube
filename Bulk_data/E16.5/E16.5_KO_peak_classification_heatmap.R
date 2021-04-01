######################task1 Cdx/Sox2 KO vs WT peaks overlap with Organ-specific peaks / narrowPeak
#list the merged peaks of WT and KO in the order of foldchange between them
setwd("/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Cdx2")
t1 <- read.table("E165_epithelialCdx2KO_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
t2 <- read.table("E165_epithelial_WTvsCdx2KO_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
#setwd("/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Sox2")
#t1 <- read.table("E16_stomSox2KO_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
#t2 <- read.table("E16_stom_WTvsSox2KO_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
rownames(t1) <- paste0(t1$V1,":",t1$V2,"-",t1$V3)
rownames(t2) <- paste0(t2$V1,":",t2$V2,"-",t2$V3)

fc <- (t1$V4+0.1)/(t2$V4+0.1)
orderPeaks <- rownames(t1)[order(fc,decreasing=T)]

#draw the heatmap of KO and WT
h1 <- read.table("E165_epithelialCdx2KO_ATAC.peak2kbNearby_100bp.bed",header=F,stringsAsFactors=F)
h2 <- read.table("E165_epithelial_WTvsCdx2KO_ATAC.peak2kbNearby_100bp.bed",header=F,stringsAsFactors=F)
#h1 <- read.table("E16_stomSox2KO_ATAC.peak2kbNearby_100bp.bed",header=F,stringsAsFactors=F)
#h2 <- read.table("E16_stom_WTvsSox2KO_ATAC.peak2kbNearby_100bp.bed",header=F,stringsAsFactors=F)
rownames(h1) <- paste0(h1$V1,":",h1$V2,"-",h1$V3)
rownames(h2) <- paste0(h2$V1,":",h2$V2,"-",h2$V3)
plotdat <- as.matrix(cbind(h1[orderPeaks,4:43],h2[orderPeaks,4:43]))

toplimit <- quantile(plotdat,0.99)

groups=c(length(fc[which(fc>=2)]),length(fc[which(fc>=0.5 & fc<2)]),length(fc[which(fc<0.5)]))
print(groups)
gaps=cumsum(groups)[1:2]

library(pheatmap)
pdf("../plots/4kbheatmap_Cdx2_WTandKO.pdf")
pheatmap(plotdat, breaks=seq(0,toplimit,length.out=101), color = colorRampPalette(colors = c('white','black'))(100), cellwidth=5, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=F, gaps_col = c(40), gaps_row = gaps)
dev.off()

library(pheatmap)
png("../plots/4kbheatmap_Cdx2_WTandKO.png",width=2000,height=2000,res=200,bg = "white")
pheatmap(plotdat, breaks=seq(0,toplimit,length.out=101), color = colorRampPalette(colors = c('white','black'))(100), cellwidth=5, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=F, gaps_col = c(40), gaps_row = gaps)
dev.off()

#get overlap peaks between merged peaks and cluster marker peaks and draw heatmap
library(GenomicRanges)
#organ=c("intestine","pancreas","lung","stomach")
organ=c("SI","CL","FS","HS")
overlap.V=c()
for(i in 1:length(organ)){
    og=organ[i]
    mergedPeak<-t1[,1:3]
    #t3<-read.table(paste0("../difPeak/E135_",og,"_difPeak.bed"),header=F,stringsAsFactors=F)
    t3<-read.table(paste0("../difPeak/E16_",og,"_difPeak_kmeans.bed"),header=F,stringsAsFactors=F)

    gr1 <- GRanges(
        seqnames = factor(mergedPeak$V1),
        ranges = IRanges(start = mergedPeak$V2,end = mergedPeak$V3)
    )
    gr2 <- GRanges(
        seqnames = factor(t3$V1),
        ranges = IRanges(start = t3$V2,end = t3$V3)
    )

    res<-countOverlaps(gr1, gr2)
    #each column has a different number for color
    res[which(res>=1)] <- i
    mergedPeak$V4<-res

    overlap.V <- c(overlap.V,mergedPeak[orderPeaks,4])
}

plotdat<-matrix(overlap.V,ncol=length(organ))
colnames(plotdat) <- organ

library(pheatmap)
png("../plots/overlapheatmap_Cdx2_E165difPeak.png",width=2000,height=2000,res=200)
#E13.5 "#D51F26","#D3C035","#89C75F","#C06CAB"
#E16.5 "#89030D","#FFA07A","#89288F","#B2026F"
pheatmap(plotdat, color=c("white","#89030D","#FFA07A","#89288F","#B2026F"), cellwidth=80, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=T, gaps_col=seq(1,(length(organ)-1),1), gaps_row = gaps,  legend=FALSE, fontsize_col=20)
dev.off()
