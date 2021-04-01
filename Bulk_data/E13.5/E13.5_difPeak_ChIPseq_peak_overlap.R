t2 <- read.table("peaks/Cdx2_combine_m14_peaks.bed",header=F)
t2 <- read.table("peaks/Cdx2_combine_m14_peaks_withmotif.bed",header=F)
t2 <- read.table("peaks/CDX2_map.bed",header=F)

#t2 <- read.table("peaks/GSM4699973_E16_forestomachSox2_ChIP_r1_m14_peaks.bed",header=F)
#t2 <- read.table("peaks/GSM4699973_E16_forestomachSox2_ChIP_r1_m14_peaks_withmotif.bed",header=F)
#t2 <- read.table("peaks/SOX2_map.bed",header=F)

library(GenomicRanges)

gr1 <- GRanges(
    seqnames = factor(t1$V1),
    ranges = IRanges(start = t1$V2,end = t1$V3)
)
gr2 <- GRanges(
    seqnames = factor(t2$V1),
    ranges = IRanges(start = t2$V2,end = t2$V3)
)
res1=countOverlaps(gr1, gr2)
res1[which(res1>1)] <- 1

library(pheatmap)
pdf("plots/E13difPeaks_E13Cdx2ChIPPeak_withmotif_overlap_heatmap.pdf")
pheatmap(res1, color=c("white","red"), cellwidth=50, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=F, gaps_row = gaps)
dev.off()

library(pheatmap)
pdf("plots/E13difPeaks_Sox2motif_overlap_heatmap.pdf")
pheatmap(res1, color=c("white","red"), cellwidth=50, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=F, gaps_row = gaps)
dev.off()

countV <- c(sum(res1[1:gaps[1]]),
sum(res1[(gaps[1]+1):gaps[2]]),
sum(res1[(gaps[2]+1):gaps[3]]),
sum(res1[(gaps[3]+1):length(res1)]))

pdf("plots/E13difPeaks_E16Sox2ChIPPeak_overlap_barplot.pdf")
par(mar=c(6,8,6,8))
barplot(rev(countV),col="grey",horiz=T,axes=F,space=c(0,1,1,1))
axis(1,at=c(0,max(countV)), labels=c(0,max(countV)), las=1, cex.axis=0.6)
axis(2,at=c(0.5,2.5,4.5,6.5), labels=rev(c("intestine","pancreas","stomach","lung")), las=2, col = NA, col.ticks = NA, cex.axis=0.6)
dev.off()
