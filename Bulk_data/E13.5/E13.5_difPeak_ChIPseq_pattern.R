setwd("/scratch/hz9fq/toronto_project/E13bulkATACpeak/difPeaks_Sox2Cdx2/")

t1 <- read.table("100bp_signal/E135_mouseCdx2_ChIP.peak2kbNearby_100bp.bed",header=F)
table(t1$V4)
gaps <- cumsum(c(4734,4828,6686))
plotdat <- t1[,5:44]

library(pheatmap)
pdf("plots/Cdx2_ChIP_pattern.pdf")
pheatmap(plotdat, color = colorRampPalette(colors = c('white','black'))(100), cluster_rows=F,cluster_cols=F,show_rownames=F,show_colnames=F,breaks=seq(0.05,quantile(as.matrix(plotdat),0.99),length.out=101),gaps_row=gaps,cellwidth=3)
dev.off()

#png
library(pheatmap)
png("plots/Cdx2_ChIP_pattern.png",width=2000,height=2000,res=200,bg = "white")
pheatmap(plotdat, color = colorRampPalette(colors = c('white','black'))(100), cluster_rows=F,cluster_cols=F,show_rownames=F,show_colnames=F,breaks=seq(0.05,quantile(as.matrix(plotdat),0.99),length.out=101),gaps_row=gaps,cellwidth=3)
dev.off()


t1 <- read.table("100bp_signal/E16_forestomachSox2_ChIP.peak2kbNearby_100bp.bed",header=F)
gaps <- cumsum(c(4734,4828,6686))
plotdat <- t1[,5:44]

library(pheatmap)
pdf("plots/Sox2_ChIP_pattern.pdf")
pheatmap(plotdat, color = colorRampPalette(colors = c('white','black'))(100), cluster_rows=F,cluster_cols=F,show_rownames=F,show_colnames=F,breaks=seq(0.1,quantile(as.matrix(plotdat),0.95),length.out=101),gaps_row=gaps,cellwidth=3)
dev.off()

#png
library(pheatmap)
png("plots/Sox2_ChIP_pattern.png",width=2000,height=2000,res=200,bg = "white")
pheatmap(plotdat, color = colorRampPalette(colors = c('white','black'))(100), cluster_rows=F,cluster_cols=F,show_rownames=F,show_colnames=F,breaks=seq(0.1,quantile(as.matrix(plotdat),0.95),length.out=101),gaps_row=gaps,cellwidth=3)
dev.off()
