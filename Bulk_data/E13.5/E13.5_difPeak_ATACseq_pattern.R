setwd("/scratch/hz9fq/toronto_project/E13bulkATACpeak/signal/difPeaks_kmeans/signal_alldifPeak")

organ=c("intestine","pancreas","stomach","lung")
colors=c("#D51F26","#D5b60A","#C06CAB","#89C75F")

for(i in 1:4){
    print(organ[i])
    t1 <- read.table(paste0("E13.5_",organ[i],"_atac.peak2kbNearby_100bp.bed"),header=F)
    gaps <- cumsum(c(4734,4828,6686,4882))[1:3]
    plotdat <- t1[,4:43]

    library(pheatmap)
    png(paste0("plots/",organ[i],"_pattern_signal.png"),width=2000,height=2000,res=200,bg = "white")
    pheatmap(plotdat, color = colorRampPalette(colors = c('white',colors[i]))(100), cluster_rows=F,cluster_cols=F,show_rownames=F,show_colnames=F,breaks=seq(0,quantile(as.matrix(plotdat),0.99),length.out=101),gaps_row=gaps,cellwidth=3)
    dev.off()
}

