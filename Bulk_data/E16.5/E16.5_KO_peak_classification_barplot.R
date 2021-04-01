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

t1$fc <- fc
t1.sort <- t1[order(t1$fc,decreasing=T),]

peaksets <- list(t1.sort[which(t1.sort$fc>=2),],t1.sort[which(t1.sort$fc<2 & t1.sort$fc>=0.5),],t1.sort[which(t1.sort$fc<0.5),])

library(GenomicRanges)
organ=c("intestine","pancreas","lung","stomach")
#organ=c("SI","CL","intestine","FS","HS","stomach")
#organ=c("SI","CL","FS","HS")

props <- c()
for(i in 1:3){
    t1 <- peaksets[[i]]
    for(j in 1:length(organ)){
        og=organ[j]
        mergedPeak<-t1[,1:3]
        #t3<-read.table(paste0("../difPeak/E95_",og,"_difPeak.bed"),header=F,stringsAsFactors=F)
        t3<-read.table(paste0("../difPeak/E135_",og,"_difPeak.bed"),header=F,stringsAsFactors=F)
        #t3<-read.table(paste0("../difPeak/E16_",og,"_difPeak_kmeans.bed"),header=F,stringsAsFactors=F)
        #t3<-read.table(paste0("../narrowPeak/",og,"_m2_peaks.narrowPeak"),header=F,stringsAsFactors=F)
        #t3<-read.table(paste0("../narrowPeak/E135_",og,"_atac_m2_peaks.narrowPeak"),header=F,stringsAsFactors=F)
        #t3<-read.table(paste0("../narrowPeak/E165_",og,"_atac_m2_peaks.narrowPeak"),header=F,stringsAsFactors=F)

        gr1 <- GRanges(
            seqnames = factor(mergedPeak$V1),
            ranges = IRanges(start = mergedPeak$V2,end = mergedPeak$V3)
        )
        gr2 <- GRanges(
            seqnames = factor(t3$V1),
            ranges = IRanges(start = t3$V2,end = t3$V3)
        )

        res<-countOverlaps(gr1, gr2)
        props <- c(props,length(res[which(res>0)])/length(res))
    }
}

#E13.5 "#D51F26","#D3C035","#89C75F","#C06CAB"
#E16.5 "#89030D","#FFA07A","#89288F","#B2026F"
col.use <- c("#D51F26","#D3C035","#89C75F","#C06CAB")

upperlimits <- c()
prt<-list()
for(i in 0:2){
    labels=organ
    props.use <- props[(i*4+1):(i*4+4)]
    par(mar=c(10,4,10,4))
    barplot(props.use,space=c(0,0,0,0),col=col.use,axes=F,ylim=c(0,max(props.use)),horiz=F)
    axis(1,at=c(0.5,1.5,2.5,3.5),labels=labels,las=2)
    p <- pretty(par("usr")[3:4])  # find nice breaks at actual y range
    l <- formatC(p, format="f", digits=12)
    axis(2,at=p,labels=p,las=1)
    dev.off()

    if(max(props.use)>=max(p)){
        tmp <- (p[length(p)]+p[2]-p[1])
        p <- c(p,tmp)
        while(tmp < max(props.use)){
            tmp <- tmp+p[2]-p[1]
            p <- c(p,tmp)
        }
        upperlimits <- c(upperlimits,tmp)
        prt[[i+1]] <- p
    }else{
        upperlimits <- c(upperlimits,max(p))
        prt[[i+1]] <- p
    }
}

pdf("barplots/Cdx2KO_E135_barplot_ver.pdf")
for(i in 0:2){
    labels=organ
    props.use <- props[(i*4+1):(i*4+4)]
    par(mar=c(10,4,10,4))
    barplot(props.use,space=c(0,0,0,0),col=col.use,axes=F,ylim=c(0,upperlimits[i+1]),horiz=F)
    axis(1,at=c(0.5,1.5,2.5,3.5),labels=labels,las=2)
    #p <- pretty(par("usr")[3:4],n=4)  # find nice breaks at actual y range
    p <- prt[[i+1]]
    l <- formatC(p, format="f", digits=12)
    axis(2,at=p,labels=p,las=1)
}
dev.off()

