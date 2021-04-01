motifoverlaps=paste0("/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1129/motif_overlap/",c("intestine","pancreas","stomach","lung"),"_peak_motifoverlap.txt")
cor_dir="/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1129/cor/"
cor_files=c("E95intestine_E135intestine_cor.txt","E95pancreas_E135pancreas_cor.txt","E95stomach_E135stomach_cor.txt","E95lung_E135lung_cor.txt")
plot_dir="/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1129/cor/plots/"
peaksets=paste0("/scratch/hz9fq/toronto _project/E95E135_CrossStages/mergedPeaks_1129/",c("intestine","pancreas","stomach","lung"),"_mergedPeak.bed")

for(set in 1:4){
    motifoverlap <- motifoverlaps[set]
    peakset <- peaksets[set]
    file <- cor_files[set]

    peak.tab <- read.table(peakset,header=F,stringsAsFactors=F)
    motif.tab <- read.table(motifoverlap,header=T,stringsAsFactors=F)
    peakNumber <- apply(motif.tab,2,sum)

    names(peakNumber) <- sapply(names(peakNumber),function(x){
        tmp <- strsplit(x,"_map")
        return(tmp[[1]][1])
    })

    print(file)
    t2 <- read.table(paste0(cor_dir,file),header=T,stringsAsFactors=F)
    rownames(t2) <- t2[,1]

    #set -300 as lower boundary
    t2$wilcox.new[which(t2$wilcox.new<(-300))] <- -300
    t2$wilcox.new[which(t2$wilcox.new>(300))] <- 300

    #wilcox.greater
    #datatype = "_wilcox.greater"
    #t2_use <- t2[,c(1,5)]
    #t2_use[,2] <- -log10(t2_use[,2])

    #wilcox.less
    #datatype = "_wilcox.less"
    #t2_use <- t2[,c(1,6)]
    #t2_use[,2] <- -log10(t2_use[,2])

    #wilcox.new
    datatype = "_wilcox.new"
    t2_use <- t2[,c(1,8)]
    
    uplimit = max(abs(t2_use[,2]))
    lowlimit = -uplimit

    t2_use$peaks<-as.numeric(peakNumber[rownames(t2_use)])
    
    blocks <- seq(min(t2_use$peaks),max(t2_use$peaks),length.out=101)
    blockctrs <- c()
    for(i in 1:(length(blocks)-1)){
        blockctrs <- c(blockctrs,mean(c(blocks[i],blocks[i+1])))
    }
    names(blockctrs) <- 1:100

    getblockctr <- function(x){
     return(names(blockctrs)[which.min(abs(blockctrs-x))])
    }
    t2_use$V4 <- sapply(t2_use$peaks,getblockctr)

    filename <- strsplit(file,"_cor")[[1]][1]
    
    #write.table(t2_use,file=file,col.names=F,row.names=F,sep='\t',quote=F)
    
    horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")
    colors <- colorRampPalette(horizon_pal)(100)
    color.V <- colors[as.numeric(t2_use$V4)]

    markers=unique(c(rownames(t2_use)[order(t2_use[,2],decreasing=T)][1:20],rownames(t2_use)[order(t2_use[,2],decreasing=F)][1:20]))

    t2_use$rank <- rank(t2_use$wilcox.new,ties.method="first")
    
    pdf(paste0(plot_dir,filename,datatype,".pdf"))
    par(mar=c(6,6,6,6))
    plot(t2_use[,2],t2_use[,5],xlab="wilcox_newValue",ylab="rank",type="p",pch=19,cex=.7,col=color.V,main=paste0(filename,datatype),xlim=c(lowlimit-10,uplimit+10))
    labels <- c(min(t2_use$peaks),round((min(t2_use$peaks)+max(t2_use$peaks))/2),max(t2_use$peaks))
    plotrix::color.legend(360,100,400,250,legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    dev.off()

    pdf(paste0(plot_dir,filename,datatype,"_ref.pdf"))
    par(mar=c(6,6,6,6))
    plot(t2_use[,2],t2_use[,5],xlab="wilcox_newValue",ylab="rank",type="p",pch=19,cex=.7,col=color.V,main=paste0(filename,datatype),xlim=c(lowlimit-10,uplimit+10))
    labels <- c(min(t2_use$peaks),round((min(t2_use$peaks)+max(t2_use$peaks))/2),max(t2_use$peaks))
    plotrix::color.legend(360,100,400,250,legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    text(t2_use[markers[20],2],t2_use[markers[20],5],"X",col="red")
    text(t2_use[markers[40],2],t2_use[markers[40],5],"X",col="red")
    dev.off()
    
    #barplot
    t2_use.order <- t2_use[order(t2_use$wilcox.new,decreasing=T),]
    candi <- t2_use.order[c(1:20,(nrow(t2_use.order)-19):nrow(t2_use.order)),]

    plotV <- candi$wilcox.new
    names(plotV) <- rownames(candi)
    plotV <- plotV[order(plotV)]
    
    pdf(paste0(plot_dir,filename,datatype,".barplot.pdf"))
    par(mar=c(4,8,4,8))
    barplot(plotV,horiz=T,col=c(rep("#00bfc4",20),rep("#f8766d",20)),las=1,xlim=c(-300,300),cex.names=0.5)
    dev.off()

}

