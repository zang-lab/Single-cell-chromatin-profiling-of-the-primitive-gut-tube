motifoverlap="/scratch/hz9fq/toronto_project/intestine_CrossStages/motif_overlap/peak_motifoverlap3.txt"
cor_dir="/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis3/cor/"
plot_dir="/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis3/cor/wilcox_plots/"
peakset="/scratch/hz9fq/toronto_project/intestine_CrossStages/peaks/mergedPeak_400bp.bed"

comparisons <- list(c("E13intestine_E16smallIntestine_cor.txt","E16smallIntestine_E16epithelialCdx2KO_cor.txt","intestine_comparison3"), c("E13intestine_E16colon_cor.txt","E16colon_E16epithelialCdx2KO_cor.txt","intestine_comparison4"), c("E13intestine_E16smallIntestine_cor.txt","E13intestine_E16colon_cor.txt","intestine_comparison5"))

motifoverlap="/scratch/hz9fq/toronto_project/stomach_CrossStages/motif_overlap/peak_motifoverlap3.txt"
cor_dir="/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis3/cor/"
plot_dir="/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis3/cor/wilcox_plots/"
peakset="/scratch/hz9fq/toronto_project/stomach_CrossStages/peaks/mergedPeak_400bp.bed"

comparisons <- list(c("E13stomach_E16forestomach_cor.txt","E16forestomach_E16stomachSox2KO_cor.txt","stomach_comparison3"), c("E13stomach_E16hindstomach_cor.txt","E16hindstomach_E16stomachSox2KO_cor.txt","stomach_comparison4"), c("E13stomach_E16forestomach_cor.txt","E13stomach_E16hindstomach_cor.txt", "stomach_comparison5"))

motif.tab <- read.table(motifoverlap,header=T,stringsAsFactors=F)
peakNumber <- apply(motif.tab,2,sum)

for(cp in 1:length(comparisons)){
	filename1 = comparisons[[cp]][1]
    filename2 = comparisons[[cp]][2]
    t1 <- read.table(paste0(cor_dir,filename1),header=T,stringsAsFactors=F)
    t2 <- read.table(paste0(cor_dir,filename2),header=T,stringsAsFactors=F)
    rownames(t1) <- t1[,1]
    rownames(t2) <- t2[,1]
    t1 <- t1[which(substr(rownames(t1),1,2)!="E1"),]
    t2 <- t2[which(substr(rownames(t2),1,2)!="E1"),]

    #wilcox.great
    #datatype = "_wilcox.great"
    #t1_use <- t1[,c(1,5)]
    #t2_use <- t2[,c(1,5)]
    #t1_use[,2] <- -log10(t1_use[,2])
    #t2_use[,2] <- -log10(t2_use[,2])

    #wilcox.less
    #datatype = "_wilcox.less"
    #t1_use <- t1[,c(1,6)]
    #t2_use <- t2[,c(1,6)]
    #t1_use[,2] <- -log10(t1_use[,2]+1e-5)
    #t2_use[,2] <- -log10(t2_use[,2]+1e-5)

    #wilcox.new
    datatype = "_wilcox.new"
    t1_use <- t1[,c(1,8)]
    t2_use <- t2[,c(1,8)]

    t1_use$V3 <- t2_use[rownames(t1_use),2]
    t1_use$V4<-as.numeric(peakNumber[rownames(t1_use)])

    blocks <- seq(min(t1_use$V4),max(t1_use$V4),length.out=101)
    blockctrs <- c()
    for(i in 1:(length(blocks)-1)){
        blockctrs <- c(blockctrs,mean(c(blocks[i],blocks[i+1])))
    }
    names(blockctrs) <- 1:100

    getblockctr <- function(x){
     return(names(blockctrs)[which.min(abs(blockctrs-x))])
    }
    t1_use$V5 <- sapply(t1_use$V4,getblockctr)

    horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")
    colors <- colorRampPalette(horizon_pal)(100)
    color.V <- colors[as.numeric(t1_use$V5)]

    uplimit = max(abs(t1_use[,2]),abs(t1_use[,3]))+5
    lowlimit = -uplimit

    pdf(paste0(plot_dir,comparisons[[cp]][3],datatype,".pdf"))
    par(mar=c(6,6,6,6))
    plot(t1_use[,2],t1_use[,3],xlab=strsplit(filename1,"_cor")[[1]][1],ylab=strsplit(filename2,"_cor")[[1]][1],type="p",pch=16,col=color.V,main="wilcox.new comparison",xlim=c(lowlimit,uplimit),ylim=c(lowlimit,uplimit))
    markers=unique(c(rownames(t1_use)[which(t1_use[,2]>=quantile(t1_use[,2],0.95))],rownames(t1_use)[which(t1_use$V2<=quantile(t1_use[,2],0.05))],rownames(t1_use)[which(t1_use[,3]>=quantile(t1_use[,3],0.95))],rownames(t1_use)[which(t1_use[,3]<=quantile(t1_use[,3],0.05))]))
    #,rownames(t1_use)[which(abs(t1_use[,3]-t1_use[,2])>=0.1)]
    abline(h=0,col="grey",lty=2,lwd=1.5)
	abline(v=0,col="grey",lty=2,lwd=1.5)
    labels <- c(min(t1_use$V4),round((min(t1_use$V4)+max(t1_use$V4))/2),max(t1_use$V4))
    plotrix::color.legend(340,100,380,250,legend=labels,rect.col=colors,align="rb",gradient="y",cex=0.7)
    for(m in markers){
        text(t1_use[m,2],t1_use[m,3]+10,m,cex=0.7,col="black")
    }
    dev.off()

    if(cp==1 | cp==2){
        print(comparisons[[cp]][3])
        print("right most")
        print(unique(c(rownames(t1_use)[which(t1_use[,2]>=quantile(t1_use[,2],0.95))],rownames(t1_use)[which(t1_use[,3]<=quantile(t1_use[,3],0.05))])))
        print("left most")
        print(unique(c(rownames(t1_use)[which(t1_use[,3]>=quantile(t1_use[,3],0.95))],rownames(t1_use)[which(t1_use$V2<=quantile(t1_use[,2],0.05))])))
    }else{
        print(comparisons[[cp]][3])
        print("right most")
        print(unique(c(rownames(t1_use)[which(t1_use[,2]>=quantile(t1_use[,2],0.95))],rownames(t1_use)[which(t1_use[,3]>=quantile(t1_use[,3],0.95))])))
        print("left most")
        print(unique(c(rownames(t1_use)[which(t1_use[,3]<=quantile(t1_use[,3],0.05))],rownames(t1_use)[which(t1_use$V2<=quantile(t1_use[,2],0.05))])))
    }
}
