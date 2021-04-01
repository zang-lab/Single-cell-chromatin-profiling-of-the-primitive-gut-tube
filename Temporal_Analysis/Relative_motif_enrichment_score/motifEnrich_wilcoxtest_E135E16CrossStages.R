#######################stomach
#signal file pair
filename.li <- list(c("GSM394318_E13_Stomach_ATAC_r1.peak_1bin.bed","GSM394320_E16_Forestomach_ATAC_r1.peak_1bin.bed"),c("GSM394318_E13_Stomach_ATAC_r1.peak_1bin.bed","GSM394321_E16_Hindstomach_ATAC_r1.peak_1bin.bed"),c("GSM394320_E16_Forestomach_ATAC_r1.peak_1bin.bed","GSM394323_E16_stomSox2KO_ATAC_r1.peak_1bin.bed"),c("GSM394321_E16_Hindstomach_ATAC_r1.peak_1bin.bed","GSM394323_E16_stomSox2KO_ATAC_r1.peak_1bin.bed"))
names(filename.li) <- c("E13stomach_E16forestomach","E13stomach_E16hindstomach","E16forestomach_E16stomachSox2KO","E16hindstomach_E16stomachSox2KO")

signal_dir="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis_1129/"
motifoverlap="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis_1129/peak_motifoverlap.txt"
cor_dir="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/stomach_CrossStages/signal_forMotifAnalysis_1129/cor/"
peakset="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/stomach_CrossStages/peaks/mergedPeak_500bp.bed"

#######################intestine
#signal file pair
filename.li <- list(c("GSM394317_E13_int_ATAC_r1.peak_1bin.bed","GSM394322_E16_int_ATAC_r1.peak_1bin.bed"),c("GSM394317_E13_int_ATAC_r1.peak_1bin.bed","GSM394319_E16_Colon_ATAC_r1.peak_1bin.bed"),c("GSM394322_E16_int_ATAC_r1.peak_1bin.bed","GSM3175495_E165_epithelialCdx2KO_ATAC_combined.peak_1bin.bed"),c("GSM394319_E16_Colon_ATAC_r1.peak_1bin.bed","GSM3175495_E165_epithelialCdx2KO_ATAC_combined.peak_1bin.bed"))
names(filename.li) <- c("E13intestine_E16smallIntestine","E13intestine_E16colon","E16smallIntestine_E16epithelialCdx2KO","E16colon_E16epithelialCdx2KO")

signal_dir="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis_1129/"
motifoverlap="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis_1129/peak_motifoverlap.txt"
cor_dir="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/intestine_CrossStages/signal_forMotifAnalysis_1129/cor/"
peakset="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/intestine_CrossStages/peaks/mergedPeak_500bp.bed"

set=4
filename <- filename.li[[set]]
outputname <- names(filename.li)[set]

motif.tab <- read.table(motifoverlap,header=T,stringsAsFactors=F)
t1 <- read.table(paste0(signal_dir,filename[1]),header=F,stringsAsFactors=F)
t2 <- read.table(paste0(signal_dir,filename[2]),header=F,stringsAsFactors=F)

dat <- as.data.frame(matrix(c(t1$V4,t2$V4),ncol=2))
dat$max <- apply(dat,1,max)
rownames(dat) <- paste0(t1$V1,":",t1$V2,"-",t1$V3)

#filter peak
#dat_f1 <- dat[which(dat$max>=0.2),]

#rownames(motif.tab) <- paste0(t1$V1,":",t1$V2,"-",t1$V3)
#motif_f1.tab <- motif.tab[rownames(dat_f1),]
dat_f1 <- dat
motif_f1.tab <- motif.tab

#fold change and motif heatmap
#get logfc
dat_f1$fc <- (dat_f1$V2+0.1)/(dat_f1$V1+0.1)
dat_f1$log2fc <- log2(dat_f1$fc)

peak.tab <- read.table(peakset,header=F,stringsAsFactors=F)
rownames(peak.tab) <- paste0(peak.tab$V1,":",peak.tab$V2,"-",peak.tab$V3)

#filter enhancer or promoter
#dat_f1.use <- dat_f1[rownames(peak.tab)[which(peak.tab[,ncol(peak.tab)]==0)],]
#motif_f1.tab.use <- motif_f1.tab[rownames(peak.tab)[which(peak.tab[,ncol(peak.tab)]==0)],]
dat_f1.use <- dat_f1
motif_f1.tab.use <- motif_f1.tab

cor.spearman <- c()
cor.pearson <- c()
cor.test <- c()
wilcox.great <- c()
wilcox.less <- c()
wilcox.two <- c()

count <- 0
for(motif in colnames(motif_f1.tab)){
    roc.table=data.frame(motif=motif_f1.tab.use[,motif],log2fc=dat_f1.use$log2fc)

    cor.spearman <- c(cor.spearman,cor(roc.table$log2fc,roc.table$motif,method='spearman'))
    res1 <- cor.test(roc.table$log2fc,roc.table$motif,method='pearson')
    cor.pearson <- c(cor.pearson,as.numeric(res1$estimate))
    cor.test <- c(cor.test,as.numeric(res1$p.value))

    motif_fc <- roc.table[which(roc.table$motif==1),2]
    nomotif_fc <- roc.table[which(roc.table$motif==0),2]
    res2 <- wilcox.test(motif_fc,nomotif_fc,alternative="greater")
    wilcox.great <- c(wilcox.great,res2$p.value)
    res3 <- wilcox.test(motif_fc,nomotif_fc,alternative="less")
    wilcox.less <- c(wilcox.less,res3$p.value)
    res4 <- wilcox.test(motif_fc,nomotif_fc,alternative="two.sided")
    wilcox.two <- c(wilcox.two,res4$p.value)

    count <- count+1
    if(count%%10==0){
        print(count)
    }
}

cor.tab <- data.frame(motif=colnames(motif_f1.tab),cor.spearman=cor.spearman,cor.pearson=cor.pearson,cor.test=cor.test,wilcox.great=wilcox.great,wilcox.less=wilcox.less,wilcox.two=wilcox.two)

wilcox.new <- c()
for(i in 1:nrow(cor.tab)){
    if(cor.tab[i,5]<0.5){
        wilcox.new <- c(wilcox.new,-log10(cor.tab[i,7]))
    }else{
        wilcox.new <- c(wilcox.new,log10(cor.tab[i,7]))
    }
}
#wilcox.new[which(wilcox.new>300)] <- 300
#wilcox.new[which(wilcox.new<(-300))] <- (-300)

cor.tab$wilcox.new <- wilcox.new
cor.tab$motif <- sapply(cor.tab$motif,function(x){
    tmp <- strsplit(x,"_map")
    return(tmp[[1]][1])
})

print(head(cor.tab))

#cor.tab <- cor.tab[order(cor.tab$V2,decreasing=T),]
write.table(cor.tab,paste0(cor_dir,outputname,"_cor.txt"),row.names=F,col.names=T,sep='\t',quote=F)
