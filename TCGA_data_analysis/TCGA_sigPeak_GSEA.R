#write out specific cancer's pan cancer peak count matrix
all.norm.cnt <- readRDS("/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/TCGA-ATAC_PanCan_Log2Norm_Counts.rds")
norm.cnt <- all.norm.cnt[,which(substr(colnames(all.norm.cnt),1,4)=="PAAD")]
colnames(norm.cnt) <- gsub("_","-", colnames(norm.cnt))
#be aware the start/end is 1 greater than bed file format (start from 1)
norm.cnt <- cbind(all.norm.cnt[,1:7],norm.cnt)
saveRDS(norm.cnt,"/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/PAAD_PanCan_Log2Norm_Counts.rds")

all.cnt <- readRDS("/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/TCGA-ATAC_PanCan_Raw_Counts.rds")
cnt <- all.cnt[,which(substr(colnames(all.norm.cnt),1,4)=="PAAD")]
colnames(cnt) <- gsub("_","-", colnames(cnt))
cnt <- cbind(all.cnt[,1:7],cnt)
saveRDS(cnt,"/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/PAAD_PanCan_Raw_Counts.rds")


#signature peak boxpot
setwd("/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/")

norm.cnt <- readRDS("/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/STAD_PanCan_Log2Norm_Counts.rds")
#read in signature peak list
#siglist <- read.table("/scratch/hz9fq/toronto_project/E16bulkATACpeak/mergedPeak_signal/DESeq2_stomach_difPeak_hg38.bed")
#siglist <- read.table("/scratch/hz9fq/toronto_project/E13bulkATACpeak/signal/difPeaks_kmeans/stomach_difPeak_kmeans_hg38.bed")
siglist <- read.table("/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/markerPeak/markerPeak_7organs/markerPeak_bysuperCluster_stomach_hg38.bed")

library(GenomicRanges)
gr1 <- GRanges(
    seqnames = factor(norm.cnt[,1]),
    ranges = IRanges(start = norm.cnt[,2]-1,end = norm.cnt[,3]-1)
)
gr2 <- GRanges(
    seqnames = factor(siglist[,1]),
    ranges = IRanges(start = siglist[,2],end = siglist[,3])
)

res<-countOverlaps(gr1, gr2)
norm.cnt.use <- norm.cnt[which(res>0),]

print(paste0("signature peak exist in pan cancer peak set : ",nrow(norm.cnt.use)))

norm.cnt.use_avg <- apply(norm.cnt.use[,8:ncol(norm.cnt.use)],2,mean)

#get the FPKM
FPKM_dat_use <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/FPKM/STAD_FPKM.txt",header=T,row.names=1,sep='\t',stringsAsFactors=F)

#get ATAC metadata
library("rjson")
meta.atac <- fromJSON(file = "Info/metadata.cart.STAD_ATAC.json")
acaseID<-c()
apatientID<-c()
for(i in 1:length(meta.atac)){
    acaseID<-c(acaseID,meta.atac[[i]]$associated_entities[[1]]$case_id)
    id <- meta.atac[[i]]$associated_entities[[1]]$entity_submitter_id
    apatientID<-c(apatientID,substr(id,1,12))
}

#get RNA metadata
meta.rna <- fromJSON(file = "/scratch/hz9fq/toronto_project/TCGA/FPKM/Info/metadata.cart.STAD.json")
rcaseID<-c()
for(i in 1:length(meta.rna)){
    rcaseID<-c(rcaseID,meta.rna[[i]]$associated_entities[[1]]$case_id)
}

k = "SOX2"

signal <- FPKM_dat_use[k,]
r1 <- colnames(FPKM_dat_use)[which(signal<0.5)]
r2 <- colnames(FPKM_dat_use)[which(signal>=3)]
print(c(length(r1),length(r2)))

#get associated atac groups using the case id
case1<-gsub("\\.","-",r1)
case1<-sapply(case1,function(x){
    if(substr(x,1,1)=="X"){
        return(substr(x,2,nchar(x)))
    }else{
        return(x)
    }
})
case1 <- unname(case1)

case2<-gsub("\\.","-",r2)
case2<-sapply(case2,function(x){
    if(substr(x,1,1)=="X"){
        return(substr(x,2,nchar(x)))
    }else{
        return(x)
    }
})
case2 <- unname(case2)

a1_case <- apatientID[which(apatientID%in%case1)]
a2_case <- apatientID[which(apatientID%in%case2)]
print(c(length(a1_case),length(a2_case)))

a1 <- c()
for(i in which(apatientID%in%a1_case)){
    libquantity <- length(meta.atac[[i]]$analysis$metadata[[1]])
    for(j in 1:libquantity){
        a1 <- c(a1,meta.atac[[i]]$analysis$metadata[[1]][[j]]$library_name)
    }
}

a2 <- c()
for(i in which(apatientID%in%a2_case)){
    libquantity <- length(meta.atac[[i]]$analysis$metadata[[1]])
    for(j in 1:libquantity){
        a2 <- c(a2,meta.atac[[i]]$analysis$metadata[[1]][[j]]$library_name)
    }
}

gvalue <- norm.cnt.use_avg[which(names(norm.cnt.use_avg)%in%c(a1,a2))]
glist <- gvalue
glist[which(names(glist)%in%a1)] <- paste0(k,"_low")
glist[which(names(glist)%in%a2)] <- paste0(k,"_high") 
glist <- factor(glist,levels=c(paste0(k,"_low"),paste0(k,"_high")))
print(table(glist))

plotdat <- data.frame(gvalue,glist)
pal=c("#7fc97f","#fdc086")
boxplot(gvalue~glist,data=plotdat,outline=F,col=pal,ylab="normalized count")
text(x=1,y=median(plotdat[which(plotdat$glist==paste0(k,"_low")),1])+0.01,labels=paste0("n = ",table(glist)[1]))
text(x=2,y=median(plotdat[which(plotdat$glist==paste0(k,"_high")),1])+0.01,labels=paste0("n = ",table(glist)[2]))

res <- t.test(gvalue[a2],gvalue[a1],alternative="greater")
print(res$p.value)


#GSEA running enrichment
library("clusterProfiler")

#read in count matrix
cnt <- readRDS("STAD_PanCan_Raw_Counts.rds")
row.names(cnt) <- paste0(cnt[,1],":",(cnt[,2]-1),"-",(cnt[,3]-1))
cnt_use <- cnt[,8:ncol(cnt)]
pairtab <- read.table("pan_cancer_peakset_mm10_mappedPair.bed",header=F)

#limma get differential genes
library(DESeq2)
library(edgeR)
set.seed(1)

#get the ranked gene list
t1.1 <- read.table("pan_cancer_peak_signal/E16_FS_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
t1.2 <- read.table("pan_cancer_peak_signal/E16_HS_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
t2.1 <- read.table("pan_cancer_peak_signal/E16_SI_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
t2.2 <- read.table("pan_cancer_peak_signal/E16_CL_ATAC.peak_1bin.bed",header=F,stringsAsFactors=F)
t1 <- data.frame(t1.1[,1:3],apply(data.frame(t1.1$V4,t1.2$V4),1,mean))
colnames(t1)[4] <- "V4"
t2 <- data.frame(t2.1[,1:3],apply(data.frame(t2.1$V4,t2.2$V4),1,mean))
colnames(t2)[4] <- "V4"

t1 <- read.table("pan_cancer_peak_signal/E135_stomach.peak_1bin.bed",header=F,stringsAsFactors=F)
t2 <- read.table("pan_cancer_peak_signal/E135_intestine.peak_1bin.bed",header=F,stringsAsFactors=F)

t1 <- read.table("pan_cancer_peak_signal/stomach.peak_1bin.bed",header=F,stringsAsFactors=F)
t2 <- read.table("pan_cancer_peak_signal/intestine.peak_1bin.bed",header=F,stringsAsFactors=F,skipNul = TRUE)

exp_dat <- data.frame(t1$V4,t2$V4)
rownames(exp_dat) <- paste0(t1$V1,":",t1$V2,"-",t1$V3)
colnames(exp_dat) <- c("stm_exp","int_exp")

exp_dat$fc <- ((exp_dat$stm_exp+0.01)/(exp_dat$int_exp+0.01))
exp_dat <- exp_dat[order(exp_dat$fc,decreasing=T),]
geneList <- exp_dat$fc
names(geneList) <- rownames(exp_dat)
#use log fold change so that there are negative values	
geneList<-log2(geneList)

k = "SOX2"
coldata <- data.frame(compare1=factor(glist, levels=c(paste0(k,"_low"),paste0(k,"_high"))))
cnt_dds <- cnt_use[,rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData = cnt_dds,
                      colData = coldata,
                      design= ~ compare1)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
#res <- results(dds, name="compare4_stomach_vs_other")
res <- lfcShrink(dds, coef=paste0("compare1_",k,"_high","_vs_",k,"_low"), type="apeglm")
nrow(res[which(res$padj<=0.1 & res$log2FoldChange>=log2(1.5)),])
#filter dif peak
res.f <- res[which(res$padj<=0.1 & res$log2FoldChange>=log2(1.5)),]

#only remain peaks exist in mouse genome
res.f2 <- res.f[which(rownames(res.f)%in%pairtab$V4),]
res.f2 <- res.f2[order(res.f2$log2FoldChange,decreasing=T),]
siglist <- rownames(res.f2)
#siglist <- rownames(res.f2)[1:200]

#transformed hg38 peaks to mm10 peaks
geneset <- pairtab$V8[which(pairtab$V4%in%siglist)]
print(length(geneset))

TTG <- data.frame(term=rep("term1",length(geneset)),gene=geneset)

#cluster profiler universe GSEA method
uni <- GSEA(
    geneList,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 5000,
    eps = 1e-100,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = TTG,
)

print(uni@result)

pdf(paste0("../result_plots/STAD_E95_peakLevel_GSEA_",k,"_high.pdf"))
p1 <- gseaplot(uni, geneSetID = 1, by = "runningScore", title=paste0(k," high signature genes"))
print(p1)
dev.off()

