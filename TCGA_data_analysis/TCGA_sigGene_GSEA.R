#TCGA input data : decompressed FPKM file, decompressed count file
#generate the FPKM table for use
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/FPKM/COAD")

#transfer ensembl gene id to gene symbol
t1 <- read.table("9b30bf20-6b80-43f8-9205-67b990a6babf/a8a9595b-e0b5-479e-8b7e-33595de6f4c6.FPKM.txt",header=F,stringsAsFactors=F)

library("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

head(listAttributes(human, page="feature_page"))

id_list <- t1$V1
Gid_list <- sapply(id_list,function(x){
    return(strsplit(x,'\\.')[[1]][1])
})
Gid_list <- unname(Gid_list)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=Gid_list,mart= human)

G_list.use <- G_list[which(G_list$hgnc_symbol!=""),]

#filter out some repeated gene id and symbol (only retain the first of the repeats)
todel <- c()
rep_Gid <- names(table(G_list.use$ensembl_gene_id))[which(table(G_list.use$ensembl_gene_id)>1)]
for(id in rep_Gid){
    tmp <- G_list.use[which(G_list.use$ensembl_gene_id==id),]
    todel <- c(todel,rownames(tmp)[2:nrow(tmp)])
}
G_list.use <- G_list.use[which(!rownames(G_list.use)%in%todel),]

todel <- c()
rep_symbol <- names(table(G_list.use$hgnc_symbol))[which(table(G_list.use$hgnc_symbol)>1)]
for(sb in rep_symbol){
    tmp <- G_list.use[which(G_list.use$hgnc_symbol==sb),]
    todel <- c(todel,rownames(tmp)[2:nrow(tmp)])
}
G_list.use <- G_list.use[which(!rownames(G_list.use)%in%todel),]

rownames(G_list.use) <- G_list.use$ensembl_gene_id

#read in all the file by metadata
library("rjson")
meta <- fromJSON(file = "../Info/metadata.cart.COAD.json")

signal <- c()
count <- 0
case_id <- c()
patient_id <- c()
to_ignore <- c()

for(i in 1:length(meta)){
    tmp.meta <- meta[[i]]
    id <- (tmp.meta$associated_entities)[[1]]$entity_submitter_id
    code <- strsplit(id,"-")[[1]][4]
    if(!substr(code,1,2)%in%c("01")){
        to_ignore <- c(to_ignore,i)
        next
    }
    
    file <- list.files(tmp.meta$file_id,pattern="*.FPKM.txt")
    dat <- read.table(paste0(tmp.meta$file_id,"/",file))
    
    #check if gene id order is right
    if(!all(dat$V1==t1$V1)){
        print(tmp.meta$file_id)
        break
    }
    
    signal <- c(signal,dat$V2)
    case_id <- c(case_id,(tmp.meta$associated_entities)[[1]]$case_id)
    patient_id <- c(patient_id,substr(id,1,12))
    count <- count+1
    print(count)
}

FPKM_dat <- as.data.frame(matrix(signal,ncol=count))
rownames(FPKM_dat) <- Gid_list
colnames(FPKM_dat) <- patient_id

FPKM_dat_use <- FPKM_dat[which(rownames(FPKM_dat)%in%G_list.use[,1]),]
rownames(FPKM_dat_use) <- G_list.use[rownames(FPKM_dat_use),2]

write.table(FPKM_dat_use,"../COAD_FPKM.txt",row.names=T,col.names=T,sep='\t',quote=F)


#generate the read count table for use
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/count/COAD")

#transfer ensembl gene id to gene symbol
t1 <- read.table("ad71220b-8916-444a-96df-47602fed31f9/e86a03f9-e821-4f91-b09c-0e624239edde.htseq.counts.gz",header=F,stringsAsFactors=F)
t1 <- t1[1:(nrow(t1)-5),]

#library("biomaRt")
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

id_list <- t1$V1
Gid_list <- sapply(id_list,function(x){
    return(strsplit(x,'\\.')[[1]][1])
})
Gid_list <- unname(Gid_list)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=Gid_list,mart= human)

G_list.use <- G_list[which(G_list$hgnc_symbol!=""),]

#filter out some repeated gene id and symbol (only retain the first of the repeats)
todel <- c()
rep_Gid <- names(table(G_list.use$ensembl_gene_id))[which(table(G_list.use$ensembl_gene_id)>1)]
for(id in rep_Gid){
    tmp <- G_list.use[which(G_list.use$ensembl_gene_id==id),]
    todel <- c(todel,rownames(tmp)[2:nrow(tmp)])
}
G_list.use <- G_list.use[which(!rownames(G_list.use)%in%todel),]

todel <- c()
rep_symbol <- names(table(G_list.use$hgnc_symbol))[which(table(G_list.use$hgnc_symbol)>1)]
for(sb in rep_symbol){
    tmp <- G_list.use[which(G_list.use$hgnc_symbol==sb),]
    todel <- c(todel,rownames(tmp)[2:nrow(tmp)])
}
G_list.use <- G_list.use[which(!rownames(G_list.use)%in%todel),]

rownames(G_list.use) <- G_list.use$ensembl_gene_id

#read in all the file by metadata
library("rjson")
meta <- fromJSON(file = "../Info/metadata.cart.COAD.json")

signal <- c()
count <- 0
case_id <- c()
patient_id <- c()
to_ignore <- c()

for(i in 1:length(meta)){
    tmp.meta <- meta[[i]]
    id <- (tmp.meta$associated_entities)[[1]]$entity_submitter_id
    code <- strsplit(id,"-")[[1]][4]
    if(!substr(code,1,2)%in%c("01")){
        to_ignore <- c(to_ignore,i)
        next
    }
    
    file <- list.files(tmp.meta$file_id,pattern="*.htseq.counts")
    dat <- read.table(paste0(tmp.meta$file_id,"/",file))
    dat <- dat[1:(nrow(dat)-5),]
    
    #check if gene id order is right
    if(!all(dat$V1==t1$V1)){
        print(tmp.meta$file_id)
        break
    }
    
    signal <- c(signal,dat$V2)
    case_id <- c(case_id,(tmp.meta$associated_entities)[[1]]$case_id)
    patient_id <- c(patient_id,substr(id,1,12))
    count <- count+1
    print(count)
}

Count_dat <- as.data.frame(matrix(signal,ncol=count))
rownames(Count_dat) <- Gid_list
colnames(Count_dat) <- patient_id

Count_dat_use <- Count_dat[which(rownames(Count_dat)%in%G_list.use[,1]),]
rownames(Count_dat_use) <- G_list.use[rownames(Count_dat_use),2]

write.table(Count_dat_use,"../COAD_Count.txt",row.names=T,col.names=T,sep='\t',quote=F)


#get the distribution of key genes
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/FPKM")

FPKM_dat_use <- read.table("COAD_FPKM.txt",header=T,row.names=1,sep='\t',stringsAsFactors=F)

#get the zscore
Zscore_dat_use <- t(scale(t(FPKM_dat_use)))
write.table(Zscore_dat_use,"COAD_Zscore.txt",row.names=T,col.names=T,sep='\t',quote=F)

keys = c("SOX2","SOX9","CDX2")

pdf("../result_plots/COAD_key_genes_distribution.pdf")
for(k in keys){
    plotV <- as.vector(Zscore_dat_use[k,])
    plotV <- as.numeric(unname(plotV))
    print(hist(plotV,breaks=200,main=paste("Histogram of",k)))
}
dev.off()


#signature genes boxpot
library("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

library("rjson")

setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/FPKM")
FPKM_dat_use <- read.table("COAD_FPKM.txt",header=T,row.names=1,sep='\t',stringsAsFactors=F)
Zscore_dat_use <- read.table("COAD_Zscore.txt",header=T,row.names=1,sep='\t',stringsAsFactors=F)

#read in signature gene list
siglist <- read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/DEG_E16/limma_DEG_E16stomach_genelist.txt")
#siglist <- read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/DEG_kmeans/stomach_DEG_genelist.txt")
#siglist <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95_scRNA/markerGene/E95_stomach_markers.txt")

#transfer mouse gene symbol to human gene symbol
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = siglist$V1 , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

HGsiglist <- genesV2$HGNC.symbol

Zscore_sig <- Zscore_dat_use[which(rownames(Zscore_dat_use)%in%HGsiglist),]
print(paste0("signature genes exist in mouse gene : ",nrow(siglist)))
print(paste0("signature genes exist in human gene : ",nrow(Zscore_sig)))

Zscore_sig_avg <- apply(Zscore_sig,2,mean)

#pdf("../result_plots/COAD_key_genes_boxplot.pdf")
k="SOX2"
signal <- FPKM_dat_use[k,]
g1 <- which(signal<0.5)
g2 <- which(signal>=3)
print(c(length(g1),length(g2)))

glist <- Zscore_sig_avg
glist[g1] <- paste0(k,"_low")
glist[g2] <- paste0(k,"_high")    
glist <- factor(glist,levels=c(paste0(k,"_low"),paste0(k,"_high")))

glist <- glist[!is.na(glist)]

plotdat <- data.frame(Zscore_sig_avg[names(glist)],glist)
colnames(plotdat) <- c("Zscore_sig_avg","glist")
pal=c("#7fc97f","#fdc086")
boxplot(Zscore_sig_avg~glist,data=plotdat,outline=F,col=pal,ylab="Z score")
text(x=1,y=median(plotdat[which(plotdat$glist==paste0(k,"_low")),1])+0.015,labels=paste0("n = ",table(glist)[1]))
text(x=2,y=median(plotdat[which(plotdat$glist==paste0(k,"_high")),1])+0.015,labels=paste0("n = ",table(glist)[2]))

res <- t.test(Zscore_sig_avg[g2],Zscore_sig_avg[g1],alternative="greater")
print(res$p.value)

#dev.off()


#GSEA running enrichment
library("clusterProfiler")

#read in count matrix
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/count/")
Count_dat <- read.table("COAD_Count.txt",header=T,stringsAsFactors=F)
#rearrange count matrix by the order of case id in the Zscore matrix
Count_dat <- Count_dat[,colnames(Zscore_dat_use)]

#limma get differential genes
library(DESeq2)
library(edgeR)
set.seed(1)

#get the ranked gene list
#16.5
dir="/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/"
t2 <- read.table(paste0(dir,"FPKM_collection.txt"),header=T,stringsAsFactors=F)
#stomach expression
left_exp <- apply(t2[,15:18],1,mean)
#intestine expression
right_exp <- apply(t2[,19:22],1,mean)
exp_dat <- data.frame(left_exp,right_exp)
rownames(exp_dat) <- t2$Gene.ID

#13.5
dir="/scratch/hz9fq/toronto_project/ALL_RNA_collection/mod_exp/"
t2 <- read.table(paste0(dir,"FPKM_collection.txt"),header=T,stringsAsFactors=F)
#stomach expression
left_exp <- apply(t2[,11:12],1,mean)
#pancreas expression
#right_exp <- apply(t2[,9:10],1,mean)
#intestine expression
right_exp <- apply(t2[,7:8],1,mean)
exp_dat <- data.frame(left_exp,right_exp)
rownames(exp_dat) <- t2$Gene.ID

#9.5
dir="/scratch/hz9fq/toronto_project/E95_scRNA/"
t2 <- read.table(paste0(dir,"E95_Expression_matrix_byOrgan.txt"),header=T,stringsAsFactors=F)
#stomach expression
left_exp <- t2$stomach
#lung expression
#left_exp <- t2$lung
#pancreas expression
#right_exp <- t2$pancreas
#intestine expression
right_exp <- t2$intestine
exp_dat <- data.frame(left_exp,right_exp)
rownames(exp_dat) <- rownames(t2)

exp_dat$fc <- ((exp_dat$left_exp+0.1)/(exp_dat$right_exp+0.1))
exp_dat <- exp_dat[order(exp_dat$fc,decreasing=T),]
geneList <- exp_dat$fc
names(geneList) <- rownames(exp_dat)
#use log fold change so that there are negative values	
geneList<-log2(geneList)

k="SOX2"

coldata <- data.frame(compare1=factor(glist, levels=c(paste0(k,"_low"),paste0(k,"_high"))))
cnt_dds <- Count_dat[,rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData = cnt_dds,
                      colData = coldata,
                      design= ~ compare1)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res <- lfcShrink(dds, coef=paste0("compare1_",k,"_high","_vs_",k,"_low"), type="apeglm")
nrow(res[which(res$padj<=0.05 & res$log2FoldChange>=log2(1.5)),])
#filter dif peak
res.f <- res[which(res$padj<=0.05 & res$log2FoldChange>=log2(1.5)),]
res.f <- res.f[order(res.f$log2FoldChange,decreasing=T),]

siglist <- rownames(res.f)
#siglist <- rownames(res.f)[1:200]

#transfer human gene symbol to mouse gene symbol
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = siglist , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
GeneSet <- unique(genesV2$MGI.symbol)
print(length(GeneSet))

GeneSet<-GeneSet[which(GeneSet%in%names(geneList))]
print(length(GeneSet))
TTG <- data.frame(term=rep("term1",length(GeneSet)),gene=GeneSet)

#cluster profiler universe GSEA method
uni <- GSEA(
    geneList,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 5000,
    eps = 1e-100,
    pvalueCutoff = 10,
    pAdjustMethod = "BH",
    TERM2GENE = TTG,
)

print(uni@result)

pdf(paste0("../result_plots/COAD_GSEA_",k,"_test.pdf"))
p1 <- gseaplot(uni, geneSetID = 1, by = "runningScore", title=paste0(k," high signature genes"))
print(p1)
dev.off()

