#############################task1 get dif gene of E13.5 vs E9.5
setwd("/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene")

E95.mtx <- read.table("E95_FPKM_matrix_4Organs.txt",header=T,row.names=1)
E135.mtx <- read.table("E135_FPKM_matrix_4Organs.txt",header=T,row.names=1)

gene.use <- intersect(rownames(E95.mtx),rownames(E135.mtx))
E95.mtx.use <- E95.mtx[gene.use,]
colnames(E95.mtx.use) <- paste0("E95_",colnames(E95.mtx.use))
E135.mtx.use <- E135.mtx[gene.use,]
colnames(E135.mtx.use) <- paste0("E135_",colnames(E135.mtx.use))

organs=c("intestine","pancreas","stomach","lung")

all.mtx <- cbind(E95.mtx.use,E135.mtx.use)
all.mtx.norm <- preprocessCore::normalize.quantiles(as.matrix(all.mtx))
all.mtx.norm <- as.data.frame(all.mtx.norm)
colnames(all.mtx.norm) <- colnames(all.mtx)
rownames(all.mtx.norm) <- rownames(all.mtx)

#get log2 fpkm
all.mtx.log <- log2(all.mtx.norm+1)

#get dif genes
for(g in organs){
    fore <- all.mtx.log[,paste0("E135_",g)]
    back <- all.mtx.log[,paste0("E95_",g)]
    fc <- (fore+0.5)/(back+0.5)
    names(fc) <- rownames(all.mtx.log)
    fc <- fc[order(fc,decreasing=T)]

    top <- fc[which(fc>2)]
    bottom <- fc[which(fc<(1/2))]

    print(c(length(top),length(bottom)))

    #write out dif genes
    out.dat <- data.frame(gene=names(top))
    write.table(out.dat,paste0("difgene_",g,"_up.genelist.txt"),row.names=F,col.names=F,quote=F,sep="\t")

    out.dat <- data.frame(gene=names(bottom))
    write.table(out.dat,paste0("difgene_",g,"_down.genelist.txt"),row.names=F,col.names=F,quote=F,sep="\t")
}

#test if pseudo count is good
test <- all.mtx.log[names(top),c(3,7)]
summary(test)
test <- all.mtx.log[names(bottom),c(3,7)]
summary(test)

#############################task2 get dif peak of E13.5 vs E9.5
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1204/")

organs=c("intestine","pancreas","stomach","lung")

sig <- c()
for(g in organs){
    E135_sig <- read.table(paste0("E135_",g,".peak_1bin.bed"),sep='\t',header=F)
    E95_sig <- read.table(paste0(g,".peak_1bin.bed"),header=F)
    sig <- c(sig,E135_sig$V4,E95_sig$V4)
}
dat <- matrix(sig,ncol=8)
#quantile normalization
dat.norm <- preprocessCore::normalize.quantiles(dat)
rownames(dat.norm) <- rownames(E135_sig)
colnames(dat.norm) <- c("E135_intestine","E95_intestine","E135_pancreas","E95_pancreas","E135_stomach","E95_stomach","E135_lung","E95_lung")

for(i in 1:length(organs)){
    g=organs[i]
    print(c(colnames(dat.norm)[2*i-1],colnames(dat.norm)[2*i]))

    fore <- log2(dat.norm[,2*i-1]+1)

    back <- log2(dat.norm[,2*i]+1)

    fc <- (fore+0.1)/(back+0.1)
    names(fc) <- rownames(E135_sig)
    fc <- fc[order(fc,decreasing=T)]

    top <- fc[which(fc>2)]
    bottom <- fc[which(fc<(1/2))]

    print(c(length(top),length(bottom)))

    #write out dif genes
    out.dat <- E135_sig[names(top),1:3]
    write.table(out.dat,paste0("difPeaks_QN/difPeak_",g,"_up.peaklist.txt"),row.names=F,col.names=F,quote=F,sep="\t")

    out.dat <- E135_sig[names(bottom),1:3]
    write.table(out.dat,paste0("difPeaks_QN/difPeak_",g,"_down.peaklist.txt"),row.names=F,col.names=F,quote=F,sep="\t")
}

##########################task3 get Association score based on distal regions of genes
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95E135_CrossStages/signal_mergedPeaks_1204/")
organs <- c("intestine","pancreas","stomach","lung")

t0 <- read.table("/nv/vol190/zanglab/hz9fq/annotation/Mouse/mm10_gene_annotation_geneID_LenOrder_TSS50kb.bed",header=F)

overlap.dat <- as.data.frame(matrix(rep(0,4*4),ncol=4))
rownames(overlap.dat) <- paste0("ATAC_",organs)
colnames(overlap.dat) <- paste0("RNA_",organs)

overlap.gene <- as.data.frame(matrix(rep(0,4*4),ncol=4))
rownames(overlap.gene) <- paste0("ATAC_",organs)
colnames(overlap.gene) <- paste0("RNA_",organs)

RNA_dir <- "/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95E135_CrossStages/RNA_difgene/"

#gene-centered confusion matrix
for(i in 1:length(organs)){
    t1 <- read.table(paste0("difPeaks_QN/difPeak_",organs[i],"_up.peaklist.txt"),header=F)
    #use top/bottom 2000 peaks
    t1 <- t1[1:2000,]
    #t1 <- t1[(nrow(t1)-1999):nrow(t1),]

    for(j in 1:length(organs)){
        print(c(i,j))

        geneList.dat <- read.table(paste0(RNA_dir,"difgene_",organs[j],"_up.genelist.txt"),header=F)
        geneList <- geneList.dat[,1]
        #use top/bottom 200 genes
        geneList <- geneList[1:200]
        #geneList <- geneList[(length(geneList)-199):length(geneList)]

        library(GenomicRanges)
        gr1 <- GRanges(
            seqnames = factor(t1[,1]),
            ranges = IRanges(start = t1[,2],end = t1[,3])
        )
        gr2 <- GRanges(
            seqnames = factor(t0[,1]),
            ranges = IRanges(start = t0[,2],end = t0[,3])
        )

        res<-countOverlaps(gr2, gr1)

        gene1 <- t0$V5[which(res>0)]
        gene2 <- geneList[which(geneList%in%gene1)]

        overlap.dat[i,j] <- length(gene2)/length(geneList)
        #overlap.dat[i,j] <- length(gene2)

        overlap.gene[i,j] <- length(gene2)

        print(c(length(gene1),length(gene2)))
        #print(c(organs[i],paste0("peak number:",nrow(t1)),paste0("gene number:",length(geneList)),paste0("over gene:",length(gene1))))

        #if(i==j){
            #write.table(gene2,paste0("up_genes/",organs[i],"_associated_up_genes.txt"),col.names=F,row.names=F,quote=F,sep="\t")
        #}
    }
}

libsize = c(2334,2195,1967,2091)
for(i in 1:4){
    overlap.dat[i,] <- (overlap.dat[i,]/libsize[i])*24528
}

overlap.norm <- overlap.dat
for(i in 1:4){
    overlap.norm[i,] <- (overlap.norm[i,]-min(overlap.norm[i,]))/(max(overlap.norm[i,])-min(overlap.norm[i,]))
}

library(pheatmap)
pheatmap(overlap.dat,cellwidth=40,cellheight=40,cluster_cols=F,cluster_rows=F,display_numbers = T,number_format = "%.3f")

display.dat <- round(overlap.dat,3)
pdf("oddsRatio_upPeak_upGene_heatmap.pdf")
pheatmap(overlap.norm,cellwidth=40,cellheight=40,cluster_cols=F,cluster_rows=F,display_numbers = display.dat,number_format = "%.3f")
dev.off()

overlap.gene.norm <- overlap.gene
for(i in 1:4){
    overlap.gene.norm[i,] <- overlap.gene.norm[i,]/max(overlap.gene.norm[i,])
}

pheatmap(overlap.gene.norm,cellwidth=40,cellheight=40,cluster_cols=F,cluster_rows=F,display_numbers = overlap.gene)
