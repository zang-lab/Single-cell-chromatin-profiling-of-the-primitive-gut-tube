#extract E9.5 DE cells from count matrix
setwd("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95_scRNA/matrics")

meta<-read.table("GSE136689_MetaFile_DE_and_SM_Cells_E8.5_E9.0_E9.5.txt",header=T,sep='\t',row.names=1,stringsAsFactors=F)
mtx<-read.table("GSE136689_Counts_Matrix_DE_and_SM_Cells_E8.5_E9.0_E9.5.txt",header=T,row.names=1,stringsAsFactors=F)

meta.use1<-meta[which(meta$Stages=="E9.5Anterior" & meta$Type=="Endoderm"),]
meta.use2<-meta[which(meta$Stages=="E9.5Posterior" & meta$Type=="Endoderm"),]
#anterior endoderm: 1512
#posterior endoderm: 1084

mtx.use <- mtx[,c(rownames(meta.use1),rownames(meta.use2))]
write.table(mtx.use,"Counts_Matrix_DE_E9.5.txt",col.names=T,row.names=T,sep='\t',quote=F)


#create seurate project
library(Seurat)
library(dplyr)
library(Matrix)

setwd("/scratch/hz9fq/toronto_project/E95_scRNA")

raw_counts<-read.table(file="matrics/Counts_Matrix_DE_E9.5.txt",sep="\t")
E95 <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features = 200, project = "E95_scRNA")

saveRDS(E95,"E95.rds")


#standard pre-processing
E95[["percent.mt"]] <- PercentageFeatureSet(E95, pattern = "^MT-")

pdf("QC_plots/matrixQC1.pdf")
par(mar=c(6,2,6,2))
VlnPlot(E95, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("QC_plots/matrixQC2.pdf")
plot1 <- FeatureScatter(E95, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E95, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1)
print(plot2)
dev.off()

#E95 <- subset(E95, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize the data
E95 <- NormalizeData(E95, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
E95 <- FindVariableFeatures(E95, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(E95), 10)

# plot variable features with and without labels
pdf("QC_plots/variableFeatures.pdf")
plot1 <- VariableFeaturePlot(E95)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1)
print(plot2)
dev.off()

#scale the data
all.genes <- rownames(E95)
E95 <- ScaleData(E95, features = all.genes)

#perform linear dimensional reduction
E95 <- RunPCA(E95, features = VariableFeatures(object = E95))

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
E95 <- JackStraw(E95, num.replicate = 100, dims = 30)
E95 <- ScoreJackStraw(E95, dims = 1:30)

pdf("QC_plots/JackStraw.pdf")
JackStrawPlot(E95, dims = 1:30)
dev.off()

pdf("QC_plots/Elbow.pdf")
ElbowPlot(E95,ndims = 30)
dev.off()

#cluster the cells
E95 <- FindNeighbors(E95, dims = 1:27)
E95 <- FindClusters(E95, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(E95), 5)

## If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
E95 <- RunUMAP(E95, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pdf("result_plots/UMAP.pdf")
DimPlot(E95, reduction = "umap")
dev.off()

saveRDS(E95,"E95.rds")

mtx<-E95@assays$RNA@scale.data
#colnames(mtx)<-gsub('\\.',"-",colnames(mtx))
umapmat<-E95@reductions$umap@cell.embeddings

stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
color_use<-stallion_pal[1:7]
colorV<-color_use[as.numeric(E95@meta.data$seurat_clusters)]

pdf("result_plots/E95_umap_cluster.pdf")
par(mar=c(8,4,6,10))
plot(umapmat,col=colorV,pch='.')
legend(11,6,legend=levels(E95@meta.data$seurat_clusters),fill=color_use,bty="n",xpd=TRUE)
dev.off()

#use organ defined by the paper to label cells
mtx<-E95@assays$RNA@scale.data
#colnames(mtx)<-gsub('\\.',"-",colnames(mtx))
umapmat<-E95@reductions$umap@cell.embeddings

meta <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95_scRNA/matrics/GSE136689_MetaFile_DE_and_SM_Cells_E8.5_E9.0_E9.5.txt",header=T,row.names=1,stringsAsFactors=F,sep='\t')

meta.use <- meta[colnames(mtx),]

color_use<-c("#FFA500","#F47D2B","#E6C2DC","#208A42","#89C75F","#89288F","#C06CAB","#90D5E4","#8A9FD1","#272E6A","#D3C035","#D51F26")
clusterV<-factor(meta.use$LineageAnnotations,levels=c("esophagus-1","esophagus-2","pharynx","respiratory-trachea","respiratory-lung","anaterior stomach","posterior stomach","liver-hepatocyte","hepatoblast","hepatopancreatic duct","dorsal pancreas","duodenum"))
colorV<-color_use[as.numeric(clusterV)]

pdf("result_plots/E95_umap_AuthorOrgan.pdf")
par(mar=c(8,4,6,10))
plot(umapmat,col=colorV,pch=16,cex=.4)
legend(11,6,legend=levels(clusterV),fill=color_use,bty="n",xpd=TRUE)
dev.off()

png("result_plots/E95_umap_AuthorOrgan_points.png")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=.6,axes=F,xlab="",ylab="")
dev.off()

pdf("result_plots/E95_umap_AuthorOrgan_axes.pdf")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=.4,type="n")
#legend(11,6,legend=levels(clusterV),fill=color_use,bty="n",xpd=TRUE)
dev.off()

pdf("result_plots/E95_umap_AuthorOrgan_legend.pdf")
par(mar=c(2,2,2,2))
plot.new()
legend("topleft",legend=levels(clusterV),fill=color_use,bty="n",xpd=TRUE)
dev.off()

#add cell type to metadata
meta <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95_scRNA/matrics/GSE136689_MetaFile_DE_and_SM_Cells_E8.5_E9.0_E9.5.txt",header=T,row.names=1,stringsAsFactors=F,sep='\t')
meta.use <- meta[colnames(mtx),]

E95 <- AddMetaData(E95, factor(meta.use$LineageAnnotations), col.name = "celltype")


#extract nomalized count matrix (by organ)
exp.mtx <- E95[["RNA"]]@data
avg.vct <- c()

for(ct in levels(Idents(E95))){
    cellname <- names(Idents(E95))[which(Idents(E95)==ct)]
    print(c(ct,length(cellname)))
    tmp.mtx <- exp.mtx[,cellname]
    tmp.vct <- apply(tmp.mtx,1,mean)
    avg.vct <- c(avg.vct,tmp.vct)
}

avg.dat <- as.data.frame(matrix(avg.vct,ncol=length(levels(Idents(E95)))))
rownames(avg.dat) <- rownames(exp.mtx)
colnames(avg.dat) <- levels(Idents(E95))

write.table(avg.dat,"E95_Expression_matrix_byLabel.txt",col.names=T,row.names=T,sep="\t",quote=F)

avg.vct <- c()
cellname <- names(Idents(E95))[which(Idents(E95)%in%c("anaterior stomach","posterior stomach"))]
print(length(cellname))
tmp.mtx <- exp.mtx[,cellname]
tmp.vct <- apply(tmp.mtx,1,mean)
avg.vct <- c(avg.vct,tmp.vct)
cellname <- names(Idents(E95))[which(Idents(E95)%in%c("duodenum"))]
print(length(cellname))
tmp.mtx <- exp.mtx[,cellname]
tmp.vct <- apply(tmp.mtx,1,mean)
avg.vct <- c(avg.vct,tmp.vct)
cellname <- names(Idents(E95))[which(Idents(E95)%in%c("dorsal pancreas"))]
print(length(cellname))
tmp.mtx <- exp.mtx[,cellname]
tmp.vct <- apply(tmp.mtx,1,mean)
avg.vct <- c(avg.vct,tmp.vct)
cellname <- names(Idents(E95))[which(Idents(E95)%in%c("respiratory-lung"))]
print(length(cellname))
tmp.mtx <- exp.mtx[,cellname]
tmp.vct <- apply(tmp.mtx,1,mean)
avg.vct <- c(avg.vct,tmp.vct)

avg.dat <- as.data.frame(matrix(avg.vct,ncol=4))
rownames(avg.dat) <- rownames(exp.mtx)
colnames(avg.dat) <- c("stomach","intestine","pancreas","lung")

write.table(avg.dat,"E95_Expression_matrix_byOrgan.txt",col.names=T,row.names=T,sep="\t",quote=F)


