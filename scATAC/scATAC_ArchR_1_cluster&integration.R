#step up input files
library(ArchR)
set.seed(1)

## Setting default number of Parallel threads to 8
addArchRThreads(threads = 8)

inputFiles <- c("/scratch/hz9fq/toronto_project/scATAC/E95_gut_endoderm_scATAC_r1/outs/fragments.tsv.gz","/scratch/hz9fq/toronto_project/scATAC/E95_gut_endoderm_scATAC_r2/outs/fragments.tsv.gz")
names(inputFiles)<-c("E9.5_r1","E9.5_r2")
#inputFiles <- getTutorialData("Hematopoiesis")

addArchRGenome("mm10")


#creating arrow files
#read in and calculate GeneScore at the same time
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles


#inferring doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)


#creating an ArchRproject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

#save and load ArchRporject
#To easily save an ArchRProject for later use or for sharing with collaborators, we use the saveArchRProject() function. This copies the current ArchRProject object and all of the Arrow files to a specified directory. If we don’t specify an output directory (as below), saveArchRProject() uses the output directory that we specified upon creation of our ArchRProject. In this case that is the folder “HemeTutorial”.
proj <- saveArchRProject(ArchRProj = proj)
#When we are ready to load this saved ArchRProject we use the loadArchRProject() object and provide the path to the folder containing the saved ArchRProject object.
proj <- loadArchRProject(path = "HemeTutorial")


#LSI dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)

# call clusters in this reduced dimension sub-space (seurat graph clustering)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")


#visualize in a 2d umap
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)

#color by sample
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#color by cluster
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

pdf("umap.pdf")
ggAlignPlots(p1, p2, type = "h")
dev.off()


#batch effect correction (harmony)
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", force = TRUE)

#color by sample
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#color by cluster
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "newClusters", embedding = "UMAP")

pdf("harmony_umap.pdf")
ggAlignPlots(p1, p2, type = "h")
dev.off()

#rename cluster if needed
library(plyr)
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)
newClusters <- mapvalues(Cdata$Clusters, from = c("C6","C8","C9","C7","C2","C3","C4","C5","C12","C10","C11","C14","C15","C13","C1"), to = paste0("C",1:15))
proj <- addCellColData(
    ArchRProj = proj,
    data = newClusters,
    name = "newClusters",
    cell = rownames(Cdata),
    force = TRUE
)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "newClusters", embedding = "UMAP")

pdf("newCluster_umap.pdf")
ggAlignPlots(p1, p2, type = "h")
dev.off()


#(Integration) define cluster identifed with scRNA-seq
library(ArchR)
set.seed(1)
proj <- loadArchRProject(path = "HemeTutorial")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

#ArchR accepts unmodified Seurat objects as input to the integration workflow
E95<-readRDS("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E95_scRNA/E95.rds")
set.seed(1)

# add the linked gene expression data to each of the Arrow files
# if you don't perform constrained integration, run this step directly
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = E95,
    addToArrow = FALSE,
    force= FALSE,
    groupRNA = "celltype",
    groupATAC = "superClusters",
    nameCell = "predictedCell_un_E95",
    nameGroup = "predictedGroup_un_E95",
    nameScore = "predictedScore_un_E95"
)

stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
#adjust the order because the order of esophagus clusters changed
pal<-stallion_pal[c(1:3,5,4,6,8:12)]

p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_un_E95", 
    pal = pal
)

pdf("plots/predictedGroup_un_E95.pdf")
p1
dev.off()

#basic plot umap
umapmat <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)
#umap by cluster
library("RColorBrewer")
#getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
#color_use<-getPalette(14)

pal<-c("#D8A767","#F47D2B","#E6C2DC","#208A42","#89C75F","#89288F","#C06CAB","#90D5E4","#8A9FD1","#FEE500","#D51F26")
Cdata$predictedGroup_un_E95 <- factor(Cdata$predictedGroup_un_E95,levels=c("esophagus-1","esophagus-2","pharynx","respiratory-trachea","respiratory-lung","anaterior stomach","posterior stomach","liver-hepatocyte","hepatoblast","dorsal pancreas","duodenum"))

colorV<-pal[as.numeric(as.factor(Cdata$predictedGroup_un_E95))]

pdf("plots/ArchR_IntegrationResult_E95byOrgan_umap.pdf")
par(mar=c(6,4,8,10))
plot(umapmat,col=colorV,pch=16,cex=0.25)
legend(13.5,7,legend=levels(as.factor(Cdata$predictedGroup_un_E95)),fill=pal,bty="n",xpd=TRUE)
dev.off()

#basic plot umap (exclude colon and unidentified)
anno.V <- as.character(Cdata$predictedGroup_un_E95)
anno.V[which(Cdata$superClusters=="colon" | Cdata$superClusters=="unidentified")] <- "colon & unidentified"
anno.V <- factor(anno.V,levels=c("esophagus-1","esophagus-2","pharynx","respiratory-trachea","respiratory-lung","anaterior stomach","posterior stomach","liver-hepatocyte","hepatoblast","dorsal pancreas","duodenum","colon & unidentified"))
pal<-c("#FFA500","#F47D2B","#E6C2DC","#208A42","#89C75F","#89288F","#C06CAB","#90D5E4","#8A9FD1","#D3C035","#D51F26","#C0C0C0")
colorV<-pal[as.numeric(as.factor(anno.V))]
pdf("plots/ArchR_IntegrationResult_E95byOrgan2_umap.pdf")
par(mar=c(6,4,8,10))
plot(umapmat,col=colorV,pch=16,cex=0.25)
legend(13.5,7,legend=levels(as.factor(anno.V)),fill=pal,bty="n",xpd=TRUE)
dev.off()

png("plots/ArchR_IntegrationResult_E95byOrgan2_umap_points.png")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.25,axes=F)
dev.off()

pdf("plots/ArchR_IntegrationResult_E95byOrgan2_umap_axes.pdf")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.25,type="n")
dev.off()

pdf("plots/ArchR_IntegrationResult_E95byOrgan2_umap_legend.pdf")
par(mar=c(2,2,2,2))
plot.new()
legend("topleft",legend=levels(as.factor(anno.V)),fill=pal,bty="n",xpd=TRUE)
dev.off()

#get metadata
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

proj <- addImputeWeights(proj)

#predctedScore scatter plot
p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = "predictedScore_un",
    pal = ArchRPalettes$coolwarm,
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

pdf("predictedScore_un_scatter.pdf")
print(p1)
dev.off()

#use cluster to integrate
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = E95,
    addToArrow = FALSE,
    force= FALSE,
    groupRNA = "seurat_clusters",
    nameCell = "predictedCell_un_E95Cluster",
    nameGroup = "predictedGroup_un_E95Cluster",
    nameScore = "predictedScore_un_E95Cluster"
)

stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
#adjust the order because the order of esophagus clusters changed
pal<-stallion_pal[c(1:7)]

p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_un_E95Cluster", 
    pal = pal
)

pdf("plots/predictedGroup_un_E95Cluster.pdf")
p1
dev.off()


#confusion matrix heatmap
cM <- confusionMatrix(proj$newClusters_re1, proj$predictedGroup_un_E95)
cM <- as.matrix(cM)
cM <- cM[paste0("C",c(1:7,"8a","8b",9:15)),]
cM <- cM[,c("duodenum","dorsal pancreas","liver-hepatocyte","hepatoblast","posterior stomach","anaterior stomach","respiratory-lung","esophagus-1","esophagus-2","pharynx","respiratory-trachea")]

library("pheatmap")
library("RColorBrewer")
pheatmap(cM,cellwidth=25,cellheight=25,display_numbers=T,cluster_rows=F,cluster_cols=F,number_format = "%d")

cM_norm <- cM
#normalize by colomn
for(i in 1:ncol(cM_norm)){
    cM_norm[,i] <- cM[,i]/max(cM[,i])
}

pdf("confusion_matrix_byCluster_normByCol.pdf")
pheatmap(cM_norm,cellwidth=25,cellheight=25,display_numbers=cM,cluster_rows=F,cluster_cols=F,number_format = "%d")
dev.off()

