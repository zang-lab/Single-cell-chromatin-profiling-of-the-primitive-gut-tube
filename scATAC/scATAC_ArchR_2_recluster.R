#get barcodes of cells in C8 and C9
library(ArchR)
set.seed(1)

proj <- loadArchRProject(path = "HemeTutorial")

#get cluster infomation
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)
cellname <- rownames(Cdata)[which(Cdata$newClusters%in%c("C8","C9"))]

cn1<-c()
cn2<-c()
for(c in cellname){
    if(substr(c,1,7)=="E9.5_r1"){
        cn1 <- c(cn1,substr(c,9,nchar(c)))
    }else{
        cn2 <- c(cn2,substr(c,9,nchar(c)))
    }
}

barcodes <- list(E9.5_r1=cn1,E9.5_r2=cn2)
saveRDS(barcodes,"reclustering_C8C9/C8C9_barcodes.rds")


#read in data and build the project
setwd("reclustering_C8C9")

## Setting default number of Parallel threads to 8
addArchRThreads(threads = 8)

inputFiles <- c("/scratch/hz9fq/toronto_project/scATAC/E95_gut_endoderm_scATAC_r1/outs/fragments.tsv.gz","/scratch/hz9fq/toronto_project/scATAC/E95_gut_endoderm_scATAC_r2/outs/fragments.tsv.gz")
names(inputFiles)<-c("E9.5_r1","E9.5_r2")

addArchRGenome("mm10")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  validBarcodes = barcodes, #select barcodes you want to use
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "C8C9reclustering",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

proj <- saveArchRProject(ArchRProj = proj)

proj <- loadArchRProject(path = "C8C9reclustering")



#reclustering
#dimensional reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# call clusters in this reduced dimension sub-space (seurat graph clustering)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

#correct batch effect
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
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

pdf("harmony_umap.pdf")
ggAlignPlots(p1, p2, type = "h")
dev.off()

#rename cluster if needed
library(plyr)
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)
newClusters <- mapvalues(Cdata$Clusters, from = c("C1","C2","C4","C5","C3","C6","C7","C8"), to = paste0("C",1:8))
proj <- addCellColData(
    ArchRProj = proj,
    data = newClusters,
    name = "newClusters",
    cell = rownames(Cdata),
    force = TRUE
)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "newClusters", embedding = "UMAP")

pdf("newClusters_umap.pdf")
ggAlignPlots(p1, p2, type = "h")
dev.off()


#GeneScore scatter plot
proj <- addImputeWeights(proj)

keys<-list(c("Sox2","Nkx2-1","Sftpc","Irx1","Foxp2","Ash2l","Elf5","Etv5"),c("Sox2","Hoxb4","Pdx1","Osr1"),
c("Pdx1","Hoxc5","Ptf1a","Neurog3","Nkx2-2","Nkx6-1","Tcf15"))
names(keys)<-c("lung","stomach","pancreas")

#umap using gene expression values from GeneScore
for(k in 1:length(keys)){
    markerGenes <- keys[[k]]
    print(markerGenes)
    
    p2 <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "GeneScoreMatrix", 
        name = markerGenes, 
        continuousSet = "horizonExtra",
        embedding = "UMAP",
        imputeWeights = getImputeWeights(proj)
    )
    
    #plot all the genes together
    p2c <- lapply(p2, function(x){
        x + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank()
        )
    })
    
    pdf(paste0("GeneScoreMatrix_",names(keys)[k],"_markerGene.pdf"))
    print(do.call(cowplot::plot_grid, c(list(ncol = 3,scale=0.8), p2c)))
    dev.off()
}


#line and dot plot
#get GeneScore matrix
GeneScore<-getMatrixFromProject(proj,"GeneScoreMatrix")
GeneScore<-GeneScore@assays@data@listData$GeneScoreMatrix
GeneScore<-as.matrix(GeneScore)
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)
rownames(GeneScore)<-features

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

for(i in 1:length(keys)){
    pdf(paste0("MarkerGene_Ryan_",names(keys)[i],"_E95_linedot.pdf"))
    par(mfrow=c(3,3),mar=c(6,4,6,2))
    for(K in keys[[i]]){
        if(K%in%rownames(GeneScore)==FALSE){
            next
        }
        # get cluster average expression
        exp<-as.vector(GeneScore[K,])
        thecluster<-Cdata$newClusters
        avgExp=c()
        for(c in paste0("C",c(1:8))){
            avgExp<-c(avgExp,mean(exp[which(thecluster==c)]))
        }
        # upper limit of confidence interval
        uppervalue<-Rmisc::CI(avgExp)[1]
        plot(avgExp,main=K,axes=T,type="o",xaxt="none",xlab="",ylab="GeneScore")
        abline(h=uppervalue, col="red")
        axis(1, at=c(1:length(avgExp)), labels=paste0("C",c(1:8)), cex.axis=0.65,las=2)
   }
    dev.off()
}


#use parentClusters C8, C9 to mark reclustering results
library(ArchR)
set.seed(1)

proj <- loadArchRProject(path = "HemeTutorial")

#get cluster infomation
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

C8cell <- rownames(Cdata[which(Cdata$newClusters=='C8'),])
C9cell <- rownames(Cdata[which(Cdata$newClusters=='C9'),])

setwd("reclustering_C8C9")

proj2 <- loadArchRProject("C8C9reclustering")
Cdata2 <- getCellColData(ArchRProj = proj2, select = NULL, drop = FALSE)

parentClusters <- c()
for(i in rownames(Cdata2)){
    if(i %in% C8cell){
        parentClusters <- c(parentClusters,"C8")
    }else if(i %in% C9cell){
        parentClusters <- c(parentClusters,"C9")
    }else{
        print("error")
    }
}

proj2 <- addCellColData(
  ArchRProj = proj2,
  data = parentClusters,
  name = "parentClusters",
  cell = rownames(Cdata2),
  force = FALSE
)

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "parentClusters", embedding = "UMAP")

pdf("parentClusters_umap.pdf")
print(p1)
dev.off()


#get reclustering results and map back to parentClusters
newClusters <- Cdata$newClusters
names(newClusters) <- rownames(Cdata)
childClusters <- Cdata2$newClusters
names(childClusters) <- rownames(Cdata2)

childClusters <- plyr::mapvalues(childClusters, from = paste0("C",1:8), to = paste0("C",c("8a",9,9,"8b",9,9,9,9)))

for(i in names(childClusters)){
    newClusters[i] <- childClusters[i]
}

proj <- addCellColData(
  ArchRProj = proj,
  data = newClusters,
  name = "newClusters_re1",
  cell = rownames(Cdata),
  force = TRUE
)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "newClusters_re1", embedding = "UMAP")
pdf("newClusters_re1_umap.pdf")
print(p2)
dev.off()

#umap by cluster
umapmat <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

library("RColorBrewer")
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
color_use<-stallion_pal
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)
colorV<-color_use[as.numeric(plyr::mapvalues(Cdata$newClusters_re1, from = paste0("C",c(1:7,"8a","8b",9:15)), to = 1:16))]

pdf("newClusters_re1_umap.pdf")
par(mar=c(6,4,6,8))
plot(umapmat,col=colorV,pch='.')
legend(14,11,legend=paste0("C",c(1:7,"8a","8b",9:15)),fill=color_use,bty="n",xpd=TRUE)
dev.off()


#plot GeneScore scatter using ggplot
library(ArchR)
set.seed(1)

umapmat <- getEmbedding(ArchRProj = proj2, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

Cdata2 <- getCellColData(ArchRProj = proj2, select = NULL, drop = FALSE)

#plot GeneScore
library(gridExtra)
library(grid)
library(ggplot2)

GeneScore<-getMatrixFromProject(proj2,"GeneScoreMatrix")
GeneScore<-GeneScore@assays@data@listData$GeneScoreMatrix
features<-getFeatures(
  ArchRProj = proj2,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)
rownames(GeneScore)<-features

imputedGeneScore <- imputeMatrix(mat = GeneScore, imputeWeights = getImputeWeights(proj2))
imputedGeneScore <- as.matrix(imputedGeneScore)

keys<-list(c("Sox2","Nkx2-1","Sftpc","Irx1","Foxp2","Ash2l","Elf5","Etv5"),c("Sox2","Hoxb4","Pdx1","Osr1"),
c("Pdx1","Hoxc5","Ptf1a","Neurog3","Nkx2-2","Nkx6-1","Tcf15"),c("Cdx1","Cdx2","Evx1","Pitx2","Tff3","Hoxa7","Hoxa9","Hoxb6","Hoxb7","Fabp1"),
c("Hoxa13","Hoxd13","Hoxb9"),c("Hhex","Msx1","Prox1","Creb3l3","Klf4","Ppy"),c("Pou5f1","Nanog","Klf4"))
names(keys)<-c("lung","stomach","pancreas","smallIntestine","largeIntestine","liver","stem")

horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")

plotdat<-as.data.frame(umapmat)
for(i in 1:length(keys)){
    p1<-list()
    for(K in keys[[i]]){
        if(K%in%rownames(imputedGeneScore)==FALSE){
            next
        }
        plotdat$exp<-as.vector(imputedGeneScore[K,])
        #transfer top genes expression to quantile(0.99)
        #uplimit<-quantile(plotdat$exp,0.95)
        #uplimit=3
        #plotdat$exp[which(plotdat$exp>=uplimit)]<-uplimit
        
        p1<-append(p1,list(ggplot(plotdat,aes(UMAP_1,UMAP_2,colour=exp))+
        scale_colour_gradientn(colours=horizon_pal)+
        geom_point(size=0.1,alpha=5/10)+
        ggtitle(K)+
        theme(aspect.ratio=1)))
    }
    m1 <- marrangeGrob(p1, ncol=2, nrow=2, as.table=F)
    ggsave(paste0("GeneScoreMatrix_",names(keys)[i],"_markerGene.pdf"), m1)
}


