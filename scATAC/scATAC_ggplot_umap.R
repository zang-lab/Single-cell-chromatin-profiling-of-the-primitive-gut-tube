library(ArchR)
set.seed(1)
proj <- loadArchRProject(path = "HemeTutorial")

#get umap coordinates
umapmat <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

#get metadata
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

#umap by cluster
library("RColorBrewer")
#getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
#color_use<-getPalette(14)
stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
color_use<-stallion_pal[1:14]
colorV<-color_use[as.numeric(as.character(Cdata$predictedGroup_un))]

pdf("ArchR_IntegrationResult_bycluster_umap.pdf")
par(mar=c(6,4,6,8))
plot(umapmat,col=colorV,pch='.')
legend(14,7,legend=paste0("cluster ",c(1:5,8:14)),fill=color_use[-c(6,7)],bty="n",xpd=TRUE)
dev.off()

#plot pseudo expression (from Integrated RNA-seq)
library(gridExtra)
library(grid)
library(ggplot2)

PseudoExp<-getMatrixFromProject(proj,"GeneIntegrationMatrix")
PseudoExp<-PseudoExp@assays@data@listData$GeneIntegrationMatrix
PseudoExp<-as.matrix(PseudoExp)
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneIntegrationMatrix",
  select = NULL,
  ignoreCase = TRUE
)
rownames(PseudoExp)<-features

keys<-list(c("Sox2","Nkx2-1","Sftpc","Irx1","Foxp2","Ash2l","Elf5","Etv5"),c("Sox2","Hoxb4","Pdx1","Osr1"),
c("Pdx1","Hoxc5","Ptf1a","Neurog3","Nkx2-2","Nkx6-1","Tcf15"),c("Cdx1","Cdx2","Evx1","Pitx2","Tff3","Hoxa7","Hoxa9","Hoxb6","Hoxb7","Fabp1"),
c("Hoxa13","Hoxd13","Hoxb9"),c("Hhex","Msx1","Prox1","Creb3l3","Klf4","Ppy"),c("Pou5f1","Nanog","Klf4"))
names(keys)<-c("lung","stomach","pancreas","smallIntestine","largeIntestine","liver","stem")

horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")

keys<-c("Otx2","Pax8","Nkx2-5","Nkx2-3","Pyy","Isl1","Msx1","Msx2","Cdx4","Chrd","Nkx2-1","Irx2","Thbs1","Sox2","Pdx1","Osr1","Neurog3","Cdx1","Cdx2","Tff3","Hoxb9","Hhex","Creb3l3")

plotdat<-as.data.frame(umapmat)

for(i in 1:length(keys)){
    p1<-list()
    for(K in keys){
        if(K%in%rownames(PseudoExp)==FALSE){
            next
        }
        plotdat$exp<-as.vector(PseudoExp[K,])
        #transfer top genes expression to quantile(0.99)
        uplimit<-quantile(plotdat$exp,0.95)
        #uplimit=3
        plotdat$exp[which(plotdat$exp>=uplimit)]<-uplimit
        
        p1<-append(p1,list(ggplot(plotdat,aes(UMAP_1,UMAP_2,colour=exp))+
        scale_colour_gradientn(colours=horizon_pal)+
        geom_point(size=0.01,alpha=5/10)+
        ggtitle(K)+
        theme(aspect.ratio=1)))
    }
    m1 <- marrangeGrob(p1, ncol=2, nrow=2, as.table=F)
    ggsave("Additional_marker_gene_pseudoExp.pdf", m1)
    ggsave(paste0("IntegrationMatrix_",names(keys)[i],"_markerGene.pdf"), m1)
}



#plot GeneScore
library(gridExtra)
library(grid)
library(ggplot2)

GeneScore<-getMatrixFromProject(proj,"GeneScoreMatrix")
GSmat<-GeneScore@assays@data@listData$GeneScoreMatrix
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)
rownames(GSmat)<-features

imputedGSmat <- imputeMatrix(mat = GSmat,imputeWeights = getImputeWeights(proj))
imputedGSmat <- as.matrix(imputedGSmat)

#keys<-list(c("Sox2","Nkx2-1","Sftpc","Irx1","Foxp2","Ash2l","Elf5","Etv5"),c("Sox2","Hoxb4","Pdx1","Osr1"), c("Pdx1","Hoxc5","Ptf1a","Neurog3","Nkx2-2","Nkx6-1","Tcf15"),c("Cdx1","Cdx2","Evx1","Pitx2","Tff3","Hoxa7","Hoxa9","Hoxb6","Hoxb7","Fabp1"), c("Hoxa13","Hoxd13","Hoxb9"),c("Hhex","Msx1","Prox1","Creb3l3","Klf4","Ppy"),c("Pou5f1","Nanog","Klf4"))
#names(keys)<-c("lung","stomach","pancreas","smallIntestine","largeIntestine","liver","stem")

keys<-c("Otx2","Pax8","Nkx2-5","Nkx2-3","Pyy","Isl1","Msx1","Msx2","Cdx4","Chrd","Nkx2-1","Irx2","Thbs1","Sox2","Pdx1","Osr1","Neurog3","Cdx1","Cdx2","Tff3","Hoxb9","Hhex","Creb3l3")

keys<-c("Krt1","Krt12","Krt14","Trp63")

keys<-c("Hhex","Creb3l3","Pdx1","Neurog3","Cdx2","Tff3","Osr1","Sox2","Nkx2-1","Nkx6-1","Ptf1a")

horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")

plotdat<-as.data.frame(umapmat)


p1<-list()
for(K in keys){
    if(K%in%rownames(imputedGSmat)==FALSE){
        next
    }
    plotdat$exp<-as.vector(imputedGSmat[K,])
    #transfer top genes expression to quantile(0.99)
    uplimit<-quantile(plotdat$exp,0.95)
    #uplimit=3
    plotdat$exp[which(plotdat$exp>=uplimit)]<-uplimit
    
    p1<-append(p1,list(ggplot(plotdat,aes(UMAP_1,UMAP_2,colour=exp))+
    scale_colour_gradientn(colours=horizon_pal)+
    geom_point(size=0.1,alpha=5/10)+
    ggtitle(K)+
    theme(aspect.ratio=1)))
}
m1 <- marrangeGrob(p1, ncol=2, nrow=2, as.table=F)
ggsave("test.pdf", m1)

#save png format
keys=c("Nkx2-5", "Pax8", "Pyy", "Isl1", "Krt14", "Trp63", "Msx1", "Msx2")
keys=c("Cdx1")
for(K in keys){
    print(K)
    if(K%in%rownames(imputedGSmat)==FALSE){
        next
    }
    plotdat$exp<-as.vector(imputedGSmat[K,])
    #transfer top genes expression to quantile(0.99)
    uplimit<-quantile(plotdat$exp,0.95)
    #uplimit=3
    plotdat$exp[which(plotdat$exp>=uplimit)]<-uplimit
    
    p1<-list(ggplot(plotdat,aes(UMAP_1,UMAP_2,colour=exp))+
    scale_colour_gradientn(colours=horizon_pal)+
    geom_point(size=0.1,alpha=5/10)+
    ggtitle(K)+
    theme(aspect.ratio=1,axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
            ))
    png(paste0("plots/MarkerUMAPpng/",K,".png"))
    print(p1)
    dev.off()
}

#save axes and legend
print(K)
if(K%in%rownames(imputedGSmat)==FALSE){
    next
}
plotdat$exp<-as.vector(imputedGSmat[K,])
#transfer top genes expression to quantile(0.99)
uplimit<-quantile(plotdat$exp,0.95)
#uplimit=3
plotdat$exp[which(plotdat$exp>=uplimit)]<-uplimit

p1<-list(ggplot(plotdat,aes(UMAP_1,UMAP_2,colour=exp))+
scale_colour_gradientn(colours=horizon_pal)+
theme(aspect.ratio=1,
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
plot.background=element_blank()))

pdf("plots/MarkerUMAPpng/UMAP_axes.pdf")
print(p1)
dev.off()


#scRNA correlation boxplot
t0 <- read.table("/scratch/hz9fq/toronto_project/ALL_RNA_collection/FPKM_collection.txt",header=T,stringsAsFactors=F)
t1 <- data.frame(apply(t0[,c(7,8)],1,mean),apply(t0[,c(9,10)],1,mean),apply(t0[,c(11,12)],1,mean),apply(t0[,c(13,14)],1,mean),apply(t0[,c(15,16)],1,mean),apply(t0[,c(17,18)],1,mean),apply(t0[,c(19,20)],1,mean))
rownames(t1)<-t0$Gene.Name
colnames(t1)<-c("E135_intestine_RNA","E135_lung_RNA","E135_pancreas_RNA","E135_stomach_RNA","E165_colon_RNA","E165_hindstomach_RNA","E165_intestine_RNA")

bulkRNA <- t1
bulkRNA<-log2(bulkRNA+1)

#get GeneScore matrix
GeneScore<-getMatrixFromProject(ArchRProj = proj,useMatrix = "GeneScoreMatrix",verbose = TRUE)
GeneScoreMTX<-GeneScore@assays@data@listData$GeneScoreMatrix
GeneScoreMTX<-as.matrix(GeneScoreMTX)
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)
rownames(GeneScoreMTX)<-features


rownames(bulkRNA)<-toupper(rownames(bulkRNA))
rownames(GeneScoreMTX)<-toupper(rownames(GeneScoreMTX))

GeneUse<-intersect(toupper(rownames(bulkRNA)),rownames(GeneScoreMTX))
markerGene<-read.table("markerGene/markerGene_bycluster_all.txt",header=T,stringsAsFactors=F)

markerGene.use<-markerGene[1,]
for(i in unique(markerGene$Clusters)){
    tmptab<-markerGene[which(markerGene$Clusters==i),]
    if(nrow(tmptab)<100){
        markerGene.use <- rbind(markerGene.use,tmptab)
    }else{
        markerGene.use <- rbind(markerGene.use,tmptab[1:100,])
    }
    
}
markerGene.use <- markerGene.use[-1,]

GeneUse<-intersect(GeneUse,toupper(markerGene$name))
bulkRNA.use<-bulkRNA[GeneUse,]
mtx.use<-GeneScoreMTX[GeneUse,]

cormtx<-matrix(rep(0,ncol(mtx.use)*ncol(bulkRNA.use)),ncol=ncol(bulkRNA.use))
rownames(cormtx)<-colnames(mtx.use)
colnames(cormtx)<-colnames(bulkRNA.use)

for(i in 1:ncol(mtx.use)){
    for(j in 1:ncol(bulkRNA.use)){
        cormtx[i,j]<-cor(mtx.use[,i],bulkRNA.use[,j],method="pearson")
    }
    if(i%%1000 == 0){
        print(i)
    }
}

clusterV<-factor(proj$newClusters_re1,levels=paste0("C",c(1:7,"8a","8b",9:15)))

pdf("scATAC_pearsonCorWithbulkRNA_difgene_boxplot.pdf")
for(lv in levels(clusterV)){
    print(lv)

    plot.table<-cormtx[which(clusterV==lv),]
    plot.vector<-as.vector(as.matrix(plot.table))
    tag.vector=c()
    for(c in colnames(plot.table)){
        tag.vector<-c(tag.vector,rep(c,nrow(plot.table)))
    }
    dat<-data.frame(cor=plot.vector,tag=tag.vector)

    par(mar=c(8,4,8,4))
    boxplot(cor~tag,data=dat,outline = FALSE,las=2,cex.axis=0.7,xlab="",ylab="",
    ylim=c(0,0.4),main=paste0(lv," pearson correlation with bulk RNA"))
    means <- aggregate(cor~tag, dat, mean)
    points(1:nrow(means), means$cor, col = "red")
    text(1:nrow(means), means$cor+0.01, labels = round(means$cor,4),cex=0.5)    
}
dev.off()

#umap by cluster
library(ArchR)
set.seed(1)
proj <- loadArchRProject(path = "HemeTutorial")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

umapmat <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
pal<-stallion_pal
Cdata$newClusters <- factor(Cdata$newClusters,levels=paste0("C",c(1:15)))
colorV<-pal[as.numeric(Cdata$newClusters)]
png("plots/recluster/Whole_before_recluster_umap_points.png")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.25,axes=F)
dev.off()

pdf("plots/recluster/Whole_before_recluster_umap_legend.pdf")
par(mar=c(2,2,2,2))
plot.new()
legend("topleft",legend=levels(as.factor(Cdata$newClusters)),fill=pal,bty="n",xpd=TRUE)
dev.off()

Cdata$newClusters_re1 <- factor(Cdata$newClusters_re1,levels=c(paste0("C",c(1:7)),"C8a","C8b",paste0("C",c(9:15))))
pal<-stallion_pal[c(1:7,16,17,9:15)]
colorV<-pal[as.numeric(Cdata$newClusters_re1)]
png("plots/recluster/Whole_after_recluster_umap_points.png")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.25,axes=F)
dev.off()

pdf("plots/recluster/Whole_after_recluster_umap_legend.pdf")
par(mar=c(2,2,2,2))
plot.new()
legend("topleft",legend=levels(as.factor(Cdata$newClusters_re1)),fill=pal,bty="n",xpd=TRUE)
dev.off()

#umap by supercluster/organ
library(ArchR)
set.seed(1)
proj <- loadArchRProject(path = "HemeTutorial")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

umapmat <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

#umap by cluster
library("RColorBrewer")
#getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
#color_use<-getPalette(14)
Cdata$superClusters <- factor(Cdata$superClusters,levels=c("colon","pharynx","esophagus","lung","stomach","liver","pancreas","intestine","unidentified"))
pal<-c("#FFA07A","#E6C2DC","#F47D2B","#89C75F","#C06CAB","#90D5E4","#D3C035","#D51F26","#C0C0C0")
colorV<-pal[as.numeric(Cdata$superClusters)]

pdf("plots/ArchR_superClusters_umap.pdf")
par(mar=c(6,4,8,10))
plot(umapmat,col=colorV,pch=16,cex=0.25)
legend(13.5,7,legend=levels(as.factor(Cdata$superClusters)),fill=pal,bty="n",xpd=TRUE)
dev.off()

png("plots/ArchR_superClusters_umap_points.png")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.25,axes=F)
dev.off()

pdf("plots/ArchR_superClusters_umap_axes.pdf")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.25,type="n")
dev.off()

pdf("plots/ArchR_superClusters_umap_legend.pdf")
par(mar=c(2,2,2,2))
plot.new()
legend("topleft",legend=levels(as.factor(Cdata$superClusters)),fill=pal,bty="n",xpd=TRUE)
dev.off()

#C8 C9 umap by cluster
library(ArchR)
set.seed(1)
proj <- loadArchRProject(path = "C8C9reclustering")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

umapmat <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
colnames(umapmat) <- c("UMAP_1","UMAP_2")

stallion_pal=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")
pal<-stallion_pal
Cdata$newClusters <- factor(Cdata$newClusters,levels=paste0("C",c(1:8)))
colorV<-pal[as.numeric(Cdata$newClusters)]
png("../plots/recluster/C8C9_recluster_umap_points.png")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.5,axes=F)
dev.off()

pdf("../plots/recluster/C8C9_recluster_umap_legend.pdf")
par(mar=c(2,2,2,2))
plot.new()
legend("topleft",legend=levels(as.factor(Cdata$newClusters)),fill=pal,bty="n",xpd=TRUE)
dev.off()

pdf("../plots/recluster/C8C9_recluster_umap_axes.pdf")
par(mar=c(2,2,2,2))
plot(umapmat,col=colorV,pch=16,cex=0.5,type="n")
dev.off()

#genescore umap for C8C9
library(gridExtra)
library(grid)
library(ggplot2)

GeneScore<-getMatrixFromProject(proj,"GeneScoreMatrix")
GSmat<-GeneScore@assays@data@listData$GeneScoreMatrix
features<-getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)
rownames(GSmat)<-features

imputedGSmat <- imputeMatrix(mat = GSmat,imputeWeights = getImputeWeights(proj))
imputedGSmat <- as.matrix(imputedGSmat)

horizon_pal=c("#000075","#2E00FF","#9408F7","#C729D6","#FA4AB5","#FF6A95","#FF8B74","#FFAC53","#FFCD32","#FFFF60")
plotdat <- umapmat

keys=c("Nkx2-1", "Pdx1")
for(K in keys){
    print(K)
    if(K%in%rownames(imputedGSmat)==FALSE){
        next
    }
    plotdat$exp<-as.vector(imputedGSmat[K,])
    #transfer top genes expression to quantile(0.99)
    uplimit<-quantile(plotdat$exp,0.95)
    #uplimit=3
    plotdat$exp[which(plotdat$exp>=uplimit)]<-uplimit
    
    p1<-list(ggplot(plotdat,aes(UMAP_1,UMAP_2,colour=exp))+
    scale_colour_gradientn(colours=horizon_pal)+
    geom_point(size=0.75,alpha=5/10)+
    ggtitle(K)+
    theme(aspect.ratio=1,axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
            ))
    png(paste0("../plots/recluster/C8C9_",K,".png"))
    print(p1)
    dev.off()
}


