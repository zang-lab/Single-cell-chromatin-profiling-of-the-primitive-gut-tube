#Identify Organ-specific genes (Marker Genes)
##build super cluster
#(1),(2,3,4),(5),(6,7),(8),(9),(10),(11),(12,13),(14),(15)
superClusters<-plyr::mapvalues(Cdata$newClusters_re1, from = paste0("C",c(1:7,"8a","8b",9:15)), to = c(rep("intestine",4),"pancreas",rep("liver",2),"stomach","lung","esophagus","unidentified","pharynx","colon","colon","unidentified","unidentified"))
proj <- addCellColData(
    ArchRProj = proj,
    data = superClusters,
    name = "superClusters",
    cell = rownames(Cdata),
    force = TRUE
)

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

#call marker gene by supercluster
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "superClusters",
    useGroups = c("intestine","pancreas","stomach","lung"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= log2(1.5)")

DEG<-markerList[[1]]
DEG$Clusters <- rep(names(markerList)[1],nrow(markerList[[1]]))
for(i in names(markerList)){
    print(paste0(i," ",nrow(markerList[[i]])))
    if(i==names(markerList)[1]){
        next
    }
    li<-markerList[[i]]
    li$Clusters <- rep(i,nrow(li))
    DEG<-rbind(DEG,li)
}

#filter out gene with GS<1 in over 50% cells
c="C15"
clusterDEG<-DEG$name[which(DEG$Clusters==c)]
print(c)
print(nrow(clusterDEG))
GeneScorecluster<-GeneScore[,which(Cdata$newClusters==c)]
GeneScorecluster<-GeneScorecluster[clusterDEG,]
remain_row=c()
for(i in 1:nrow(GeneScorecluster)){
    tmp<-GeneScorecluster[i,]
    if(length(tmp[which(tmp>=1)])/length(tmp)>0.5){
        remain_row=c(remain_row,i)
    }
}
print(length(remain_row))
DEGC15<-DEG[which(DEG$name%in%clusterDEG[remain_row] & DEG$Clusters==c),]

DEGother<-DEG[which(DEG$Clusters!="C15"),]
DEGf1<-rbind(DEGother,DEGC15)
DEGf1<-as.data.frame(DEGf1)

markers<-c(c("Sox2","Nkx2-1","Sftpc","Irx1","Foxp2","Ash2l","Elf5","Etv5"),c("Sox2","Hoxb4","Pdx1","Osr1"),
c("Pdx1","Hoxc5","Ptf1a","Neurog3","Nkx2-2","Nkx6-1","Tcf15"),c("Cdx1","Cdx2","Evx1","Pitx2","Tff3","Hoxa7","Hoxa9","Hoxb6","Hoxb7","Fabp1"),
c("Hoxa13","Hoxd13","Hoxb9"),c("Hhex","Msx1","Prox1","Creb3l3","Klf4","Ppy"))
print(table(markers%in%DEG$name))

#write out DEG
#orig.DEGC15<-DEG[which(DEG$Clusters=="C15"),]
#write.table(orig.DEGC15,"markerGene/markerGene_bysupercluster_C15orig.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(DEG,"markerGene/markerGene_bysupercluster_all.txt",row.names=F,col.names=T,quote=F,sep='\t')
for(i in names(markerList)){
    li<-DEG[which(DEG$Clusters==i),]
    write.table(li,paste0("markerGene/markerGene_bysupercluster_",i,".txt"),row.names=F,col.names=T,quote=F,sep='\t')
}

#marker gene heatmap
keys<-list(c("Sox2","Nkx2-1","Sftpc","Irx1","Foxp2","Ash2l","Elf5","Etv5"),c("Sox2","Hoxb4","Pdx1","Osr1"),
c("Pdx1","Hoxc5","Ptf1a","Neurog3","Nkx2-2","Nkx6-1","Tcf15"),c("Cdx1","Cdx2","Evx1","Pitx2","Tff3","Hoxa7","Hoxa9","Hoxb6","Hoxb7","Fabp1"),
c("Hoxa13","Hoxd13","Hoxb9"),c("Hhex","Msx1","Prox1","Creb3l3","Klf4","Ppy"))
names(keys)<-c("lung","stomach","pancreas","smallIntestine","largeIntestine","liver")

for(i in names(keys)){
    heatmapGS <- plotMarkerHeatmap(
        seMarker = markersGS,
        cutOff = "FDR <= 0.05 & Log2FC >= log2(1.5)",
        labelMarkers = keys[[i]],
        transpose = TRUE
    )
    if(class(heatmapGS)=="list"){
        heatmapGS <- plotMarkerHeatmap(
        	seMarker = markersGS,
        	cutOff = "FDR <= 0.05 & Log2FC >= log2(1.5)",
        	transpose = TRUE
    	)
    }
    pdf(paste0("markerGene/markerGene_bysupercluster_",i,"_heatmap.pdf"))
    ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
    #ComplexHeatmap::draw(Heatmap(heatmapGS[["mat"]], col=heatmapGS[["color"]], cluster_rows=F,cluster_columns=F, show_column_names=F,heatmap_legend_param=list(direction="horizontal",title=heatmapGS[["name"]],at=heatmapGS[["limits"]],legend_width=unit(3,"cm"))), heatmap_legend_side = "bot", annotation_legend_side = "bot")
    dev.off()
}


#all markers in the DEG heatmap
features<-getFeatures(proj,"GeneScoreMatrix")
markersGS<-markersGS[which(features%in%DEG$name),]
heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.05 & Log2FC >= log2(1.5)",
    labelMarkers = markers,
    transpose = TRUE
)

pdf(paste0("markerGene/markerGene_bysupercluster_all_heatmap.pdf"))
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#ComplexHeatmap::draw(Heatmap(heatmapGS[["mat"]], col=heatmapGS[["color"]], cluster_rows=F,cluster_columns=F, show_column_names=F,heatmap_legend_param=list(direction="horizontal",title=heatmapGS[["name"]],at=heatmapGS[["limits"]],legend_width=unit(3,"cm"))), heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#all markers heatmap
GeneScore_mk<-GeneScore[markers,]
GeneScore_mk_cluster<-GeneScore_mk[,1]
for(c in unique(Cdata$superClusters)){
    GeneScore_mk_cluster<-cbind(GeneScore_mk_cluster,apply(GeneScore_mk[,which(Cdata$superClusters==c)],1,mean))
}
GeneScore_mk_cluster<-GeneScore_mk_cluster[,-1]
colnames(GeneScore_mk_cluster)<-unique(Cdata$superClusters)
tmpV<-as.vector(as.matrix(GeneScore_mk_cluster))
plotM<-t(matrix(tmpV,ncol=ncol(GeneScore_mk_cluster)))
rownames(plotM)<-colnames(GeneScore_mk_cluster)
colnames(plotM)<-markers
plotM<-plotM[paste0("C",c(15,5,67,8,1,234,10,9,11,14,1213)),]
#log2 norm
plotM <- log2(t(t(plotM)/colSums(plotM)) * 10^4 + 1)
#scale row
plotM <- sweep(plotM - rowMeans(plotM), 1, matrixStats::rowSds(plotM), `/`)
limits=c(-2,2)
plotM[plotM > max(limits)] <- max(limits)
plotM[plotM < min(limits)] <- min(limits)

heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGSf1,
    cutOff = "FDR <= 0.01 & Log2FC >= 0.8",
    labelMarkers = c("Cdx2"),
    transpose = TRUE
)
pdf(paste0("markerGene/markerGeneMK_bysupercluster_all_heatmap.pdf"))
ComplexHeatmap::draw(Heatmap(plotM, col=heatmapGS[["color"]], cluster_rows=F,cluster_columns=T, show_column_names=T,heatmap_legend_param=list(direction="horizontal",title=heatmapGS[["name"]],at=heatmapGS[["limits"]],legend_width=unit(3,"cm"))), heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


#build pseudo-bulk replicates and call peaks
table(proj$superClusters)

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "superClusters",force=T)

#find macs2 path
pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "superClusters", 
    pathToMacs2 = pathToMacs2,
    additionalParams = "--nomodel --nolambda -q 0.01 --extsize 100",
    force=TRUE
)

#add peak matrix
proj <- addPeakMatrix(proj,force=TRUE)


#call Organ-specific peaks (Marker peaks)
table(proj$superClusters)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "superClusters",
    useGroups = c("intestine","pancreas","liver","stomach","lung","esophagus","pharynx"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markersPeaks
saveRDS(markersPeaks,"markerPeak/markerPeak_7organs/markersPeaks_bysuperClusters.rds")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= log2(2)")

difPeaks<-markerList[[1]]
#difPeaks$seqnames<-factor(difPeaks$seqnames,levels=c(paste0("chr",1:19),"chrX"))
#difPeaks<-difPeaks[order(difPeaks$seqnames,difPeaks$start),]
difPeaks$Clusters <- rep(names(markerList)[1],nrow(markerList[[1]]))
for(i in names(markerList)){
    print(paste0(i," ",nrow(markerList[[i]])))
    if(i==names(markerList)[1]){
        next
    }
    li<-markerList[[i]]
    #li$seqnames<-factor(li$seqnames,levels=c(paste0("chr",1:19),"chrX"))
	#li<-difPeaks[order(li$seqnames,li$start),]
    li$Clusters <- rep(i,nrow(li))
    difPeaks<-rbind(difPeaks,li)
}

#write out difPeaks
write.table(difPeaks,"markerPeak/markerPeak_7organs/markerPeak_bysuperCluster_all.txt",row.names=F,col.names=T,quote=F,sep='\t')
t1<-difPeaks[,c(1,3,4,8)]
write.table(t1,"markerPeak/markerPeak_7organs/markerPeak_bysuperCluster_all.bed",row.names=F,col.names=F,quote=F,sep='\t')
for(i in names(markerList)){
    li<-difPeaks[which(difPeaks$Clusters==i),]
    li<-li[,c(1,3,4)]
    write.table(li,paste0("markerPeak/markerPeak_7organs/markerPeak_bysuperCluster_",i,".bed"),row.names=F,col.names=F,quote=F,sep='\t')
}

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= log2(2)",
  transpose = TRUE,
)

pdf("markerPeak/markerPeak_7organs/markerPeak_bysuperCluster_heatmap.pdf")
ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


#

