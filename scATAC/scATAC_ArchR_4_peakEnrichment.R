#custom enrichment
#read in organ-specific peaks (marker peaks) defined before
library(ArchR)
set.seed(1)
proj <- loadArchRProject(path = "HemeTutorial")
markersPeaks<-readRDS("markerPeak/markerPeak_7organs/markersPeaks_bysuperClusters.rds")

#read in E13.5 organ-specific peaks
peakdir="/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/E13bulkATACpeak/signal/difPeaks_kmeans/"
customPeaks <- c(
	lung_specificPeak = paste0(peakdir,"lung_difPeak_kmeans.bed"),
	pancreas_specificPeak = paste0(peakdir,"pancreas_difPeak_kmeans.bed"),
    intestine_specificPeak = paste0(peakdir,"intestine_difPeak_kmeans.bed"),
    stomach_specificPeak = paste0(peakdir,"stomach_difPeak_kmeans.bed")
)

proj <- addPeakAnnotations(ArchRProj = proj, regions = customPeaks, name = "custom_peak",force = TRUE)

enrichRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "custom_peak",
    cutOff = "FDR <= 0.1 & Log2FC >= log2(2)"
)

enrichRegions

mlog10Padj<-enrichRegions@assays@data@listData$mlog10Padj

heatmapEM <- plotEnrichHeatmap(enrichRegions, n = 10, cutOff = 0, transpose = TRUE)
pdf("markerPeak/markerPeak_7organs/customePeak4_heatmap_bysuperCluster.pdf")
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#limits=c(0,500)
#mlog10Padj[mlog10Padj > max(limits)] <- max(limits)
#mlog10Padj[mlog10Padj < min(limits)] <- min(limits)

#normalize for visualization
plotdat<-t(mlog10Padj)
plotdat <- plotdat[c("intestine","pancreas","stomach","lung","liver","pharynx","esophagus"),]
plotdat.norm <- plotdat
for(i in 1:ncol(plotdat.norm)){
    plotdat.norm[,i] <- plotdat.norm[,i]/max(plotdat.norm[,i])
}

library(pheatmap)
pdf("markerPeak/markerPeak_7organs/customePeak4_heatmap2_bysuperCluster.pdf")
pheatmap(plotdat.norm,breaks=seq(0,1,length.out=101),color=colorRampPalette(colors = c("#E6E7E8","#3A97FF","#8816A7","black"))(100),display_numbers = round(plotdat,2),cluster_rows=F,cluster_cols=F)
dev.off()
