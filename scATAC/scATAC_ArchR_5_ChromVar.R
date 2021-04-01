#ChromVar of Sox2/Cdx2 WT/KO ATAC peaks groups
library(GenomicRanges)

#read in peak sets
t1 <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Cdx2/Peak_group1.bed")
gr1 <- GRanges(
    seqnames = factor(t1$V1),
    ranges = IRanges(start = t1$V2,end = t1$V3)
)

t2 <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Cdx2/Peak_group2.bed")
gr2 <- GRanges(
    seqnames = factor(t2$V1),
    ranges = IRanges(start = t2$V2,end = t2$V3)
)

t3 <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Cdx2/Peak_group3.bed")
gr3 <- GRanges(
    seqnames = factor(t3$V1),
    ranges = IRanges(start = t3$V2,end = t3$V3)
)

t4 <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Sox2/Peak_group1.bed")
gr4 <- GRanges(
    seqnames = factor(t4$V1),
    ranges = IRanges(start = t4$V2,end = t4$V3)
)

t5 <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Sox2/Peak_group2.bed")
gr5 <- GRanges(
    seqnames = factor(t5$V1),
    ranges = IRanges(start = t5$V2,end = t5$V3)
)

t6 <- read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/Sox2Cdx2_peak/signal_Sox2/Peak_group3.bed")
gr6 <- GRanges(
    seqnames = factor(t6$V1),
    ranges = IRanges(start = t6$V2,end = t6$V3)
)

gr <- GRangesList("Cdx2_g1" = gr1, "Cdx2_g2" = gr2, "Cdx2_g3" = gr3, "Sox2_g1" = gr4, "Sox2_g2" = gr5, "Sox2_g3" = gr6)

proj <- addPeakAnnotations(
  ArchRProj = proj,
  regions = gr,
  name = "WTKOPeakSet",
  force = TRUE,
  logFile = createLogFile("addPeakAnnotations")
)

proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "WTKOPeakSet",
  matrixName = "WTKOPeakMatrix",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, plot = TRUE, name = "WTKOPeakMatrix")

markerMotif <- getFeatures(proj, useMatrix = "WTKOPeakMatrix")
markerMotif <- sort(grep("z:", markerMotif, value = TRUE))

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "WTKOPeakMatrix", 
    name = markerMotif, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

pdf("WTKOPeak_ChromVar.pdf")
for(i in 1:length(p)){
    print(p[[i]])
}
dev.off()
