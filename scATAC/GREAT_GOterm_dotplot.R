setwd("/scratch/hz9fq/toronto_project/scATAC/ArchR_workspace/markerPeak/markerPeak_7organs/GREAT/")
files=c("greatExportAll_intestine.tsv","greatExportAll_pancreas.tsv","greatExportAll_liver.tsv","greatExportAll_stomach.tsv","greatExportAll_lung.tsv","greatExportAll_esophagus.tsv","greatExportAll_pharynx.tsv")
names=c("intestine","pancreas","liver","stomach","lung","esophagus","pharynx")

great_res <- list()
term_use <- c()
for(i in 1:length(files)){
    t1 <- read.table(files[i],sep="\t")
    t2 <- t1[which(t1$V1=="GO Biological Process"),]
    t2 <- t2[order(t2$V7),c(1,2,3,7,21)]
    t2.use <- t2[1:5,]
    #print(nrow(t1.use))
    term_use <- unique(c(term_use,t2.use$V3))
    great_res[[names[i]]] <- t2
}

df1 <- data.frame(Term=0,ratio=0,FDR=0,organ=0)
for(i in 1:length(files)){
    tmp <- great_res[[i]]
    tmp2 <- tmp[which(tmp$V3%in%term_use),c("V3","V21","V7")]
    tmp2$organ <- names[i]
    colnames(tmp2) <- colnames(df1)
    df1 <- rbind(df1,tmp2)
}

df1 <- df1[-1,]

df1$Term <- factor(df1$Term,levels=rev(term_use))
df1$organ <- factor(df1$organ,levels=names)
df1$mlogFDR <- -log10(df1$FDR)
df1$ratio <- df1$ratio

library(ggplot2)
solarExtra_pal<-c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D")

#get FDR rank
organ_names <- names(table(df1$organ))
rank_V <- c()
for(o in organ_names){
    rank_V <- c(rank_V,rank(-df1[which(df1$organ==o),5]))
}
df1$rank <- rank_V

p1 <- ggplot(df1) + aes(x=organ,y=Term,colour=rank,size=ratio) + scale_colour_gradientn(colours=rev(solarExtra_pal),oob=scales::squish,name="rank by FDR") + scale_size(range = c(3, 12)) + geom_point() + theme(axis.text.y = element_text(size = 12))
ggsave("test.pdf",height=12,width=12,p1)




