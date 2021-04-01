#bart results scatter plot atac vs gene
setwd("D:/projects/toronto_project/updatefiles/update_210315")

t1 <- read.table("COAD_fc_bart_results.txt",header=T)
t2 <- read.table("COAD_gene_bart_results.txt",header=T)

rownames(t1) <- t1$TF
rownames(t2) <- t2$TF
t2 <- t2[rownames(t1),]

plotdat <- data.frame(ATAC_pvalue=t1$irwin_hall_pvalue,RNA_pvalue=t2$irwin_hall_pvalue)
rownames(plotdat) <- rownames(t1)

#-log10(pvalue)
plotdat[,1] <- -log10(plotdat[,1])
plotdat[,2] <- -log10(plotdat[,2])

#select tf
col.V <- rep("grey",nrow(plotdat))
col.V[which(plotdat[,1]>=2&plotdat[,2]>=2)] <- "red"

pdf("COAD_signature_gene_vs_peak.pdf")
par(mar=c(6,6,6,6))
plot(ATAC_pvalue ~ RNA_pvalue, plotdat,pch=16,cex=0.7,col=col.V,xlab="COAD signature genes -log10(p.value)",ylab="COAD signature peaks -log10(p.value)")
abline(h=2,col="grey",lty=2)
abline(v=2,col="grey",lty=2)
for(i in which(plotdat[,1]>=2&plotdat[,2]>=2)){
    text(plotdat[i,2],plotdat[i,1]+0.1,rownames(plotdat)[i],cex=0.7)
}
dev.off()


#bart results scatter plot pvalue vs rank
t1 <- read.table("COAD_fc_rev_bart_results.txt",header=T)
t1$logp <- -log10(t1$irwin_hall_pvalue)
t1$rank <- rank(-t1$logp,ties.method="first")
rownames(t1) <- t1$TF

#select tf
col.V <- rep("grey",nrow(t1))
col.V[which(t1$logp>=3)] <- "red"

pdf("COAD_fc_rev_bart_results.pdf")
par(mar=c(5,5,5,5))
plot(logp ~ rank,t1,pch=16,cex=0.7,col=col.V,ylab="-log10(p.value)",xlab="rank from high to low",axes=F)
coords<-par("usr")
axis(1,at=c(coords[1],1,nrow(t1),coords[2]),labels=c("",1,nrow(t1),""),cex.axis=0.7)
p <- pretty(coords[3:4])
p <- p[2:(length(p)-1)]
axis(2,at=c(coords[3],p,coords[4]),labels=c("",p,""),cex.axis=0.7)
axis(3,at=c(coords[1],coords[2]),labels=c("",""),cex.axis=0.7,col.ticks=NA)
axis(4,at=c(coords[3],coords[4]),labels=c("",""),cex.axis=0.7,col.ticks=NA)
abline(h=3,col="grey",lty=2)
for(i in 1:5){
    text(t1[i,9],t1[i,8]+0.1,rownames(t1)[i],cex=0.7)
}
dev.off()
