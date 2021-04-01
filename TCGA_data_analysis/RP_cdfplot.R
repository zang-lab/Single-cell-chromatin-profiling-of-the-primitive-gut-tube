#Rp cdf
set.seed(1)

setwd("/scratch/hz9fq/toronto_project/TCGA/count")
#/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/RP/Sox2motif_with_20k_RP.bed
t0<-read.table("/sfs/lustre/bahamut/scratch/hz9fq/toronto_project/TCGA/ATAC-seq/Sox2motif_with_pan_cancer_RP.bed",header=F)
t1<-read.table("COAD_SOX2high_signature_hg38.txt",header=F)

t0$V7 <- as.numeric(t0$V6%in%t1$V1)
t2 <- t0[,c(5,7)]
t2$V7[which(t2$V7==0)] <- "other"
t2$V7[which(t2$V7==1)] <- "SOX2 high signature"
t2$V7 <- factor(t2$V7,levels=c("SOX2 high signature","other"))
colnames(t2) <- c("rp","group")

rp_high <- t2$rp[which(t2$group=="SOX2 high signature")]
rp_other <- t2$rp[which(t2$group=="other")]

ks.res <-ks.test(rp_high,rp_other,alternative="less")

pdf("COAD_rp_cdf.pdf")
par(mar=c(10,6,6,10))
ecdf1 <- ecdf(rp_high)
ecdf2 <- ecdf(rp_other)
plot(ecdf1, verticals=TRUE, do.points=FALSE,col="orange",main=NA,xlab="Regulatory Potential",lwd=1.5,xaxs="i",yaxs="i",xlim=c(0,max(t2$rp)),col.01line="black")
lines(ecdf2, verticals=TRUE, do.points=FALSE,col='blue',lwd=1.5,xaxs="i",yaxs="i",col.01line="black")
coords<-par("usr")
P <- ks.res$p.value
P_final <- "1.280x10^-4"
legend(coords[2]+0.2,(coords[3]+coords[4])/2,c("SOX2 high signature genes","other genes"),col=c("orange","blue"),lty=1,lwd=1.5,cex=0.6,xpd=T)
legend((coords[1]+coords[2])/2-1.5,(coords[3]+coords[4])/2,legend=paste0("P-value = ",P_final),bty="n",cex=0.8)
dev.off()
