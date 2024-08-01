library(stats)
library(gplots)
library(ggplot2)
#library(rgl)
library(multcomp)
library(car)
library(multtest)
library(edgeR)
library(effects)
library(limma)
library(edgeR)
library(pheatmap)
library(stringr)

#------------------------------------------------------------------------------
# Import Data
#------------------------------------------------------------------------------

lll <- read.table("data/F_data_phospho.txt",header=T,row.names=1)
# lll <- lll[,c(4,1,5,2,6,3)]
col.2.use <- c(3,3,3,3,3,2,2,2,2,2,4,4,4,4,4)

#------------------------------------------------------------------------------
# Perform pedestal and transformation
#------------------------------------------------------------------------------

temp1 <- lll
temp1 <- log(temp1+2,2)
sample.data <- temp1

#------------------------------------------------------------------------------
# Pre-Normalization Box Plot
#------------------------------------------------------------------------------

WorkingList <- vector("list", dim(sample.data)[2])
for ( i in 1:dim(sample.data)[2]) {
  WorkingList[[i]] <- as.numeric(sample.data[,i])
}
pdf("phospho/pre_norm_box_plot.pdf")
boxplot(WorkingList,col=col.2.use,outline=FALSE,xlab="Sample Index",ylab="Log2(Detection Value+2)")
# legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(15,15))
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red", "blue"),pch=c(15,15))
dev.off()

#------------------------------------------------------------------------------
# Perform cross sample normalization
#------------------------------------------------------------------------------

running.normalized <- normalizeBetweenArrays(sample.data,"cyclicloess")

#------------------------------------------------------------------------------
# Post normaliation box plot
#------------------------------------------------------------------------------

WorkingList <- vector("list", dim(running.normalized)[2])
for ( i in 1:dim(running.normalized)[2]) {
  WorkingList[[i]] <- as.numeric(running.normalized[,i])
}

pdf("phospho/post_norm_box_plot.pdf")
boxplot(WorkingList,col=col.2.use,outline=FALSE,xlab="Sample Index",ylab="CyclicLoess(Log2(Detection Value+2))")
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red", "blue"),pch=c(15,15))
dev.off()

#------------------------------------------------------------------------------
# Post normalization PCA plot
#------------------------------------------------------------------------------

pca.results <- princomp(running.normalized,cor=F)
pca.loadings <- loadings(pca.results)
phosphovar <- (pca.results$sdev^2)
variancePer <- round(phosphovar/sum(phosphovar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")

pdf("phospho/post_norm_PCA.pdf")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
dev.off()

#-------------------------------------------------------------------------------
# Perform Noise Modeling
#-------------------------------------------------------------------------------

group1 <- running.normalized[,col.2.use==3]
g1.mean <- apply(group1,1,mean)
g1.sd <- apply(group1,1,sd)
g1.cv <- (g1.sd/g1.mean)*100
lowess.g1 <- lowess(g1.cv~g1.mean,f=1/32)

group2 <- running.normalized[,col.2.use==2]
g2.mean <- apply(group2,1,mean)
g2.sd <- apply(group2,1,sd)
g2.cv <- (g2.sd/g2.mean)*100
lowess.g2 <- lowess(g2.cv~g2.mean,f=1/32)

group3 <- running.normalized[,col.2.use==4]
g3.mean <- apply(group3,1,mean)
g3.sd <- apply(group3,1,sd)
g3.cv <- (g3.sd/g3.mean)*100
lowess.g3 <- lowess(g3.cv~g3.mean,f=1/32)

#-------------------------------------------------------------------------------
# Generate Noise Model Plot
#-------------------------------------------------------------------------------

xrange <- range(c(1:14))
yrange <- range(c(0:50))
plot(xrange,yrange,type='n',xlab="Mean(CyclicLoess(Log2(Detection Value+2)))",ylab="Coefficient of Variation")
points(g1.mean,g1.cv,pch=21,col=8)
points(g2.mean,g2.cv,pch=21,col=8)
points(g3.mean,g3.cv,pch=21,col=8)
lines(lowess.g1,col=3,lwd=2)
lines(lowess.g2,col=2,lwd=2)
lines(lowess.g3,col=4,lwd=2)
abline(h=0)
abline(v=4.5,lty=2,lwd=2)
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
pdf("phospho/Noise_model_plot.pdf")
plot(xrange,yrange,type='n',xlab="Mean(CyclicLoess(Log2(Detection Value+2)))",ylab="Coefficient of Variation")
points(g1.mean,g1.cv,pch=21,col=8)
points(g2.mean,g2.cv,pch=21,col=8)
points(g3.mean,g3.cv,pch=21,col=8)
lines(lowess.g1,col=3,lwd=2)
lines(lowess.g2,col=2,lwd=2)
lines(lowess.g3,col=4,lwd=2)
abline(h=0)
abline(v=4.5,lty=2,lwd=2)
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
dev.off()

#-------------------------------------------------------------------------------
# Remove Noise & Floor
#-------------------------------------------------------------------------------

noise.pick <- 4.5
dim(running.normalized)
rrr <- running.normalized
rrr <- ifelse(rrr<noise.pick,noise.pick,rrr)
temp <- apply(rrr,1,mean)>noise.pick
Final.working.data <- rrr[temp,]
dim(Final.working.data)
write.table(Final.working.data,"phospho/Final_working_data_phospho.txt",sep="\t")

#------------------------------------------------------------------------------
# Post noise PCA plot
#------------------------------------------------------------------------------

pca.results <- princomp(Final.working.data,cor=F)
pca.loadings <- loadings(pca.results)
phosphovar <- (pca.results$sdev^2)
variancePer <- round(phosphovar/sum(phosphovar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
pdf("phospho/post_noise_PCA.pdf")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
dev.off()


#-------------------------------------------------------------------------------
# Significance Testing
#-------------------------------------------------------------------------------

#R17 vs GFP
running.pvs <- NULL
for (i in 1:dim(Final.working.data)[1]) {
  # print(i)
  temp1 <- as.numeric(Final.working.data[i,c(6,7,8,9,10)])
  temp2 <- as.numeric(Final.working.data[i,c(11,12,13,14,15)])
  temp3 <- t.test(temp1,temp2,paired = T)
  temp4 <- temp3$p.value
  temp5 <- temp3$estimate
  temp6 <- 2^temp5
  temp7 <- ifelse(temp6<1,-1/temp6,temp6)
  temp7 <- temp7*(-1) # reverse difference to be CD vs Ctrl
  running.pvs <- rbind(running.pvs,c(temp7,abs(temp7),temp4))
}
dimnames(running.pvs)[[1]] <- dimnames(Final.working.data)[[1]]

#-------------------------------------------------------------------------------
# FDR MCC (alpha=0.05)
#-------------------------------------------------------------------------------

temp1 <- mt.rawp2adjp(running.pvs[,3], proc="BH")
temp2 <- temp1$adj[order(temp1$index),]
temp2 <- as.data.frame(temp2)
running.pvs <- cbind(running.pvs,temp2[,2])

#-------------------------------------------------------------------------------
# Significance Testing Using Limma
#-------------------------------------------------------------------------------
# dir.create('phospho/limma/')
design <- read.csv('data/design.csv')
fit <- lmFit(Final.working.data, design)
contrast.matrix <- makeContrasts(GFP-R17, GFP-Fing, R17-Fing, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

GFP_vs_r17 <- topTable(fit2, coef=1, adjust="BH", number = length(rownames(lll)))
write.csv(GFP_vs_r17, 'phospho/limma/GFP_vs_r17_phospho.csv')

GFP_vs_fing <- topTable(fit2, coef=2, adjust="BH", number = length(rownames(lll)))
write.csv(GFP_vs_fing, 'phospho/limma/GFP_vs_fing_phospho.csv')

r17_vs_fing <- topTable(fit2, coef=3, adjust="BH", number = length(rownames(lll)))
write.csv(r17_vs_fing, 'phospho/limma/r17_vs_fing_phospho.csv')

results <- decideTests(fit2)
pdf('phospho/limma/limma_venn.pdf')
vennDiagram(results)
dev.off()

sig_norm <- Final.working.data[fit2$F.p.value < 0.0001,]
sig_raw <- lll[fit2$F.p.value < 0.0001,]


#-------------------------------------------------------------------------------
# Filter select
#-------------------------------------------------------------------------------

temp1 <- running.pvs[,2]>1.25
temp2 <- running.pvs[,4]<0.05
temp3 <- (temp1+temp2)==2
passing.mcc.frags <- names(temp3[temp3==TRUE])
length(passing.mcc.frags)
running.pvs <- cbind(running.pvs,as.logical(temp3))
dimnames(running.pvs)[[2]] <- c("R17_vs_Fing_Difference","R17_vs_Fing_Difference","Paired_Ttest_uncorrected_P","Paired_Ttest_Corrected_P","Pass_Selection_0_No_1_Yes")
write.table(running.pvs, "phospho/R17_vs_Fing_differential_testing_results.txt", sep = "\t")

coded.keepers <- dimnames(running.pvs)[[1]][running.pvs[,5]==1]
coded.keepers <- na.omit(coded.keepers)
grab.keeper.info <- running.pvs[coded.keepers,]
write.table(grab.keeper.info, "phospho/R17_vs_Fing_selected_differential.txt", sep = "\t")

#-------------------------------------------------------------------------------
# Post Selection Visuals - PCA
#-------------------------------------------------------------------------------

# pca.results <- princomp(sig,cor=F)
pca.results <- princomp(Final.working.data[coded.keepers,],cor=F)


pca.loadings <- loadings(pca.results)
phosphovar <- (pca.results$sdev^2)
variancePer <- round(phosphovar/sum(phosphovar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
pdf("phospho/post_selection_PCA.pdf")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
dev.off()

# 
# par(mfrow=c(1,2))
# x.range <- range(pca.loadings[,1])
# y.range <- range(pca.loadings[,2])
# xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
# ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
# plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
# points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
# points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
# legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
# abline(h=0,lty=2)
# x.range <- range(pca.loadings[,2])
# y.range <- range(pca.loadings[,3])
# xlab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")	
# ylab <- paste(c("Principal Component 3 (",variancePer[3],"%)"),collapse="")
# plot(pca.loadings[,2],pca.loadings[,3],type='n',xlab=xlab,ylab=ylab,bg="gray")
# points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=col.2.use,cex=2,pch=19)
# points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col="black",cex=2,pch=21)
# legend("topright",c("GFP","R17", "Fingolimod"),cex=.8,col=c("green","red","blue"),pch=c(19,19))
# abline(v=0,lty=2)
# dev.off()
# 

#-------------------------------------------------------------------------------
# Save data for heatmap visualization
#-------------------------------------------------------------------------------
write.table(sig_raw, '../vis_data/f_phospho_limma_selection.txt', sep = '\t')
write.table(lll[coded.keepers,], '../vis_data/f_selected_phospho.txt', sep = '\t')

