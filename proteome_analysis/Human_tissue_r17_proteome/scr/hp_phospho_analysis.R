#------------------------------------------------------------------------------
# Source Functions
#------------------------------------------------------------------------------

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

lll <- read.table("data/hp_data_phospho.txt",header=T,row.names=1)
lll <- lll[,c(4,1,5,2,6,3)]
col.2.use <- c(rep(c(3,2),3))

#------------------------------------------------------------------------------
# Impute Missing Data
#------------------------------------------------------------------------------

new.lll <- NULL
for(i in 1:dim(lll)[1]) {
  temp1 <- as.numeric(lll[i,])
  temp2 <- sum(is.na(temp1))
  if(temp2>0) {
    print(i)
    temp3 <- mean(temp1,na.rm=T)
    temp4 <- ifelse(is.na(temp1),temp3,temp1)
    new.lll <- rbind(new.lll,temp4)
  } else {
    new.lll <- rbind(new.lll,temp1)
  }
}
dimnames(new.lll) <- dimnames(lll)

#------------------------------------------------------------------------------
# Perform pedestal and transformation
#------------------------------------------------------------------------------

temp1 <- new.lll
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
legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(15,15))
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
legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(15,15))
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
lines(c(pca.loadings[1,1],pca.loadings[2,1]),c(pca.loadings[1,2],pca.loadings[2,2]),lty=1,col=1)
lines(c(pca.loadings[3,1],pca.loadings[4,1]),c(pca.loadings[3,2],pca.loadings[4,2]),lty=1,col=1)
lines(c(pca.loadings[5,1],pca.loadings[6,1]),c(pca.loadings[5,2],pca.loadings[6,2]),lty=1,col=1)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(19,19))
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

#-------------------------------------------------------------------------------
# Generate Noise Model Plot
#-------------------------------------------------------------------------------

xrange <- range(c(1:14))
yrange <- range(c(0:50))
plot(xrange,yrange,type='n',xlab="Mean(CyclicLoess(Log2(Detection Value+2)))",ylab="Coefficient of Variation")
points(g1.mean,g1.cv,pch=21,col=8)
points(g2.mean,g2.cv,pch=21,col=8)
lines(lowess.g1,col=3,lwd=2)
lines(lowess.g2,col=2,lwd=2)
abline(h=0)
abline(v=4.5,lty=2,lwd=2)
legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(19,19))
pdf("phospho/Noise_model_plot.pdf")
plot(xrange,yrange,type='n',xlab="Mean(CyclicLoess(Log2(Detection Value+2)))",ylab="Coefficient of Variation")
points(g1.mean,g1.cv,pch=21,col=8)
points(g2.mean,g2.cv,pch=21,col=8)
lines(lowess.g1,col=3,lwd=2)
lines(lowess.g2,col=2,lwd=2)
abline(h=0)
abline(v=4.5,lty=2,lwd=2)
legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(19,19))
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
lines(c(pca.loadings[1,1],pca.loadings[2,1]),c(pca.loadings[1,2],pca.loadings[2,2]),lty=1,col=1)
lines(c(pca.loadings[3,1],pca.loadings[4,1]),c(pca.loadings[3,2],pca.loadings[4,2]),lty=1,col=1)
lines(c(pca.loadings[5,1],pca.loadings[6,1]),c(pca.loadings[5,2],pca.loadings[6,2]),lty=1,col=1)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("topright",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(19,19))
dev.off()


#-------------------------------------------------------------------------------
# Significance Testing
#-------------------------------------------------------------------------------

running.pvs <- NULL
for (i in 1:dim(Final.working.data)[1]) {
  print(i)
  temp1 <- as.numeric(Final.working.data[i,c(1,3,5)])
  temp2 <- as.numeric(Final.working.data[i,c(2,4,6)])
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
# Filter select
#-------------------------------------------------------------------------------

temp1 <- running.pvs[,2]>1.25
temp2 <- running.pvs[,3]<0.10
temp3 <- (temp1+temp2)==2
passing.mcc.frags <- names(temp3[temp3==TRUE])
length(passing.mcc.frags)
running.pvs <- cbind(running.pvs,as.logical(temp3))
dimnames(running.pvs)[[2]] <- c("CD_vs_Ctrl_Difference","Abs_CD_vs_Ctrl_Difference","Paired_Ttest_uncorrected_P","Paired_Ttest_Corrected_P","Pass_Selection_0_No_1_Yes")
write.table(running.pvs, "phospho/differential_testing_results.txt", sep = "\t")

coded.keepers <- dimnames(running.pvs)[[1]][running.pvs[,5]==1]
grab.keeper.info <- running.pvs[coded.keepers,]
write.table(grab.keeper.info, "phospho/selected_differential.txt", sep = "\t")

#-------------------------------------------------------------------------------
# Generate Volcano Plot and PCA depicting Gene Selection for Class Differences
#-------------------------------------------------------------------------------

nochange <- setdiff(dimnames(running.pvs)[[1]],dimnames(grab.keeper.info)[[1]])
up <- dimnames(grab.keeper.info)[[1]][grab.keeper.info[,1]>0]
down <- dimnames(grab.keeper.info)[[1]][grab.keeper.info[,1]<0]
xrange <- range(running.pvs[,1])
yrange <- range(c(0,max(-log(running.pvs[,3])+10)))

pdf("phospho/post_selection_volcano_plot.pdf")
plot(xrange,yrange,xlab="Linear Fold Change",ylab="-log(Paired T-test Uncorrected P)",type='n',main="CD vs Ctrl")
points(running.pvs[,1][nochange],-log(running.pvs[,3][nochange]),col=8,cex=1.5,pch=19)
points(running.pvs[,1][up],-log(running.pvs[,3][up]),col=2,cex=1.5,pch=19)
points(running.pvs[,1][down],-log(running.pvs[,3][down]),col=3,cex=1.5,pch=19)
abline(v=-1)
abline(v=1)
abline(h=-log(0.10),lty=2)
abline(v=-1.25,lty=2)
abline(v=1.25,lty=2)
temp <- legend("topright", legend = c(" ", " ", " "),text.width=strwidth("100000000000"),pch=c(19,19,19),xjust=1,yjust=1,title="Difference",col=c(2,8,3),cex=1.5)
a <- paste(c("Up (n=",length(up),")"),collapse="")
b <- paste(c("None (n=",length(nochange),")"),collapse="")
c <- paste(c("Down (n=",length(down),")"),collapse="")
text(temp$rect$left + temp$rect$w, temp$text$y,c(a,b,c), pos = 2)
dev.off()

#-------------------------------------------------------------------------------
# Post Selection Visuals - PCA
#-------------------------------------------------------------------------------

pca.results <- princomp(Final.working.data[coded.keepers,],cor=F)
pca.loadings <- loadings(pca.results)
phosphovar <- (pca.results$sdev^2)
variancePer <- round(phosphovar/sum(phosphovar)*100,1)
variancePer <- variancePer[1:3]
palette("default")

pdf("phospho/post_selection_PCA.pdf")
par(mfrow=c(1,2))
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
lines(c(pca.loadings[1,1],pca.loadings[2,1]),c(pca.loadings[1,2],pca.loadings[2,2]),lty=1,col=1)
lines(c(pca.loadings[3,1],pca.loadings[4,1]),c(pca.loadings[3,2],pca.loadings[4,2]),lty=1,col=1)
lines(c(pca.loadings[5,1],pca.loadings[6,1]),c(pca.loadings[5,2],pca.loadings[6,2]),lty=1,col=1)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("top",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(19,19))
abline(h=0,lty=2)
x.range <- range(pca.loadings[,2])
y.range <- range(pca.loadings[,3])
xlab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")	
ylab <- paste(c("Principal Component 3 (",variancePer[3],"%)"),collapse="")
plot(pca.loadings[,2],pca.loadings[,3],type='n',xlab=xlab,ylab=ylab,bg="gray")
lines(c(pca.loadings[1,2],pca.loadings[2,2]),c(pca.loadings[1,3],pca.loadings[2,3]),lty=1,col=1)
lines(c(pca.loadings[3,2],pca.loadings[4,2]),c(pca.loadings[3,3],pca.loadings[4,3]),lty=1,col=1)
lines(c(pca.loadings[5,2],pca.loadings[6,2]),c(pca.loadings[5,3],pca.loadings[6,3]),lty=1,col=1)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,2],pca.loadings[,3],lwd=0.5,col="black",cex=2,pch=21)
legend("right",c("Ctrl","CD"),cex=.8,col=c("green","red"),pch=c(19,19))
abline(v=0,lty=2)
dev.off()

#------------------------------------------------------------------------------
# Post Selection Visuals - HeatMap (uncluseted)
#------------------------------------------------------------------------------
# 
# cdat <- Final.working.data[coded.keepers,c(1,3,5,2,4,6)]
# cdat <- cor(cdat)
# write.table(cdat, "hp_phospho_cor.txt", sep = '\t')
# 
# my_colour = list(sample = c('Ctrl_1' = "#00FF00",
#                             'Ctrl_2' = "#00FF00",
#                             'Ctrl_3' = "#00FF00",
#                             'CD_1' = "#FF0000",
#                             'CD_2' = "#FF0000",
#                             'CD_3' = "#FF0000"))
# 
# my_sample_col <- data.frame(sample = dimnames(cdat)[[1]])
# row.names(my_sample_col) <- dimnames(cdat)[[1]]
# my_sample_col$sample <- factor(my_sample_col$sample, levels = dimnames(cdat)[[1]])
# 
# pheatmap(cdat,annotation_colors=my_colour,annotation_col=my_sample_col,cellwidth=40,cellheight=40)
# 
# svg("phospho/post_selection_cor_heatmap.svg")
# pheatmap(cdat,annotation_colors=my_colour,annotation_col=my_sample_col,cellwidth=40,cellheight=40)
# dev.off()
# 
# temp1 <- Final.working.data[coded.keepers,]
# write.table(temp1, "hp_phospho_HM.txt", sep = '\t')
# cc <- colorpanel(64, "green", "red")
# heatmap.2(as.matrix(temp1), trace="none", density.info="none",col=cc,labRow="",margins=c(15,5),scale="none")
# pdf("phospho/post_selection_heatmap.pdf")
# heatmap.2(as.matrix(temp1), trace="none", density.info="none",col=cc,labRow="",margins=c(15,5),scale="none")
# dev.off()

#------------------------------------------------------------------------------
# Save Session (if needed)
#------------------------------------------------------------------------------

# save.image("hp_phospho.RData")

