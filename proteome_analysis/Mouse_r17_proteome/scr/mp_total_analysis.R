#------------------------------------------------------------------------------
# Source Functions
#------------------------------------------------------------------------------

library(stats)
library(gplots)
library(ggplot2)
# library(rgl)
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

lll <- read.table("data/mp_data_total.txt",header=T,row.names=1)
col.2.use <- c(rep(3,5),rep(2,5))

#------------------------------------------------------------------------------
# Impute Missing Data
#------------------------------------------------------------------------------

# new.lll <- NULL
# for(i in 1:dim(lll)[1]) {
# 	temp1 <- as.numeric(lll[i,])
# 	temp2 <- sum(is.na(temp1))
# 	if(temp2>0) {
# 		print(i)
# 		temp3 <- mean(temp1,na.rm=T)
# 		temp4 <- ifelse(is.na(temp1),temp3,temp1)
# 		new.lll <- rbind(new.lll,temp4)
# 	} else {
# 		new.lll <- rbind(new.lll,temp1)
# 	}
# }
# dimnames(new.lll) <- dimnames(lll)

#------------------------------------------------------------------------------
# Perform pedestal and transformation
#------------------------------------------------------------------------------

# temp1 <- new.lll
# temp1 <- log(temp1+2,2)
temp1 <- log(lll+2,2)
sample.data <- temp1

#------------------------------------------------------------------------------
# Pre-Normalization Box Plot
#------------------------------------------------------------------------------

pdf('total/pre_normalization_boxplot.pdf')
WorkingList <- vector("list", dim(sample.data)[2])
for ( i in 1:dim(sample.data)[2]) {
	WorkingList[[i]] <- as.numeric(sample.data[,i])
}
boxplot(WorkingList,col=col.2.use,outline=FALSE,xlab="Sample Index",ylab="Log2(Detection Value+2)")
legend("topright",c("GFP","R17"),cex=.8,col=c("green","red"),pch=c(15,15))
dev.off()

#------------------------------------------------------------------------------
# Pre normalization PCA plot
#------------------------------------------------------------------------------

pdf('total/pre_norm_pca.pdf')
pca.results <- princomp(running.normalized,cor=F)
pca.loadings <- loadings(pca.results)
totalvar <- (pca.results$sdev^2)
variancePer <- round(totalvar/sum(totalvar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("right",c("GFP","R17"),cex=.8,col=c("green","red"),pch=c(19,19))
dev.off()

#------------------------------------------------------------------------------
# Perform cross sample normalization
#------------------------------------------------------------------------------

running.normalized <- normalizeBetweenArrays(sample.data,"cyclicloess")

#------------------------------------------------------------------------------
# Post normaliation box plot
#------------------------------------------------------------------------------

pdf('total/post_normalization_boxplot.pdf')
WorkingList <- vector("list", dim(running.normalized)[2])
for ( i in 1:dim(running.normalized)[2]) {
	WorkingList[[i]] <- as.numeric(running.normalized[,i])
}
boxplot(WorkingList,col=col.2.use,outline=FALSE,xlab="Sample Index",ylab="CyclicLoess(Log2(Detection Value+2))")
legend("topright",c("GFP","R17"),cex=.8,col=c("green","red"),pch=c(15,15))
dev.off()
running.normalized <- na.omit(running.normalized)

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

pdf('total/noise_model_plot.pdf')
xrange <- range(c(1:12))
yrange <- range(c(0:25))
plot(xrange,yrange,type='n',xlab="Mean(CyclicLoess(Log2(Detection Value+2)))",ylab="Coefficient of Variation")
points(g1.mean,g1.cv,pch=21,col=8)
points(g2.mean,g2.cv,pch=21,col=8)
lines(lowess.g1,col=3,lwd=2)
lines(lowess.g2,col=2,lwd=2)
abline(h=0)
abline(v=4.5,lty=2,lwd=2)
legend("topright",c("GFP","R17"),cex=.8,col=c("green","red"),lty=1,lwd=2)
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
write.table(Final.working.data,"total/Final_working_data_txt",sep="\t")

#------------------------------------------------------------------------------
# Post normalization PCA plot
#------------------------------------------------------------------------------

pdf('total/post_norm_pca.pdf')
pca.results <- princomp(Final.working.data,cor=F)
pca.loadings <- loadings(pca.results)
totalvar <- (pca.results$sdev^2)
variancePer <- round(totalvar/sum(totalvar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("right",c("GFP","R17"),cex=.8,col=c("green","red"),pch=c(19,19))
dev.off()

#-------------------------------------------------------------------------------
# Significance Testing
#-------------------------------------------------------------------------------

running.pvs <- NULL
for (i in 1:dim(Final.working.data)[1]) {
	print(i)
	temp1 <- as.numeric(Final.working.data[i,c(1:5)])
	temp2 <- as.numeric(Final.working.data[i,c(6:10)])
	temp3 <- t.test(temp1,temp2, paired = TRUE)
	temp4 <- temp3$p.value
	temp5 <- mean(temp2)-mean(temp1)
	temp6 <- 2^temp5
	temp7 <- ifelse(temp6<1,-1/temp6,temp6)
	temp7 <- temp7*(-1) # reverse difference to be R17 vs GFP
	running.pvs <- rbind(running.pvs,c(temp7,abs(temp7),temp4))
}
dimnames(running.pvs)[[1]] <- dimnames(Final.working.data)[[1]]

running.means <- NULL
for (i in 1:dim(Final.working.data)[1]) {
	print(i)
	temp1 <- mean(as.numeric(Final.working.data[i,c(1:5)]))
	temp2 <- mean(as.numeric(Final.working.data[i,c(6:10)]))
	running.means <- rbind(running.means,c(temp1,temp2))
}
dimnames(running.means)[[1]] <- dimnames(Final.working.data)[[1]]

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

temp1 <- running.pvs[,2]>1.5
temp2 <- running.pvs[,4]<0.05
temp3 <- (temp1+temp2)==2
passing.mcc.frags <- names(temp3[temp3==TRUE])
length(passing.mcc.frags)
running.pvs <- cbind(running.pvs,as.logical(temp3))
dimnames(running.pvs)[[2]] <- c("R17_vs_GFP_Difference","R17_vs_GFP_Difference","Paired_Ttest_uncorrected_P","Paired_Ttest_Corrected_P","Pass_Selection_0_No_1_Yes")
write.table(running.pvs, "total/R17_vs_GFP_differential_testing_results.txt", sep = "\t")

coded.keepers <- dimnames(running.pvs)[[1]][running.pvs[,5]==1]
grab.keeper.info <- running.pvs[coded.keepers,]
write.table(grab.keeper.info, "total/R17_vs_GFP_selected_differential.txt", sep = "\t")

#-------------------------------------------------------------------------------
# Post Selection Visuals - PCA
#-------------------------------------------------------------------------------

pca.results <- princomp(Final.working.data[coded.keepers,],cor=F)
pca.loadings <- loadings(pca.results)
totalvar <- (pca.results$sdev^2)
variancePer <- round(totalvar/sum(totalvar)*100,1)
variancePer <- variancePer[1:3]
palette("default")
x.range <- range(pca.loadings[,1])
y.range <- range(pca.loadings[,2])
xlab <- paste(c("Principal Component 1 (",variancePer[1],"%)"),collapse="")	
ylab <- paste(c("Principal Component 2 (",variancePer[2],"%)"),collapse="")
pdf("total/post_selection_PCA.pdf")
plot(pca.loadings[,1],pca.loadings[,2],type='n',xlab=xlab,ylab=ylab,bg="gray")
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col=col.2.use,cex=2,pch=19)
points(pca.loadings[,1],pca.loadings[,2],lwd=0.5,col="black",cex=2,pch=21)
legend("right",c("GFP","R17"),cex=.8,col=c("green","red"),pch=c(19,19))
dev.off()

#-------------------------------------------------------------------------------
# Save data for heatmap visualization
#-------------------------------------------------------------------------------

write.table(lll[coded.keepers,], '../vis_data/mp_selected_total.txt', sep = '\t')
