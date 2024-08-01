library(BiocManager)
library(hdf5r)
library(missMethyl)
library(minfi)
library(limma)
library(DMRcate)
library(ggplot2)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


# setwd("~/OneDrive - National Institutes of Health/Pit_r17_final/methylation/updated_analysis")
# list.files("data1")
# list.files("data2")

basedir1 <- "data1"
targets1 <- read.metharray.sheet(basedir1)

rgSet1 <- read.metharray.exp(targets = targets1, force = TRUE)

basedir2 <- "data2"
targets2 <- read.metharray.sheet(basedir2)

rgSet2 <- read.metharray.exp(targets = targets2, force = TRUE)

z <- minfi:::.harmonizeDataFrames(minfi:::.pDataFix(colData(rgSet1)), minfi:::.pDataFix(colData(rgSet2)))
targets <- rbind(targets1, targets2)

# combine arrays
rgSet <- combineArrays(object1 = rgSet1, object2 = rgSet2,
                       outType= "IlluminaHumanMethylationEPIC")

# get annotations
annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

## Perform subset quantile normalization. SWAN corrects for the technical differences between the Infinium I and II assay designs and produces a smoother overall β
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=TRUE)

png(filename = "results/densitybyprobetype.png", width = 8, height = 4, units = "in", res = 330)
par(mfrow=c(1,2))
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")
dev.off()

## Filter low quality probes
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]

## Extract beta and M values
set.seed(10)
mset_reduced <- mSetSw[sample(1:nrow(mSetSw), 20000),]
meth <- getMeth(mset_reduced)
unmeth <- getUnmeth(mset_reduced)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mset_reduced)

png("results/MDS_plot.png", height = 4, width = 5, units = 'in', res = 330)
par(mfrow=c(1,1))
plotMDS(Mval, labels=targets$sample_name, col=as.integer(factor(targets$disease_type)))
dev.off()

## Create Design Matrix
group <- factor(targets$disease_type,levels=c("normal", "CD"))
# id <- factor(targets$sample_name)
# design <- model.matrix(~id + group)
design <- model.matrix(~group)

fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)

summary(decideTests(fit.reduced))

beta_sig <- beta[which(fit.reduced$p.value[,'groupCD'] < 0.05),]
write.csv(beta_sig, file = 'results/beta_sig_05.csv')

 

#----------------------------------------------------#
## Remove unwanted variation with RUV
# get M-values for ALL probes
meth <- getMeth(mSet)
unmeth <- getUnmeth(mSet)
M <- log2((meth + 100)/(unmeth + 100))
# setup the factor of interest
grp <- factor(targets$disease_type, labels=c(0,1))
# extract Illumina negative control data
INCs <- getINCs(rgSet)

# add negative control data to M-values
Mc <- rbind(M,INCs)
# create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)
# table(ctl1)
rfit1 <- RUVfit(Y = Mc, X = grp, ctl = ctl1) # Stage 1 analysis
rfit2 <- RUVadj(Y = Mc, fit = rfit1)

top1 <- topRUV(rfit2, num=Inf, p.BH = 1)

ctl2 <- rownames(M) %in% rownames(top1[top1$p.BH_X1.1 > 0.5,])
table(ctl2)

# Perform RUV adjustment and fit
rfit3 <- RUVfit(Y = M, X = grp, ctl = ctl2) # Stage 2 analysis
rfit4 <- RUVadj(Y = M, fit = rfit3)
# Look at table of top results
topRUV(rfit4)

# Perform RUV adjustment and fit with specific k value
rfit5 <- RUVfit(Y = M, X = grp, ctl = ctl2, method = "ruv4", k=1) # Stage 2 analysis
rfit6 <- RUVadj(Y = M, fit = rfit5)
# Perform RUV adjustment and fit with specific k value
rfit7 <- RUVfit(Y = M, X = grp, ctl = ctl2, method = "ruv4", k=2) # Stage 2 analysis
rfit8 <- RUVadj(Y = M, fit = rfit5)

## Visualize effects of RUV adjustment
Madj1 <- getAdj(M, rfit3) # get adjusted values
Madj2 <- getAdj(M, rfit5)
Madj3 <- getAdj(M, rfit7)

png(filename = "results/MDS_unadj_vs_ruv.png", width = 8, height = 8, units = "in", res = 330)
par(mfrow=c(2,2))
plotMDS(M, labels=targets$sample_name, col=as.integer(factor(targets$disease_type)),
        main="Unadjusted", gene.selection = "common"
        )

plotMDS(Madj1, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Adjusted: RUV-inverse", gene.selection = "common"
        )
plotMDS(Madj2, labels=targets$sample_name, col=as.integer(factor(targets$disease_type)),
        main="RUV4 K=1", gene.selection = "common"
)

plotMDS(Madj3, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Adjusted: RUV4 K=2", gene.selection = "common"
)

dev.off()

## Maybe k=1 / Madj2 / rfit5 is best?
#----------------------------------------------------#


## Region level analysis and GO
myAnnotation <- cpg.annotate(object = M, datatype = "array", what = "M", 
                             arraytype = c("EPIC"), 
                             analysis.type = "differential", design = design, 
                             coef = 2)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
write.table(x = data.frame(results.ranges), file = "results/DMRs.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )

png('results/probe_number_bias.png', width = 6, height = 6, units = "in", res = 330)
gst.region <- goregion(results.ranges, all.cpg=rownames(M), 
                       collection="GO", array.type="EPIC", plot.bias=TRUE)
dev.off()

cols <- c(2,4)[group]
names(cols) <-group
beta <- getBeta(mSet)

png('results/top_DMR.png', width = 10, height = 12, units = "in", res = 330)
par(mfrow=c(1,1))
DMR.plot(ranges=results.ranges, dmr=93, CpGs=M, phen.col=cols, 
         what="M", arraytype="EPIC", genome="hg19")
dev.off()

topGSA(gst.region, n=100)
