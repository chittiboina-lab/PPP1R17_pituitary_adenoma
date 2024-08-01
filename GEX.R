#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
library(stringr)
library(sctree)
library(scSorter)
library(SeuratDisk)
library(SoupX)
library(CellChat)

paths <- c('data/SN_MO/P1_T_outs', 
           'data/SN_MO/P1_N_outs',	
           'data/SN_MO/P2_T_outs',	
           'data/SN_MO/P2_N_outs',	
           'data/SN_MO/P3_T_outs',	
           'data/SN_MO/P4_T_outs')
names <- c()

# Loop over input paths
for (path in paths) {
  
  name <- substr(path, 12, 15)
  names <- c(names, name)

  ## Load data
  obj <- Read10X_h5(stringr::str_interp('${path}/${name}_outs_filtered_feature_bc_matrix.h5'))
  rna_counts <- obj$`Gene Expression`
  
  ## Load Raw Data
  data_raw <- Read10X_h5(stringr::str_interp('${path}/${name}_outs_raw_feature_bc_matrix.h5'))
  data_raw <- data_raw$`Gene Expression`
  
  ##-------------------------------------------------------##
  ## Run SoupX
  ## Adapted from Patrick Fletcher and
  ## https://cellgeni.readthedocs.io/en/latest/notebooks.html
  ##-------------------------------------------------------##
  # SoupChannel(raw, filtered)
  sc <- SoupChannel(data_raw, rna_counts)
  
  ## Preliminary clustering
  prelim <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = stringr::str_interp("${name}")
  )
  
  prelim <- NormalizeData(prelim)

  prelim <- FindVariableFeatures(prelim)
  prelim <- ScaleData(prelim)
  prelim <- RunPCA(prelim)
  prelim <- RunUMAP(prelim, dims = 1:30)
  prelim <- FindNeighbors(prelim, dims = 1:30)
  prelim <- FindClusters(prelim, resolution = 0.5)
  
  meta <- prelim@meta.data
  umap <- prelim@reductions$umap@cell.embeddings

  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc <- setDR(sc, umap)
  
  ## Estimate Contamination
  sc = autoEstCont(sc)
  
  ## Remove Contamination
  adjusted_counts = adjustCounts(sc)
  
  ## -------------------------------------------------------##
  ## Continue Pipeline
  ## -------------------------------------------------------##
  ## create a Seurat object containing the adjusted RNA data
  
  obj <- CreateSeuratObject(
    counts = adjusted_counts,
    assay = "RNA",
    project = stringr::str_interp("${name}")
  )
  
  if (dir.exists('results/')) { } else {
    dir.create('results/')
  }
  
  if (dir.exists(stringr::str_interp('results/${name}/'))) { } else {
    dir.create(stringr::str_interp('results/${name}/'))
  }
  
  # add surgical annotation to metadata
  if (endsWith(name, "T")) {
    obj$anno = "core"
  } 
  else if (endsWith(name, "N")) {
    obj$anno = "margin"
  } 
  
  ##-------------------------------------------------------##
  ## Quality Control & Pre-Processing
  ##-------------------------------------------------------##
  ## Calculate and visualize QC Stats
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('results/${name}/${name}_QC_stats.png'), plot = QC_stats)
  
  ## Remove empty droplets and low quality cells
  obj <- subset(obj, nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 25)
  
  post_cutoff_QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('results/${name}/${name}_post_cutoff_QC_stats.png'), plot = post_cutoff_QC_stats)
  
  ## Remove erythrocytes
  obj <- subset(obj, HBB == 0 & HBA1 == 0 & HBA2 == 0)
  
  ## Log Normalize
  obj <- NormalizeData(obj)
  
  ## Find Variable Features
  obj <- FindVariableFeatures(obj)
  
  ## Scale
  obj <- ScaleData(obj, features = rownames(data))
  
  ## Run PCA
  obj <- RunPCA(obj)
  
  ## Run UMAP
  obj <- RunUMAP(obj, dims = 1:30)
  
  ## Find and Remove Doublets
  # pK Identification (no ground-truth)
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- Idents(obj)
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- obj@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) ## Unsure of most appropriate parameters, these are from the example on the doublet finder github
  
  doublet_col <- rev(names(obj@meta.data))[1] # this column in meta.data have a different name for every sample and/or parameter setting, so I store the name of the column to refer to it later without hard coding
  doublets <- DimPlot(obj, reduction = 'umap', group.by =doublet_col)
  ggsave(stringr::str_interp('results/${name}/${name}_doublets.png'), plot = doublets)
  
  ## actually remove doublets
  # obj <- subset(x = obj, subset = obj@meta.data[[doublet_col]] == 'Singlet') # Doesn't work, the alternative approach below does
  obj <- obj[, obj@meta.data[, doublet_col] == "Singlet"]
  
  ##-------------------------------------------------------##
  ## Classify Individual Cells
  ##-------------------------------------------------------##
  
  ## Using scSorter
  topgenes <- head(VariableFeatures(obj), 2000)
  
  obj_exp = GetAssayData(obj)
  topgene_filter = rowSums(as.matrix(obj_exp)[topgenes, ]!=0) > ncol(obj_exp)*.1
  topgenes = topgenes[topgene_filter]
  
  ## load annotation file
  anno <- read.csv("data/anno.csv")
  
  picked_genes = unique(c(anno$Marker, topgenes))
  obj_exp = obj_exp[rownames(obj_exp) %in% picked_genes, ]
  
  ## run scSorter
  rts <- scSorter(obj_exp, anno)
  
  obj$pred.type <- rts$Pred_Type
 
  ##-------------------------------------------------------##
  ## Clustering
  ##-------------------------------------------------------##
  ## Cluster Cells
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5)
  
  ## Visualize before analysis
  clusters <- DimPlot(obj, reduction = 'umap')
  ggsave(stringr::str_interp('results/${name}/${name}_dimplot.png'), plot = clusters)
  cell_types <- DimPlot(obj, reduction = 'umap', group.by = 'pred.type')
  ggsave(stringr::str_interp('results/${name}/${name}_cell_types.png'), plot = cell_types)
  
  ##-------------------------------------------------------##
  ## Save .h5ad for future visualization in python
  ##-------------------------------------------------------##
  # SaveH5Seurat(obj, filename = stringr::str_interp("data/${name}_GEX.h5Seurat"), overwrite = TRUE)
  i <- sapply(obj@meta.data, is.factor)
  obj@meta.data[i] <- lapply(obj@meta.data[i], as.character)
  as.h5Seurat(obj, stringr::str_interp("results/${name}/${name}_GEX.h5Seurat"))
  Convert(stringr::str_interp("results/${name}/${name}_GEX.h5Seurat"), dest = "h5ad")
  
  ## Export Corticotroph Barcodes
  Idents(obj) <- obj$pred.type
  write.csv(WhichCells(obj, idents = 'corticotrophs'), stringr::str_interp('results/${name}/${name}_cort_barcodes.csv'))
  
  ## Tabulate Cell Types
  write.csv(table(Idents(obj)), stringr::str_interp('results/${name}/${name}_cell_types.csv'))
  
  ## Find Marker genes
  markers <- FindAllMarkers(obj)
  write.csv(markers, stringr::str_interp("results/${name}/${name}_markers.csv"))
  
  #save object to unique name in the environment
  assign(name, obj)
}

names <- c('P1_N', 'P1_T', 'P2_N', 'P2_T', 'P3_T', 'P4_T')

# rm(list = ls(clusters, doublets, bcmvn_obj, post_cutoff_QC_stats, QC_stats, ref, rna_counts, sweep.res.list_obj, sweep.stats_obj,))
dir.create('results/integrated/')

obj.list <- mget(names)

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated, features = rownames(data))
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = 'pca', dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = 'pca', dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

Idents(integrated) <- 'pred.type'

DefaultAssay(integrated) <- "RNA"
cort_DEGs <- FindMarkers(integrated, ident.1 = 'corticotrophs')
write.csv(cort_DEGs, 'results/integrated/panCort_markers.csv')

DEGs <- FindAllMarkers(integrated)
write.csv(DEGs, 'results/integrated/all_markers.csv')

p1 <- DimPlot(integrated)
p2 <- DimPlot(integrated, group.by = 'seurat_clusters')
plot <- p1 + p2 
ggsave('results/integrated/int_celltype_clusters .png', plot, width = 16, height = 16)

# Secondary Analysis (reclustering)
corts <- subset(integrated, subset = pred.type == 'corticotrophs')
DefaultAssay(corts) <- 'RNA'
Idents(corts) <- corts$anno
T_cort_markers <- FindMarkers(corts, ident.1 = 'core')
write.csv(T_cort_markers, 'results/integrated/Core_vs_Margin_cort_markers.csv')

corts <- FindVariableFeatures(corts)
corts <- ScaleData(corts, features = row.names(corts))
corts <- RunPCA(corts, npcs = 30)
corts <- RunUMAP(corts, reduction = 'pca', dims = 1:30)
corts <- FindNeighbors(corts, reduction = 'pca', dims = 1:30)
corts <- FindClusters(corts, resolution = 0.4)

p1 <- DimPlot(corts)
p2 <- DimPlot(corts, group.by = 'anno')
p3 <- DimPlot(corts, group.by = 'orig.ident')
combined <- p1+p2+p3
ggsave('results/integrated/corts_clust_anno_origident.png', combined, width = 24, height = 8)

p1 <- FeaturePlot(corts, features = c('TBX19', 'POMC','PMAIP1', 'PPP1R17'))
ggsave('results/integrated/corts_markers.png', p1, width = 16, height = 16)

cort_clust_markers <- FindAllMarkers(corts)
write.csv(cort_clust_markers, 'results/integrated/int_corts_clust_markers.csv')

## Save integrated corts object in anndata format
i <- sapply(corts@meta.data, is.factor)
corts@meta.data[i] <- lapply(corts@meta.data[i], as.character)
as.h5Seurat(corts, "results/integrated/integrated_corts_GEX.h5Seurat")
Convert("results/integrated/integrated_corts_GEX.h5Seurat", dest = "h5ad")

## Save integrated object in anndata format
DefaultAssay(integrated) <- "RNA"
i <- sapply(integrated@meta.data, is.factor)
integrated@meta.data[i] <- lapply(integrated@meta.data[i], as.character)
as.h5Seurat(integrated, "results/integrated/integrated_GEX.h5Seurat")
Convert("results/integrated/integrated_GEX.h5Seurat", dest = "h5ad")

# Secondary Analysis (reclustering)
corts <- subset(integrated, subset = pred.type == 'corticotrophs')
DefaultAssay(corts) <- 'RNA'
Idents(corts) <- corts$anno
T_cort_markers <- FindMarkers(corts, ident.1 = 'core')
write.csv(T_cort_markers, 'results/integrated/Core_vs_Margin_cort_markers.csv')

corts <- FindVariableFeatures(corts)
corts <- ScaleData(corts, features = row.names(corts))
corts <- RunPCA(corts, npcs = 30)
corts <- RunUMAP(corts, reduction = 'pca', dims = 1:30)
corts <- FindNeighbors(corts, reduction = 'pca', dims = 1:30)
corts <- FindClusters(corts, resolution = 0.4)

p1 <- DimPlot(corts)
p2 <- DimPlot(corts, group.by = 'anno')
p3 <- DimPlot(corts, group.by = 'orig.ident')
combined <- p1+p2+p3
ggsave('results/integrated/corts_clust_anno_origident.png', combined, width = 24, height = 8)

p1 <- FeaturePlot(corts, features = c('TBX19', 'POMC','PMAIP1', 'PPP1R17'))
ggsave('results/integrated/corts_markers.png', p1, width = 16, height = 16)

cort_clust_markers <- FindAllMarkers(corts)
write.csv(cort_clust_markers, 'results/integrated/int_corts_clust_markers.csv')

## Save integrated corts object in anndata format
i <- sapply(corts@meta.data, is.factor)
corts@meta.data[i] <- lapply(corts@meta.data[i], as.character)
as.h5Seurat(corts, "results/integrated/integrated_corts_GEX.h5Seurat")
Convert("results/integrated/integrated_corts_GEX.h5Seurat", dest = "h5ad")

## Save integrated object in anndata format
DefaultAssay(integrated) <- "RNA"
i <- sapply(integrated@meta.data, is.factor)
integrated@meta.data[i] <- lapply(integrated@meta.data[i], as.character)
as.h5Seurat(integrated, "results/integrated/integrated_GEX.h5Seurat")
Convert("results/integrated/integrated_GEX.h5Seurat", dest = "h5ad")

## Integrated Cell Chat analysis, with corts separated by tumor/normal

# integrated <- LoadH5Seurat('results/integrated/integrated_GEX.h5Seurat')

integrated$cell_type <- ifelse(integrated$pred.type == 'corticotrophs' & integrated$anno == 'core', "Tumor_Corts", integrated$pred.type)
integrated$cell_type <- ifelse(integrated$pred.type == 'corticotrophs' & integrated$anno == 'margin', "Normal_Corts", integrated$cell_type)

dir.create('results/cellchat/')
dir.create('results/cellchat/celltypes/')
dir.create('results/cellchat/pathways/')

cellchat <- createCellChat(object = integrated, group.by = "cell_type")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

PathwayNet <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(PathwayNet, stringr::str_interp('results/cellchat/pathways.csv'))

AggCountsNet <- subsetCommunication(cellchat, slot.name = "counts")
write.csv(AggCountsNet, stringr::str_interp('results/cellchat/counts.csv'))

AggWeightsNet <- subsetCommunication(cellchat, slot.name = "weights")
write.csv(AggWeightsNet, stringr::str_interp('results/cellchat/weights.csv'))

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
png(stringr::str_interp('results/cellchat/interactions_count.png'))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = c('#9f4800','#8b2be2','#1ac938','#ff7c00','#a3a3a3','#023eff','#f14cc1','#00d7ff','#e8000b','#ffc400'))
dev.off()
png(stringr::str_interp('results/cellchat/interactions_weight.png'))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = c('#9f4800','#8b2be2','#1ac938','#ff7c00','#a3a3a3','#023eff','#f14cc1','#00d7ff','#e8000b','#ffc400'))
dev.off()


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  png(stringr::str_interp('results/cellchat/celltypes/${rownames(mat)[i]}_interactions.png'))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = c('#9f4800','#8b2be2','#1ac938','#ff7c00','#a3a3a3','#023eff','#f14cc1','#00d7ff','#e8000b','#ffc400'))
  dev.off()
}

pathways <- cellchat@netP[1]$pathways
for (pathway in pathways) {
  
  pathways.show <- c(pathway) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # Circle plot
  par(mfrow=c(1,1))
  png(stringr::str_interp('results/cellchat/pathways/${pathway}.png'))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use = c('#9f4800','#8b2be2','#1ac938','#ff7c00','#a3a3a3','#023eff','#f14cc1','#00d7ff','#e8000b','#ffc400'))
  dev.off()
}

png('results/cellchat/pathways_from_T.png', width = 9, height = 9, units = 'in', res = 600)
netVisual_chord_gene(cellchat, sources.use = 9, lab.cex = 0.5,legend.pos.y = 30, legend.pos.x = 0, #scale = TRUE
                     color.use = c('#9f4800','#8b2be2','#1ac938','#ff7c00','#a3a3a3','#023eff','#f14cc1','#00d7ff','#e8000b','#ffc400'))
dev.off()
png('results/cellchat/pathways_from_FSC.png', width = 9, height = 9, units = 'in', res = 600)
netVisual_chord_gene(cellchat, sources.use = 2, lab.cex = 0.5,legend.pos.y = 30, legend.pos.x = 0, #scale = TRUE
                     color.use = c('#9f4800','#8b2be2','#1ac938','#ff7c00','#a3a3a3','#023eff','#f14cc1','#00d7ff','#e8000b','#ffc400'))
dev.off()

saveRDS(cellchat, file = stringr::str_interp('results/cellchat/cellchat.RDS'))

## Analysis of previously published single cell datasets
paths <- c('data/SC_GEX/p1_cd_M.h5',	'data/SC_GEX/p2_cd_T.h5',	'data/SC_GEX/p3_cd_T.h5',	'data/SC_GEX/p4_gh_T.h5',
           'data/SC_GEX/p2_cd_N.h5',	'data/SC_GEX/p3_cd_N.h5',	'data/SC_GEX/p4_gh_N.h5',	'data/SC_GEX/p5_nf_M.h5')
names <- c()

# Loop over input paths
for (path in paths) {
  
  name <- substr(path, 13, 19)
  names <- c(names, name)
  
  ## Load data
  obj <- Read10X_h5(path)
  
  ## -------------------------------------------------------##
  ## Continue Pipeline
  ## -------------------------------------------------------##
  ## create a Seurat object containing the adjusted RNA data
  
  obj <- CreateSeuratObject(obj,
                            assay = "RNA",
                            project = stringr::str_interp('${name}')
  )
  
  if (dir.exists('results/SC_GEX/')) { } else {
    dir.create('results/SC_GEX/')
  }
  
  if (dir.exists(stringr::str_interp('results/SC_GEX/${name}/'))) { } else {
    dir.create(stringr::str_interp('results/SC_GEX/${name}/'))
  }
  
  # add surgical annotation to metadata
  if (endsWith(name, "T")) {
    obj$anno = "core"
  } 
  else if (endsWith(name, "N")) {
    obj$anno = "margin"
  } 
  
  ##-------------------------------------------------------##
  ## Quality Control & Pre-Processing
  ##-------------------------------------------------------##
  ## Calculate and visualize QC Stats
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('results/SC_GEX/${name}/${name}_QC_stats.png'), plot = QC_stats)
  
  ## Remove empty droplets and low quality cells
  obj <- subset(obj, nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 25)
  
  post_cutoff_QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('results/SC_GEX/${name}/${name}_post_cutoff_QC_stats.png'), plot = post_cutoff_QC_stats)
  
  ## Remove erythrocytes
  obj <- subset(obj, HBB == 0 & HBA1 == 0 & HBA2 == 0)
  
  ## Log Normalize
  obj <- NormalizeData(obj)
  
  ## Find Variable Features
  obj <- FindVariableFeatures(obj)
  
  ## Scale
  obj <- ScaleData(obj, features = rownames(data))
  
  ## Run PCA
  obj <- RunPCA(obj)
  
  ## Run UMAP
  obj <- RunUMAP(obj, dims = 1:30)
  
  ## Find and Remove Doublets
  # pK Identification (no ground-truth)
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- Idents(obj)
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- obj@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) ## Unsure of most appropriate parameters, these are from the example on the doublet finder github
  
  doublet_col <- rev(names(obj@meta.data))[1] # this column in meta.data have a different name for every sample and/or parameter setting, so I store the name of the column to refer to it later without hard coding
  doublets <- DimPlot(obj, reduction = 'umap', group.by =doublet_col)
  ggsave(stringr::str_interp('results/SC_GEX/${name}/${name}_doublets.png'), plot = doublets)
  
  ## actually remove doublets
  # obj <- subset(x = obj, subset = obj@meta.data[[doublet_col]] == 'Singlet') # Doesn't work, the alternative approach below does
  obj <- obj[, obj@meta.data[, doublet_col] == "Singlet"]
  
  ##-------------------------------------------------------##
  ## Classify Individual Cells
  ##-------------------------------------------------------##
  
  ## Using scSorter
  topgenes <- head(VariableFeatures(obj), 2000)
  
  obj_exp = GetAssayData(obj)
  topgene_filter = rowSums(as.matrix(obj_exp)[topgenes, ]!=0) > ncol(obj_exp)*.1
  topgenes = topgenes[topgene_filter]
  
  ## load annotation file
  anno <- read.csv("data/anno.csv")
  
  picked_genes = unique(c(anno$Marker, topgenes))
  obj_exp = obj_exp[rownames(obj_exp) %in% picked_genes, ]
  
  ## run scSorter
  rts <- scSorter(obj_exp, anno)
  
  obj$pred.type <- rts$Pred_Type
  
  ##-------------------------------------------------------##
  ## Clustering
  ##-------------------------------------------------------##
  ## Cluster Cells
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5)
  
  ## Visualize before analysis
  clusters <- DimPlot(obj, reduction = 'umap')
  ggsave(stringr::str_interp('results/SC_GEX/${name}/${name}_dimplot.png'), plot = clusters)
  cell_types <- DimPlot(obj, reduction = 'umap', group.by = 'pred.type')
  ggsave(stringr::str_interp('results/SC_GEX/${name}/${name}_cell_types.png'), plot = cell_types)
  
  ##-------------------------------------------------------##
  ## Save .h5ad for future visualization in python
  ##-------------------------------------------------------##
  # SaveH5Seurat(obj, filename = stringr::str_interp("data/${name}_GEX.h5Seurat"), overwrite = TRUE)
  i <- sapply(obj@meta.data, is.factor)
  obj@meta.data[i] <- lapply(obj@meta.data[i], as.character)
  as.h5Seurat(obj, stringr::str_interp("results/SC_GEX/${name}/${name}_GEX.h5Seurat"))
  Convert(stringr::str_interp("results/SC_GEX/${name}/${name}_GEX.h5Seurat"), dest = "h5ad")
  
  ## Tabulate Cell Types
  write.csv(table(Idents(obj)), stringr::str_interp('results/SC_GEX/${name}/${name}_cell_types.csv'))
  
  ## Find Marker genes
  markers <- FindAllMarkers(obj)
  write.csv(markers, stringr::str_interp("results/SC_GEX/${name}/${name}_markers.csv"))
  
  #save object to unique name in the environment
  assign(name, obj)
}
# rm(list = ls(clusters, doublets, bcmvn_obj, post_cutoff_QC_stats, QC_stats, ref, rna_counts, sweep.res.list_obj, sweep.stats_obj,))
dir.create('results/SC_GEX/integrated/')

obj.list <- mget(names)

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated, features = rownames(data))
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = 'pca', dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = 'pca', dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

Idents(integrated) <- 'pred.type'

DEGs <- FindAllMarkers(integrated)
write.csv(DEGs, 'results/SC_GEX/integrated/all_markers.csv')

DefaultAssay(integrated) <- "RNA"
i <- sapply(integrated@meta.data, is.factor)
integrated@meta.data[i] <- lapply(integrated@meta.data[i], as.character)
as.h5Seurat(integrated, "results/SC_GEX/integrated/integrated_SC.h5Seurat")
Convert("results/SC_GEX/integrated/integrated_SC.h5Seurat", dest = "h5ad")

