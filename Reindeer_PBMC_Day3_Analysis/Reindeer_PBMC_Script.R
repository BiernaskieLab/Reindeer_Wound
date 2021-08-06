setwd("C:/Users/jaffe/Desktop/Reindeer_Blood_Seurat")

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/jaffe/Desktop/Reindeer_Blood_Seurat/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "Reindeer_PBMC", min.cells = 3, min.features = 200)
pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 4 & nCount_RNA < 26000)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
write.csv(pbmc@meta.data, file = "metadata_reindeer.csv")
meta_data = read.csv(file = "metadata_reindeer.csv")
pbmc$sample_ident = meta_data$sample_ident
View(pbmc@meta.data)
save(pbmc, file = "pbmc.Robj")

Idents(object = pbmc) <- "sample_ident"
just_pbmc = subset(pbmc, idents = c("5XD_Back_PBMC", "5XD_Antler_PBMC", "4XE_Back_PBMC", "4XE_Antler_PBMC"))
save(just_pbmc, file = "just_pbmc.Robj")

load(file = "pbmc.Robj")
load(file = "just_pbmc.Robj")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
just_pbmc <- NormalizeData(just_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

just_pbmc <- FindVariableFeatures(just_pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(just_pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(just_pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

###Do both objects separately now cause memory
load(file = "just_pbmc.Robj")

all.genes <- rownames(just_pbmc)
just_pbmc <- ScaleData(just_pbmc, features = all.genes)

just_pbmc <- RunPCA(just_pbmc, features = VariableFeatures(object = just_pbmc))
# Examine and visualize PCA results a few different ways
print(just_pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(just_pbmc, dims = 1:20, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(just_pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(just_pbmc, ndims = 30)
just_pbmc = pbmc

#B_plasma_cells <- JackStraw(B_plasma_cells, num.replicate = 100, dims = 1:50)
#B_plasma_cells <- ScoreJackStraw(B_plasma_cells, dims = 1:50)
#JackStrawPlot(B_plasma_cells, dims = 1:50)
save(just_pbmc, file = "just_pbmc.Robj")
rm(list = ls())

load(file = "just_pbmc.Robj")

just_pbmc <- FindNeighbors(just_pbmc, dims = 1:20)
just_pbmc <- FindClusters(just_pbmc, resolution = 0.6)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
just_pbmc <- RunUMAP(just_pbmc, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(just_pbmc, reduction = "umap")

######### Cell Types
#Blood cells. Blood contains many types of cells: white blood cells (monocytes, lymphocytes, neutrophils, eosinophils, basophils, and macrophages), red blood cells (erythrocytes), and platelets.

#T cells
FeaturePlot(just_pbmc, features = c("CD3D", "CD3E", "CD8A", "CD4", "TRAC", "TRBC1", "TRDC"), min.cutoff = "q10")

#B cells
FeaturePlot(just_pbmc, features = c("CD19", "MS4A1", "CD79A"), min.cutoff = "q10")

#Monocytes
FeaturePlot(just_pbmc, features = c("CD14", "ITGAM", "CCR2", "CD68"), min.cutoff = "q10")

#Neutrophil
#Neutrophils	CD15, CD16, CD49d(-)
FeaturePlot(just_pbmc, features = c("S100A8", "S100A9", "ITGAM", "FCGR3A", "ITGB2", "FCGR2A", "FCGR2B", "CD44", "CD55"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("FCGR2A", "ITGAM", "CD44", "S100A8", "S100A9", "CD55"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("FUT4", "FCGR3A", "ITGA4"), min.cutoff = "q10")

#Plasma (SDC1)
FeaturePlot(just_pbmc, features = c("IGHG1", "IGHG3", "IGHG4", "IGLC3", "JCHAIN", "IGHGP", "IGHG2", "IGLL5", "MZB1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(just_pbmc, features = c("MROH7", "ST6GALNAC4", "SPAG4", "SIK1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(just_pbmc, features = c("SDC1"), min.cutoff = "q10")

#Platelets (resting)	CD42b; Platelets (activated)	CD62P
FeaturePlot(just_pbmc, features = c("ITGA2B", "SELP", "ITGB3"), min.cutoff = "q10")

#Red Blood Cells/Eythrocyte; Red Blood Cells	CD235a
FeaturePlot(just_pbmc, features = c("CD34", "PTPRC", "GYPA", "TFRC", "CD47"), min.cutoff = "q10")

#Dendritic; CD1c, CD83, CD141, CD209, MHC II
FeaturePlot(just_pbmc, features = c("PTGDS", "GZMB", "IRF7", "CST3", "LILRA4", "ITGAX", "HLA-DRA", "CD1B", "ZBTB46", "FLT3", "CD209"))
FeaturePlot(just_pbmc, features = c("FSCN1", "CCL19", "CCL17", "LAMP3", "LGALS2", "TMEM176A", "EBI3", "GSN", "TMEM176B", "TFPI2"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("CD1C", "CD83", "CD141", "CD209", "HLA-DR"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("XCR1", "CD1C", "CCR7", "LAMP3", "LAMP", "CLEC9A"), min.cutoff = "q10")

#NK Cells
FeaturePlot(just_pbmc, features = c("KLRB1", "KLRD1", "NCR1", "IL2RB", "IL12RB2"), min.cutoff = "q10")

#Stem cells
FeaturePlot(just_pbmc, features = c("CD34", "THY1"), min.cutoff = "q10")

just_pbmc_all_markers = FindAllMarkers(just_pbmc)
save(just_pbmc_all_markers, file = "just_pbmc_all_markers.csv")


####pbmc + tissue
load(file = "pbmc.Robj")
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(just_pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 20:30, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(just_pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc, ndims = 50)
just_pbmc = pbmc

pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.6)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:30)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

######### Cell Types
#Blood cells. Blood contains many types of cells: white blood cells (monocytes, lymphocytes, neutrophils, eosinophils, basophils, and macrophages), red blood cells (erythrocytes), and platelets.

#T cells
FeaturePlot(pbmc, features = c("CD3D", "CD3E", "CD8A", "CD4", "TRAC", "TRBC1", "TRDC"), min.cutoff = "q10")

#B cells
FeaturePlot(pbmc, features = c("CD19", "MS4A1", "CD79A"), min.cutoff = "q10")

#Monocytes
FeaturePlot(pbmc, features = c("CD14", "ITGAM", "CCR2", "CD68"), min.cutoff = "q10")

#Neutrophil
#Neutrophils	CD15, CD16, CD49d(-)
FeaturePlot(pbmc, features = c("S100A8", "S100A9", "ITGAM", "FCGR3A", "ITGB2", "FCGR2A", "FCGR2B", "CD44", "CD55"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("FCGR2A", "ITGAM", "CD44", "S100A8", "S100A9", "CD55"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("FUT4", "FCGR3A", "ITGA4"), min.cutoff = "q10")

#Plasma (SDC1)
FeaturePlot(pbmc, features = c("IGHG1", "IGHG3", "IGHG4", "IGLC3", "JCHAIN", "IGHGP", "IGHG2", "IGLL5", "MZB1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(pbmc, features = c("MROH7", "ST6GALNAC4", "SPAG4", "SIK1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(pbmc, features = c("SDC1"), min.cutoff = "q10")

#Platelets (resting)	CD42b; Platelets (activated)	CD62P
FeaturePlot(pbmc, features = c("ITGA2B", "SELP", "ITGB3"), min.cutoff = "q10")

#Red Blood Cells/Eythrocyte; Red Blood Cells	CD235a
FeaturePlot(pbmc, features = c("CD34", "PTPRC", "GYPA", "TFRC", "CD47"), min.cutoff = "q10")

#Dendritic; CD1c, CD83, CD141, CD209, MHC II
FeaturePlot(pbmc, features = c("PTGDS", "GZMB", "IRF7", "CST3", "LILRA4", "ITGAX", "HLA-DRA", "CD1B", "ZBTB46", "FLT3", "CD209"))
FeaturePlot(pbmc, features = c("FSCN1", "CCL19", "CCL17", "LAMP3", "LGALS2", "TMEM176A", "EBI3", "GSN", "TMEM176B", "TFPI2"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("CD1C", "CD83", "CD141", "CD209", "HLA-DR"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("XCR1", "CD1C", "CCR7", "LAMP3", "LAMP", "CLEC9A"), min.cutoff = "q10")

#NK Cells
FeaturePlot(pbmc, features = c("KLRB1", "KLRD1", "NCR1", "IL2RB", "IL12RB2"), min.cutoff = "q10")

#Stem cells
FeaturePlot(pbmc, features = c("CD34", "THY1"), min.cutoff = "q10")

pbmc_all_markers = FindAllMarkers(pbmc)
save(pbmc_all_markers, file = "pbmc_all_markers.csv")

jpeg(file = "reindeer_ectopic_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("KRT14", "KRT17", "KRT15", "MT4", "KRT5", "KRT14"), min.cutoff = 'q2')
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("TH", "KIT"))
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("CD68", "FOLR2", "FOLR3", "NOS2", "CSF1R", "DUSP2"), min.cutoff = 'q2')
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_CD45_IL1B_Myeloid.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("IL1B", "SRGN", "PTPRC"))
dev.off()

jpeg(file = "reindeer_ectopic_selected_markers_neutrophil.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("S100A8", "S100A9"))
dev.off()


###Integration Script
setwd("C:/Users/jaffe/Desktop/Reindeer_Blood_Seurat/Integration")

library(Seurat)
library(patchwork)

load(file = "just_pbmc.Robj")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(just_pbmc, split.by = "sample_ident")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
ElbowPlot(immune.combined, ndims = 30)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample_ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

just_pbmc = immune.combined
DefaultAssay(just_pbmc) <- "integrated"

######### Cell Types
#Blood cells. Blood contains many types of cells: white blood cells (monocytes, lymphocytes, neutrophils, eosinophils, basophils, and macrophages), red blood cells (erythrocytes), and platelets.

#T cells
FeaturePlot(just_pbmc, features = c("CD3D", "CD3E", "CD8A", "CD4", "TRAC", "TRBC1", "TRDC"), min.cutoff = "q10")

#B cells
FeaturePlot(just_pbmc, features = c("CD19", "MS4A1", "CD79A"), min.cutoff = "q10")

#Monocytes
FeaturePlot(just_pbmc, features = c("CD14", "ITGAM", "CCR2", "CD68"), min.cutoff = "q10")

#Neutrophil
#Neutrophils	CD15, CD16, CD49d(-)
FeaturePlot(just_pbmc, features = c("S100A8", "S100A9", "ITGAM", "FCGR3A", "ITGB2", "FCGR2A", "FCGR2B", "CD44", "CD55"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("FCGR2A", "ITGAM", "CD44", "S100A8", "S100A9", "CD55"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("FUT4", "FCGR3A", "ITGA4"), min.cutoff = "q10")

#Plasma (SDC1)
FeaturePlot(just_pbmc, features = c("IGHG1", "IGHG3", "IGHG4", "IGLC3", "JCHAIN", "IGHGP", "IGHG2", "IGLL5", "MZB1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(just_pbmc, features = c("MROH7", "ST6GALNAC4", "SPAG4", "SIK1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(just_pbmc, features = c("SDC1"), min.cutoff = "q10")

#Platelets (resting)	CD42b; Platelets (activated)	CD62P
FeaturePlot(just_pbmc, features = c("ITGA2B", "SELP", "ITGB3"), min.cutoff = "q10")

#Red Blood Cells/Eythrocyte; Red Blood Cells	CD235a
FeaturePlot(just_pbmc, features = c("CD34", "PTPRC", "GYPA", "TFRC", "CD47"), min.cutoff = "q10")

#Dendritic; CD1c, CD83, CD141, CD209, MHC II
FeaturePlot(just_pbmc, features = c("PTGDS", "GZMB", "IRF7", "CST3", "LILRA4", "ITGAX", "HLA-DRA", "CD1B", "ZBTB46", "FLT3", "CD209"))
FeaturePlot(just_pbmc, features = c("FSCN1", "CCL19", "CCL17", "LAMP3", "LGALS2", "TMEM176A", "EBI3", "GSN", "TMEM176B", "TFPI2"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("CD1C", "CD83", "CD141", "CD209", "HLA-DR"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("XCR1", "CD1C", "CCR7", "LAMP3", "LAMP", "CLEC9A"), min.cutoff = "q10")

#NK Cells
FeaturePlot(just_pbmc, features = c("KLRB1", "KLRD1", "NCR1", "IL2RB", "IL12RB2"), min.cutoff = "q10")

#Stem cells
FeaturePlot(just_pbmc, features = c("CD34", "THY1"), min.cutoff = "q10")

just_pbmc_all_markers = FindAllMarkers(just_pbmc)
save(just_pbmc_all_markers, file = "just_pbmc_all_markers.csv")
save(just_pbmc, file = "just_pbmc.Robj")


####
rm(list = ls())
load(file = "pbmc.Robj")
ifnb.list <- SplitObject(pbmc, split.by = "sample_ident")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 40, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
ElbowPlot(immune.combined, ndims = 40)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample_ident")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

pbmc = immune.combined
DefaultAssay(pbmc) <- "RNA"

######### Cell Types
#Blood cells. Blood contains many types of cells: white blood cells (monocytes, lymphocytes, neutrophils, eosinophils, basophils, and macrophages), red blood cells (erythrocytes), and platelets.

#T cells
FeaturePlot(pbmc, features = c("CD3D", "CD3E", "CD8A", "CD4", "TRAC", "TRBC1", "TRDC"), min.cutoff = "q10")

#B cells
FeaturePlot(pbmc, features = c("CD19", "MS4A1", "CD79A"), min.cutoff = "q10")

#Monocytes
FeaturePlot(pbmc, features = c("CD14", "ITGAM", "CCR2", "CD68"), min.cutoff = "q10")

#Neutrophil
#Neutrophils	CD15, CD16, CD49d(-)
FeaturePlot(pbmc, features = c("S100A8", "S100A9", "ITGAM", "FCGR3A", "ITGB2", "FCGR2A", "FCGR2B", "CD44", "CD55"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("FCGR2A", "ITGAM", "CD44", "S100A8", "S100A9", "CD55"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("FUT4", "FCGR3A", "ITGA4"), min.cutoff = "q10")

#Plasma (SDC1)
FeaturePlot(pbmc, features = c("IGHG1", "IGHG3", "IGHG4", "IGLC3", "JCHAIN", "IGHGP", "IGHG2", "IGLL5", "MZB1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(pbmc, features = c("MROH7", "ST6GALNAC4", "SPAG4", "SIK1"), min.cutoff = "q10", raster=FALSE)
FeaturePlot(pbmc, features = c("SDC1"), min.cutoff = "q10")

#Platelets (resting)	CD42b; Platelets (activated)	CD62P
FeaturePlot(pbmc, features = c("ITGA2B", "SELP", "ITGB3"), min.cutoff = "q10")

#Red Blood Cells/Eythrocyte; Red Blood Cells	CD235a
FeaturePlot(pbmc, features = c("CD34", "PTPRC", "GYPA", "TFRC", "CD47"), min.cutoff = "q10")

#Dendritic; CD1c, CD83, CD141, CD209, MHC II
FeaturePlot(pbmc, features = c("PTGDS", "GZMB", "IRF7", "CST3", "LILRA4", "ITGAX", "HLA-DRA", "CD1B", "ZBTB46", "FLT3", "CD209"))
FeaturePlot(pbmc, features = c("FSCN1", "CCL19", "CCL17", "LAMP3", "LGALS2", "TMEM176A", "EBI3", "GSN", "TMEM176B", "TFPI2"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("CD1C", "CD83", "CD141", "CD209", "HLA-DR"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("XCR1", "CD1C", "CCR7", "LAMP3", "LAMP", "CLEC9A"), min.cutoff = "q10")

#NK Cells
FeaturePlot(pbmc, features = c("KLRB1", "KLRD1", "NCR1", "IL2RB", "IL12RB2"), min.cutoff = "q10")

#Stem cells
FeaturePlot(pbmc, features = c("CD34", "THY1"), min.cutoff = "q10")

pbmc_all_markers = FindAllMarkers(pbmc)
save(pbmc_all_markers, file = "pbmc_all_markers.csv")
save(pbmc, file = "pbmc.Robj")


DefaultAssay(just_pbmc) <- "RNA"
FeaturePlot(just_pbmc, features = c("S100A8", "S100A9", "CSF3R", "CSF1R"), min.cutoff = "q10")

#CSF3R - neutrophils

DimPlot(just_pbmc, label = TRUE)
table(just_pbmc$seurat_clusters)
FeaturePlot(just_pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

FeaturePlot(just_pbmc, features = c("S100A8", "S100A9"))

# Cust 0 = S100A8/A9 Myeloid
# S100A9
# S100A8

# Cust 1 = CD4+T
# CD4+ T Cells

# Cust 2 = CD8+T
# CD8+ T Cells

# Cust 3 = CD247+T

# Cust 4 = MS4A1+CD79B+ B Cells

# Cust 5 = UK

# Cust 6 = Megakaryocyte
# GNG11+ Megakaryocyte

# Cust 7 = NKG2A.1+ NK Cells

just_pbmc_all_markers = FindAllMarkers(just_pbmc)
write.csv(just_pbmc_all_markers, file = "just_pbmc_all_markers.csv")
FeaturePlot(just_pbmc, features = c("S100A8", "S100A9"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("CD4", "CD8A", "CD247", "CD3E", "CD3D"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("MS4A1", "CD79B", "CD19"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("GNG11", "ITGA2B", "SELP", "ITGB3"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("NKG2A.1"), min.cutoff = "q10")
FeaturePlot(just_pbmc, features = c("NKG2A.1", "KLRB1", "KLRD1", "NCR1", "IL2RB", "IL12RB2"), min.cutoff = "q10")


load(file = "pbmc.Robj")
DimPlot(pbmc, label = TRUE)
DimPlot(pbmc, label = TRUE, split.by = "sample_ident")
FeaturePlot(pbmc, features = c("S100A8", "S100A9", "CSF3R", "CSF1R", "CD68", "CD14", "FCGR3A"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("CD4", "CD8A", "CD247", "CD3E", "CD3D"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("MS4A1", "CD79B", "CD19"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("GNG11", "ITGA2B", "SELP", "ITGB3"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("NKG2A.1", "KLRB1", "KLRD1", "NCR1", "IL2RB", "IL12RB2"), min.cutoff = "q10")
FeaturePlot(pbmc, features = c("IRF8", "BLA-DQB", "BOLA-DQB"), min.cutoff = "q10")

cluster_14 = FindMarkers(pbmc, ident.1 = "14")
write.csv(cluster_14, file = "cluster_14.csv")
cluster_13 = FindMarkers(pbmc, ident.1 = "13")
write.csv(cluster_13, file = "cluster_13.csv")



load(file = "just_pbmc.Robj")
DimPlot(just_pbmc)
just_pbmc = subset(just_pbmc, idents = c("19", "20"), invert = TRUE)
new.cluster.ids <- c("CD8_T", "CD4_T", "Myeloid", "Myeloid", "CD8_T", "CD4_T", "Myeloid", "CD247_T", "CD8_T", "CD4_T", "B_cell", "Other",
                     "Megakaryocyte", "NK_cell", "Other", "CD8_T", "CD247_T", "DC", "Myeloid")
names(new.cluster.ids) <- levels(just_pbmc)
just_pbmc <- RenameIdents(just_pbmc, new.cluster.ids)
DimPlot(just_pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
save(just_pbmc, file = "just_pbmc.Robj")
just_pbmc@meta.data$annotations = just_pbmc@active.ident
table(just_pbmc@active.ident)

Idents(just_pbmc) = "sample_ident"
Antler = subset(just_pbmc, idents = c("5XD_Antler_PBMC", "4XE_Antler_PBMC"))
DimPlot(Antler)
table(Antler@meta.data$annotations)
Back = subset(just_pbmc, idents = c("5XD_Antler_PBMC", "4XE_Antler_PBMC"), invert = TRUE)
DimPlot(Back)
table(Back@meta.data$annotations)


load(file = "pbmc.Robj")
DefaultAssay(pbmc) <- "RNA"
FeaturePlot(pbmc, features = c("S100A8", "S100A9", "CSF3R", "CSF1R", "CD68", "CD14", "FCGR3A"))



cluster_7220 = FindMarkers(pbmc, ident.1 = c("7", "2", "20"))
write.csv(cluster_7220, file = "cluster_7220.csv")
FeaturePlot(pbmc, features = c("KRT14", "KRT5"))
FeaturePlot(pbmc, features = c("CSF3R", "CSF1R", "S100A8", "S100A9"))

cluster_15 = FindMarkers(pbmc, ident.1 = "15")
write.csv(cluster_15, file = "cluster_15.csv")
cluster_11 = FindMarkers(pbmc, ident.1 = "11")
write.csv(cluster_11, file = "cluster_11.csv")

#https://www.dropbox.com/sh/gdtnor459o91bhx/AABmSYvPIXEXrejc3tZ0Qufya?dl=0

write.csv(pbmc@meta.data, file = "metadata_reindeer.csv")
meta_data = read.csv(file = "metadata_reindeer.csv")
pbmc$sample_orig = meta_data$sample_orig
save(pbmc, file = "pbmc.Robj")
DimPlot(pbmc)
new.cluster.ids <- c("Myeloid", "T_cells", "Keratinocytes", "T_cells", "Myeloid", "Tissue", "Myeloid", "Keratinocytes", "Tissue", "T_cells", "Myeloid", "Myeloid",
                     "T_cells", "Unknown", "Unknown", "Myeloid", "Myeloid", "Tissue", "Tissue", "B_cells", "Keratinocytes", "DC")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
save(pbmc, file = "pbmc.Robj")
pbmc@meta.data$annotations = pbmc@active.ident

install.packages("ape")
library(ape)
library(phytools)
plot.phylo(pbmc)
pbmc@meta.data$tip.label = pbmc@meta.data$sample_orig
pbmc@meta.data$edge = pbmc@meta.data$sample_orig
pbmc@meta.data$Nnode = pbmc@meta.data$sample_orig

install.packages('Rcpp')
library(Rcpp)

setwd("C:/Users/jaffe/Desktop/Reindeer_Blood_Seurat/Integration")
load(file = "pbmc.Robj")
library(Seurat)
DefaultAssay(pbmc) <- "integrated"

Idents(pbmc) = "seurat_clusters"
Idents(pbmc) = "sample_orig"
Idents(pbmc) = "annotations"
pbmc <- BuildClusterTree(object = pbmc)
PlotClusterTree(object = pbmc)

Myeloid = subset(pbmc, idents = c("Myeloid"))
DefaultAssay(Myeloid) <- "integrated"
Idents(Myeloid) = "sample_orig"
Idents(Myeloid) = "seurat_clusters"
Myeloid <- BuildClusterTree(object = Myeloid)
PlotClusterTree(object = Myeloid)


T_cells = subset(pbmc, idents = c("T_cells"))
DefaultAssay(T_cells) <- "integrated"
Idents(T_cells) = "sample_orig"
Idents(T_cells) = "seurat_clusters"
T_cells <- BuildClusterTree(object = T_cells)
PlotClusterTree(object = T_cells)
