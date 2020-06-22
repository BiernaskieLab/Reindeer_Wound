#### Reindeer Single-Cell Analysis May 3

tmux
bsub -Is bash

export PATH="/home/ajaffer/anaconda3/bin:$PATH"
conda env list
source activate Seurat_v3.1
cd /home/ajaffer/Reindeer_Single_Cell_Analysis_Bovine/
  R
getwd()
sessionInfo()

library(dplyr)
library(Seurat)
library(patchwork)

### Day 0 Sample Redo

# Load the dataset
R3D0.data <- Read10X(data.dir = "/home/ajaffer/Reindeer_Round3_FASTQs/Round3_Agg/Round3_Agg_Day0_Bovine/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data):
R3D0 <- CreateSeuratObject(counts = R3D0.data, project = "reindeer_RNA", min.cells = 2, min.features = 100)
R3D0


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
R3D0[["percent.mt"]] <- PercentageFeatureSet(R3D0, pattern = "^mt-")
# Visualize QC metrics as a violin plot
jpeg(file = "VlnPlot_R3D0_Redo.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
### The bovine reference genome either is not annotated with mitochondrial DNA or the reindeer sample did not align to them
### The min cutoff was not included because we wanted to ganulocytes which have low transcriptional expression
R3D0 <- subset(R3D0, subset = nFeature_RNA < 4000)
R3D0 <- subset(R3D0, subset = nCount_RNA < 30000)
jpeg(file = "VlnPlot_R3D0_Redo_After_Subset.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


write.csv(R3D0@meta.data, file = "R3D0@meta.data.csv")
### Modify CSV on Excel (-1 = back, -2 = antler)
meta_data = read.csv(file = "R3D0@meta.data.csv")
new_col = data.frame(meta_data$sample_ident)
R3D0$sample_ident = meta_data$sample_ident

R3D0.list <- SplitObject(R3D0, split.by = "sample_ident")
R3D0.list <- R3D0.list[c("Back", "Antler")]

for (i in 1:length(R3D0.list)) {
  R3D0.list[[i]] <- NormalizeData(R3D0.list[[i]], verbose = FALSE)
  R3D0.list[[i]] <- FindVariableFeatures(R3D0.list[[i]], selection.method = "vst",
                                         nfeatures = 2000, verbose = FALSE)
}

save(R3D0, file = "R3D0.Robj")

R3D0.list <- R3D0.list[c("Back", "Antler")]

########### Example with D0 - dims 1:10,  ###########
R3D0.anchors_PC10 <- FindIntegrationAnchors(object.list = R3D0.list, dims = 1:10)
R3D0.integrated_PC10 <- IntegrateData(anchorset = R3D0.anchors_PC10, dims = 1:10)
DefaultAssay(R3D0.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D0.integrated_PC10 <- ScaleData(R3D0.integrated_PC10, verbose = FALSE)
R3D0.integrated_PC10 <- RunPCA(R3D0.integrated_PC10, npcs = 10, verbose = FALSE)
R3D0.integrated_PC10 <- RunUMAP(R3D0.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "R3D0.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D0 - dims 1:20,  ###########
R3D0.anchors_PC20 <- FindIntegrationAnchors(object.list = R3D0.list, dims = 1:20)
R3D0.integrated_PC20 <- IntegrateData(anchorset = R3D0.anchors_PC20, dims = 1:20)
DefaultAssay(R3D0.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D0.integrated_PC20 <- ScaleData(R3D0.integrated_PC20, verbose = FALSE)
R3D0.integrated_PC20 <- RunPCA(R3D0.integrated_PC20, npcs = 20, verbose = FALSE)
R3D0.integrated_PC20 <- RunUMAP(R3D0.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "R3D0.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D0 - dims 1:30,  ###########
R3D0.anchors_PC30 <- FindIntegrationAnchors(object.list = R3D0.list, dims = 1:30)
R3D0.integrated_PC30 <- IntegrateData(anchorset = R3D0.anchors_PC30, dims = 1:30)
DefaultAssay(R3D0.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D0.integrated_PC30 <- ScaleData(R3D0.integrated_PC30, verbose = FALSE)
R3D0.integrated_PC30 <- RunPCA(R3D0.integrated_PC30, npcs = 30, verbose = FALSE)
R3D0.integrated_PC30 <- RunUMAP(R3D0.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "R3D0.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D0 - dims 1:40,  ###########
R3D0.anchors_PC40 <- FindIntegrationAnchors(object.list = R3D0.list, dims = 1:40)
R3D0.integrated_PC40 <- IntegrateData(anchorset = R3D0.anchors_PC40, dims = 1:40)
DefaultAssay(R3D0.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D0.integrated_PC40 <- ScaleData(R3D0.integrated_PC40, verbose = FALSE)
R3D0.integrated_PC40 <- RunPCA(R3D0.integrated_PC40, npcs = 40, verbose = FALSE)
R3D0.integrated_PC40 <- RunUMAP(R3D0.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "R3D0.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:50,  ###########
R3D0.anchors_PC50 <- FindIntegrationAnchors(object.list = R3D0.list, dims = 1:50)
R3D0.integrated_PC50 <- IntegrateData(anchorset = R3D0.anchors_PC50, dims = 1:50)
DefaultAssay(R3D0.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D0.integrated_PC50 <- ScaleData(R3D0.integrated_PC50, verbose = FALSE)
R3D0.integrated_PC50 <- RunPCA(R3D0.integrated_PC50, npcs = 50, verbose = FALSE)
R3D0.integrated_PC50 <- RunUMAP(R3D0.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "R3D0.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 0 PC 50
jpeg(file = "R3D0.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(R3D0.integrated_PC50, ndims = 50)
dev.off()
# We are using PC 10 for this sample

#### Day 0 Annotations

R3D0.integrated_PC10 <- FindNeighbors(R3D0.integrated_PC10, dims = 1:10)
R3D0.integrated_PC10 <- FindClusters(R3D0.integrated_PC10, resolution = 0.6)

write.csv(R3D0.integrated_PC10@meta.data, file = "R3D0.integrated_PC10@meta.data.csv")
Idents(R3D0.integrated_PC10) <- "seurat_clusters"

jpeg(file = "R3D0.integrated_PC10_0.6_res_clustering.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC10, reduction = "umap", label = T)
dev.off()

save(R3D0.integrated_PC10, file = "R3D0.integrated_PC10.Robj")

jpeg(file = "R3D0_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                               "CRABP1"))
dev.off()

jpeg(file = "R3D0_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("KRT14", "KRT17", "KRT15"))
dev.off()

jpeg(file = "R3D0_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "R3D0_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("TH", "KIT"))
dev.off()

jpeg(file = "R3D0_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "R3D0_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("CD68", "FOLR2", "FOLR3"))
dev.off()

jpeg(file = "R3D0_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "R3D0_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "R3D0_selected_markers_Unknown.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("AIF1", "NR4A3", "PTPRC", "IL17A", "CD74", "S100A8", "S100A9"))
dev.off()

jpeg(file = "R3D0_selected_markers_CSFR1.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("CSFR1"))
dev.off()

jpeg(file = "R3D0_selected_markers_Endothelial_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("PECAM1", "ACKR1", "PTN", "APOA1", "MT4"))
dev.off()

jpeg(file = "R3D0_selected_markers_Keratinocytes_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D0.integrated_PC10, features = c("KRT14", "KRT17", "KRT15", "MT4"))
dev.off()


Idents(R3D0.integrated_PC10) <- "seurat_clusters"
new.cluster.ids <- c("Fibroblast", "Basal Keratinocyte", "Basal Keratinocyte", "Endothelial", "Suprabasal Keratinocyte", "CD45+IL1B+Myeloid", "Fibroblast", "Fibroblast", "Melanocyte", "Basal Keratinocyte", "T-cell", "Macrophage", "VSM", "Macrophage", "Fibroblast", "Schwann", "Macrophage")
names(new.cluster.ids) <- levels(R3D0.integrated_PC10)
R3D0.integrated_PC10 <- RenameIdents(R3D0.integrated_PC10, new.cluster.ids)

jpeg(file = "R3D0.integrated_PC10_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC10, reduction = "umap", label = T)
dev.off()

jpeg(file = "R3D0.integrated_PC10_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC10, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

jpeg(file = "R3D0.integrated_PC10_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D0.integrated_PC10, reduction = "umap", label = F)
dev.off()

save(R3D0.integrated_PC10, file = "R3D0.integrated_PC10.Robj")


### Day 3 Sample Redo

# Load the dataset
R3D3.data <- Read10X(data.dir = "/home/ajaffer/Reindeer_Round3_FASTQs/Round3_Agg/Round3_Agg_Day3_Bovine/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data):
R3D3 <- CreateSeuratObject(counts = R3D3.data, project = "reindeer_RNA", min.cells = 2, min.features = 100)
R3D3


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
R3D3[["percent.mt"]] <- PercentageFeatureSet(R3D3, pattern = "^mt-")
# Visualize QC metrics as a violin plot
jpeg(file = "VlnPlot_R3D3_Redo.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
### The bovine reference genome either is not annotated with mitochondrial DNA or the reindeer sample did not align to them
### The min cutoff was not included because we wanted to ganulocytes which have low transcriptional expression
R3D3 <- subset(R3D3, subset = nFeature_RNA < 5900)
R3D3 <- subset(R3D3, subset = nCount_RNA < 50000)
jpeg(file = "VlnPlot_R3D3_Redo_After_Subset.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


write.csv(R3D3@meta.data, file = "R3D3@meta.data.csv")
### Modify CSV on Excel (-1 = back, -2 = antler)
meta_data = read.csv(file = "R3D3@meta.data.csv")
new_col = data.frame(meta_data$sample_ident)
R3D3$sample_ident = meta_data$sample_ident

R3D3.list <- SplitObject(R3D3, split.by = "sample_ident")
R3D3.list <- R3D3.list[c("Back", "Antler")]

for (i in 1:length(R3D3.list)) {
  R3D3.list[[i]] <- NormalizeData(R3D3.list[[i]], verbose = FALSE)
  R3D3.list[[i]] <- FindVariableFeatures(R3D3.list[[i]], selection.method = "vst",
                                         nfeatures = 2000, verbose = FALSE)
}

save(R3D3, file = "R3D3.Robj")

load(file = "R3D3.Robj")

R3D3.list <- R3D3.list[c("Back", "Antler")]

########### Example with D0 - dims 1:10,  ###########
R3D3.anchors_PC10 <- FindIntegrationAnchors(object.list = R3D3.list, dims = 1:10)
R3D3.integrated_PC10 <- IntegrateData(anchorset = R3D3.anchors_PC10, dims = 1:10)
DefaultAssay(R3D3.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D3.integrated_PC10 <- ScaleData(R3D3.integrated_PC10, verbose = FALSE)
R3D3.integrated_PC10 <- RunPCA(R3D3.integrated_PC10, npcs = 10, verbose = FALSE)
R3D3.integrated_PC10 <- RunUMAP(R3D3.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "R3D3.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D0 - dims 1:20,  ###########
R3D3.anchors_PC20 <- FindIntegrationAnchors(object.list = R3D3.list, dims = 1:20)
R3D3.integrated_PC20 <- IntegrateData(anchorset = R3D3.anchors_PC20, dims = 1:20)
DefaultAssay(R3D3.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D3.integrated_PC20 <- ScaleData(R3D3.integrated_PC20, verbose = FALSE)
R3D3.integrated_PC20 <- RunPCA(R3D3.integrated_PC20, npcs = 20, verbose = FALSE)
R3D3.integrated_PC20 <- RunUMAP(R3D3.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "R3D3.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D0 - dims 1:30,  ###########
R3D3.anchors_PC30 <- FindIntegrationAnchors(object.list = R3D3.list, dims = 1:30)
R3D3.integrated_PC30 <- IntegrateData(anchorset = R3D3.anchors_PC30, dims = 1:30)
DefaultAssay(R3D3.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D3.integrated_PC30 <- ScaleData(R3D3.integrated_PC30, verbose = FALSE)
R3D3.integrated_PC30 <- RunPCA(R3D3.integrated_PC30, npcs = 30, verbose = FALSE)
R3D3.integrated_PC30 <- RunUMAP(R3D3.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "R3D3.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D0 - dims 1:40,  ###########
R3D3.anchors_PC40 <- FindIntegrationAnchors(object.list = R3D3.list, dims = 1:40)
R3D3.integrated_PC40 <- IntegrateData(anchorset = R3D3.anchors_PC40, dims = 1:40)
DefaultAssay(R3D3.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D3.integrated_PC40 <- ScaleData(R3D3.integrated_PC40, verbose = FALSE)
R3D3.integrated_PC40 <- RunPCA(R3D3.integrated_PC40, npcs = 40, verbose = FALSE)
R3D3.integrated_PC40 <- RunUMAP(R3D3.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "R3D3.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:50,  ###########
R3D3.anchors_PC50 <- FindIntegrationAnchors(object.list = R3D3.list, dims = 1:50)
R3D3.integrated_PC50 <- IntegrateData(anchorset = R3D3.anchors_PC50, dims = 1:50)
DefaultAssay(R3D3.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D3.integrated_PC50 <- ScaleData(R3D3.integrated_PC50, verbose = FALSE)
R3D3.integrated_PC50 <- RunPCA(R3D3.integrated_PC50, npcs = 50, verbose = FALSE)
R3D3.integrated_PC50 <- RunUMAP(R3D3.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "R3D3.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 3 PC 50
jpeg(file = "R3D3.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(R3D3.integrated_PC50, ndims = 50)
dev.off()
# We are using PC 20 for this sample?

########### Day 3 Cell Type Annotations ###############

R3D3.integrated_PC20 <- FindNeighbors(R3D3.integrated_PC20, dims = 1:20)
R3D3.integrated_PC20 <- FindClusters(R3D3.integrated_PC20, resolution = 0.6)

write.csv(R3D3.integrated_PC20@meta.data, file = "R3D3.integrated_PC20@meta.data.csv")
Idents(R3D3.integrated_PC20) <- "seurat_clusters"

jpeg(file = "R3D3.integrated_PC20_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20, reduction = "umap", label = T)
dev.off()
###### ARZINA REMEMBER TO SAVE THIS, save image (line 669) ######
save(R3D3.integrated_PC20, file = "R3D3.integrated_PC20.Robj")

jpeg(file = "R3D3_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT", "CRABP1"))
dev.off()

jpeg(file = "R3D3_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("KRT14", "KRT17", "KRT15"))
dev.off()

jpeg(file = "R3D3_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "R3D3_selected_markers_Immune_Cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("AIF1", "NR4A3"))
dev.off()

jpeg(file = "R3D3_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("TH", "KIT"))
dev.off()

jpeg(file = "R3D3_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "R3D3_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("CD68", "FOLR2", "FOLR3"))
dev.off()

jpeg(file = "R3D3_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "R3D3_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "R3D3_selected_markers_Real_Endo.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("PECAM1", "FLT1", "LYVE1"))
dev.off()

jpeg(file = "R3D3_selected_markers_PTPRC.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("PTPRC"))
dev.off()

jpeg(file = "R3D3_selected_markers_neutrophil.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("S100A8", "S100A9"))
dev.off()

jpeg(file = "R3D3_selected_markers_IL17_CD74.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("IL17A", "CD74"))
dev.off()

jpeg(file = "R3D3_selected_markers_Macro_Myeloid.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("CD68", "NOS2", "DUSP2", "NCF1"))
dev.off()

jpeg(file = "R3D3_selected_markers_CSFR1_NOS2.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("CSF1R", "NOS2"))
dev.off()

jpeg(file = "R3D3_selected_markers_Endothelial_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("PECAM1", "ACKR1", "PTN", "APOA1", "MT4"))
dev.off()

jpeg(file = "R3D3_selected_markers_Keratinocytes_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D3.integrated_PC20, features = c("KRT14", "KRT17", "KRT15", "MT4"))
dev.off()


Idents(R3D3.integrated_PC20) <- "seurat_clusters"
new.cluster.ids <- c("CD68+DUSP2+Macrophage", "CD68+NOS+Macrophage", "Arterial Endothelial", "Basal Keratinocyte", "CD68+CSF1R+Macrophage", "Basal Keratinocyte", "Fibroblast", "Suprabasal Keratinocyte", "CD68+NOS+Macrophage", "Fibroblast", "Melanocyte", "Suprabasal Keratinocyte", "IL17+T-cell", "CD74+T-cell", "Lymphatic Endothelial", "CD68+NCF1+Myeloid", "VSM", "UK", "Arterial Endothelial")
names(new.cluster.ids) <- levels(R3D3.integrated_PC20)
R3D3.integrated_PC20 <- RenameIdents(R3D3.integrated_PC20, new.cluster.ids)

jpeg(file = "R3D3.integrated_PC20_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20, reduction = "umap", label = T)
dev.off()

jpeg(file = "R3D3.integrated_PC20_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

R3D3.integrated_PC20_subset = subset(R3D3.integrated_PC20, idents = c("CD68+DUSP2+Macrophage", "CD68+NOS+Macrophage", "Arterial Endothelial", "Basal Keratinocyte", "CD68+CSF1R+Macrophage", "Basal Keratinocyte", "Fibroblast", "Suprabasal Keratinocyte", "CD68+NOS+Macrophage", "Fibroblast", "Melanocyte", "Suprabasal Keratinocyte", "IL17+T-cell", "CD74+T-cell", "Lymphatic Endothelial", "CD68+NCF1+Myeloid", "VSM", "Arterial Endothelial"))

R3D3.integrated_PC20_subset = subset(R3D3.integrated_PC20, idents != c("UK"))


jpeg(file = "R3D3.integrated_PC20_annotated_labels_subset.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20_subset, reduction = "umap", label = T)
dev.off()

jpeg(file = "R3D3.integrated_PC20_annotated_split_Back_Antler_subset.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20_subset, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

jpeg(file = "R3D3.integrated_PC20_annotated_labels_subset.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D3.integrated_PC20_subset, reduction = "umap", label = F)
dev.off()

save(R3D3.integrated_PC20, file = "R3D3.integrated_PC20.Robj")
save(R3D3.integrated_PC20_subset, file = "R3D3.integrated_PC20_subset.Robj")




### Day 7 Sample Redo

# Load the dataset
R3D7.data <- Read10X(data.dir = "/home/ajaffer/NovaSeq_Aligned_Outputs/Combined_Reindeer_Bovine_Ref_Agg/Agg_Day7_Bovine/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data):
R3D7 <- CreateSeuratObject(counts = R3D7.data, project = "reindeer_RNA", min.cells = 2, min.features = 100)
R3D7


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
R3D7[["percent.mt"]] <- PercentageFeatureSet(R3D7, pattern = "^MT-")
# Visualize QC metrics as a violin plot
jpeg(file = "VlnPlot_R3D7_Redo.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
### The bovine reference genome either is not annotated with mitochondrial DNA or the reindeer sample did not align to them
### The min cutoff was not included because we wanted to ganulocytes which have low transcriptional expression
R3D7 <- subset(R3D7, subset = nFeature_RNA < 6000)
R3D7 <- subset(R3D7, subset = nCount_RNA < 65000)
jpeg(file = "VlnPlot_R3D7_Redo_After_Subset.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


write.csv(R3D7@meta.data, file = "R3D7@meta.data.csv")
### Modify CSV on Excel (-1 = back, -2 = antler)
meta_data = read.csv(file = "R3D7@meta.data.csv")
new_col = data.frame(meta_data$sample_ident)
R3D7$sample_ident = meta_data$sample_ident

R3D7.list <- SplitObject(R3D7, split.by = "sample_ident")
R3D7.list <- R3D7.list[c("Back", "Antler")]

for (i in 1:length(R3D7.list)) {
  R3D7.list[[i]] <- NormalizeData(R3D7.list[[i]], verbose = FALSE)
  R3D7.list[[i]] <- FindVariableFeatures(R3D7.list[[i]], selection.method = "vst",
                                         nfeatures = 2000, verbose = FALSE)
}

save(R3D7, file = "R3D7.Robj")
load(file = "R3D7.Robj")

R3D7.list <- R3D7.list[c("Back", "Antler")]

########### Example with D7 - dims 1:10,  ###########
R3D7.anchors_PC10 <- FindIntegrationAnchors(object.list = R3D7.list, dims = 1:10)
R3D7.integrated_PC10 <- IntegrateData(anchorset = R3D7.anchors_PC10, dims = 1:10)
DefaultAssay(R3D7.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D7.integrated_PC10 <- ScaleData(R3D7.integrated_PC10, verbose = FALSE)
R3D7.integrated_PC10 <- RunPCA(R3D7.integrated_PC10, npcs = 10, verbose = FALSE)
R3D7.integrated_PC10 <- RunUMAP(R3D7.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "R3D7.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D7 - dims 1:20,  ###########
R3D7.anchors_PC20 <- FindIntegrationAnchors(object.list = R3D7.list, dims = 1:20)
R3D7.integrated_PC20 <- IntegrateData(anchorset = R3D7.anchors_PC20, dims = 1:20)
DefaultAssay(R3D7.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D7.integrated_PC20 <- ScaleData(R3D7.integrated_PC20, verbose = FALSE)
R3D7.integrated_PC20 <- RunPCA(R3D7.integrated_PC20, npcs = 20, verbose = FALSE)
R3D7.integrated_PC20 <- RunUMAP(R3D7.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "R3D7.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D7 - dims 1:30,  ###########
R3D7.anchors_PC30 <- FindIntegrationAnchors(object.list = R3D7.list, dims = 1:30)
R3D7.integrated_PC30 <- IntegrateData(anchorset = R3D7.anchors_PC30, dims = 1:30)
DefaultAssay(R3D7.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D7.integrated_PC30 <- ScaleData(R3D7.integrated_PC30, verbose = FALSE)
R3D7.integrated_PC30 <- RunPCA(R3D7.integrated_PC30, npcs = 30, verbose = FALSE)
R3D7.integrated_PC30 <- RunUMAP(R3D7.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "R3D7.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D7 - dims 1:40,  ###########
R3D7.anchors_PC40 <- FindIntegrationAnchors(object.list = R3D7.list, dims = 1:40)
R3D7.integrated_PC40 <- IntegrateData(anchorset = R3D7.anchors_PC40, dims = 1:40)
DefaultAssay(R3D7.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D7.integrated_PC40 <- ScaleData(R3D7.integrated_PC40, verbose = FALSE)
R3D7.integrated_PC40 <- RunPCA(R3D7.integrated_PC40, npcs = 40, verbose = FALSE)
R3D7.integrated_PC40 <- RunUMAP(R3D7.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "R3D7.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D7 - dims 1:50,  ###########
R3D7.anchors_PC50 <- FindIntegrationAnchors(object.list = R3D7.list, dims = 1:50)
R3D7.integrated_PC50 <- IntegrateData(anchorset = R3D7.anchors_PC50, dims = 1:50)
DefaultAssay(R3D7.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D7.integrated_PC50 <- ScaleData(R3D7.integrated_PC50, verbose = FALSE)
R3D7.integrated_PC50 <- RunPCA(R3D7.integrated_PC50, npcs = 50, verbose = FALSE)
R3D7.integrated_PC50 <- RunUMAP(R3D7.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "R3D7.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 7 PC 50
jpeg(file = "R3D7.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(R3D7.integrated_PC50, ndims = 50)
dev.off()
# We are using PC 30 for this sample?

########### Day 7 Cell Type Annotations ###############

R3D7.integrated_PC30 <- FindNeighbors(R3D7.integrated_PC30, dims = 1:30)
R3D7.integrated_PC30 <- FindClusters(R3D7.integrated_PC30, resolution = 0.6)

write.csv(R3D7.integrated_PC30@meta.data, file = "R3D7.integrated_PC30@meta.data.csv")
Idents(R3D7.integrated_PC30) <- "seurat_clusters"
#sample_ident (it said this command was not found)

jpeg(file = "R3D7.integrated_PC30_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC30, reduction = "umap", label = T)
dev.off()
###### ARZINA REMEMBER TO SAVE THIS, save image (line 669) ######
save(R3D7.integrated_PC30, file = "R3D7.integrated_PC30.Robj")

jpeg(file = "R3D7_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT", "CRABP1"))
dev.off()

jpeg(file = "R3D7_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("KRT14", "KRT17", "KRT15"))
dev.off()

jpeg(file = "R3D7_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "R3D7_selected_markers_Immune_Cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("AIF1", "NR4A3"))
dev.off()

jpeg(file = "R3D7_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("TH", "KIT"))
dev.off()

jpeg(file = "R3D7_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "R3D7_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("CD68", "FOLR2", "FOLR3"))
dev.off()

jpeg(file = "R3D7_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "R3D7_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "R3D7_selected_markers_Real_Endo.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("PECAM1", "FLT1", "LYVE1"))
dev.off()

jpeg(file = "R3D7_selected_markers_PTPRC.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("PTPRC"))
dev.off()

jpeg(file = "R3D7_selected_markers_neutrophil.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("S100A8", "S100A9"))
dev.off()

jpeg(file = "R3D7_selected_markers_IL17_CD74.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("IL17A", "CD74"))
dev.off()

jpeg(file = "R3D7_selected_markers_Macro_Myeloid.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("CD68", "NOS2", "DUSP2", "NCF1"))
dev.off()

jpeg(file = "R3D7_selected_markers_CSFR1_NOS2.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("CSF1R", "NOS2"))
dev.off()

jpeg(file = "R3D7_selected_markers_Keratinocytes_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D7.integrated_PC30, features = c("KRT14", "KRT17", "KRT15", "MT4"))
dev.off()


Idents(R3D7.integrated_PC30) <- "seurat_clusters"
new.cluster.ids <- c("DUSP2+Macrophage", "CD68+NOS+Macrophage", "DUSP2+Macrophage", "Keratinocyte", "CD68+NOS+Macrophage", "Fibroblast", "Fibroblast", "Arterial Endothelial", "DUSP2+Macrophage", "T-cell", "CD68+CSF1R+Macrophage", "CD68+CCR7+Macrophage", "VSM", "DUSP2+Macrophage", "Fibroblast", "CD68+CSF1R+Macrophage", "Melanocyte", "Putative T-cell", "Lymphatic Endothelial", "DUSP2+Macrophage")
names(new.cluster.ids) <- levels(R3D7.integrated_PC30)
R3D7.integrated_PC30 <- RenameIdents(R3D7.integrated_PC30, new.cluster.ids)

jpeg(file = "R3D7.integrated_PC30_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC30, reduction = "umap", label = T)
dev.off()

jpeg(file = "R3D7.integrated_PC30_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC30, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

save(R3D7.integrated_PC30, file = "R3D7.integrated_PC30.Robj")

jpeg(file = "R3D7.integrated_PC30_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D7.integrated_PC30, reduction = "umap", label = F)
dev.off()



### Day 14 Sample Redo

# Load the dataset
R3D14.data <- Read10X(data.dir = "/home/ajaffer/Reindeer_Round3_FASTQs/Round3_Agg/Round3_Agg_Day14_Bovine/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data):
R3D14 <- CreateSeuratObject(counts = R3D14.data, project = "reindeer_RNA", min.cells = 2, min.features = 100)
R3D14


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
R3D14[["percent.mt"]] <- PercentageFeatureSet(R3D14, pattern = "^mt-")
# Visualize QC metrics as a violin plot
jpeg(file = "VlnPlot_R3D14_Redo.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
### The bovine reference genome either is not annotated with mitochondrial DNA or the reindeer sample did not align to them
### The min cutoff was not included because we wanted to ganulocytes which have low transcriptional expression
R3D14 <- subset(R3D14, subset = nFeature_RNA < 4500)
R3D14 <- subset(R3D14, subset = nCount_RNA < 40000)
jpeg(file = "VlnPlot_R3D14_Redo_After_Subset.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3D14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


write.csv(R3D14@meta.data, file = "R3D14@meta.data.csv")
### Modify CSV on Excel (-1 = back, -2 = antler)
meta_data = read.csv(file = "R3D14@meta.data.csv")
new_col = data.frame(meta_data$sample_ident)
R3D14$sample_ident = meta_data$sample_ident

R3D14.list <- SplitObject(R3D14, split.by = "sample_ident")
R3D14.list <- R3D14.list[c("Back", "Antler")]

for (i in 1:length(R3D14.list)) {
  R3D14.list[[i]] <- NormalizeData(R3D14.list[[i]], verbose = FALSE)
  R3D14.list[[i]] <- FindVariableFeatures(R3D14.list[[i]], selection.method = "vst",
                                          nfeatures = 2000, verbose = FALSE)
}

save(R3D14, file = "R3D14.Robj")
load(file = "R3D14.Robj")

R3D14.list <- R3D14.list[c("Back", "Antler")]

########### Example with D14 - dims 1:10,  ###########
R3D14.anchors_PC10 <- FindIntegrationAnchors(object.list = R3D14.list, dims = 1:10)
R3D14.integrated_PC10 <- IntegrateData(anchorset = R3D14.anchors_PC10, dims = 1:10)
DefaultAssay(R3D14.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D14.integrated_PC10 <- ScaleData(R3D14.integrated_PC10, verbose = FALSE)
R3D14.integrated_PC10 <- RunPCA(R3D14.integrated_PC10, npcs = 10, verbose = FALSE)
R3D14.integrated_PC10 <- RunUMAP(R3D14.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "R3D14.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D14 - dims 1:20,  ###########
R3D14.anchors_PC20 <- FindIntegrationAnchors(object.list = R3D14.list, dims = 1:20)
R3D14.integrated_PC20 <- IntegrateData(anchorset = R3D14.anchors_PC20, dims = 1:20)
DefaultAssay(R3D14.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D14.integrated_PC20 <- ScaleData(R3D14.integrated_PC20, verbose = FALSE)
R3D14.integrated_PC20 <- RunPCA(R3D14.integrated_PC20, npcs = 20, verbose = FALSE)
R3D14.integrated_PC20 <- RunUMAP(R3D14.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "R3D14.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D14 - dims 1:30,  ###########
R3D14.anchors_PC30 <- FindIntegrationAnchors(object.list = R3D14.list, dims = 1:30)
R3D14.integrated_PC30 <- IntegrateData(anchorset = R3D14.anchors_PC30, dims = 1:30)
DefaultAssay(R3D14.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D14.integrated_PC30 <- ScaleData(R3D14.integrated_PC30, verbose = FALSE)
R3D14.integrated_PC30 <- RunPCA(R3D14.integrated_PC30, npcs = 30, verbose = FALSE)
R3D14.integrated_PC30 <- RunUMAP(R3D14.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "R3D14.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()


########### Example with D14 - dims 1:40,  ###########
R3D14.anchors_PC40 <- FindIntegrationAnchors(object.list = R3D14.list, dims = 1:40)
R3D14.integrated_PC40 <- IntegrateData(anchorset = R3D14.anchors_PC40, dims = 1:40)
DefaultAssay(R3D14.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D14.integrated_PC40 <- ScaleData(R3D14.integrated_PC40, verbose = FALSE)
R3D14.integrated_PC40 <- RunPCA(R3D14.integrated_PC40, npcs = 40, verbose = FALSE)
R3D14.integrated_PC40 <- RunUMAP(R3D14.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "R3D14.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D14 - dims 1:50,  ###########
R3D14.anchors_PC50 <- FindIntegrationAnchors(object.list = R3D14.list, dims = 1:50)
R3D14.integrated_PC50 <- IntegrateData(anchorset = R3D14.anchors_PC50, dims = 1:50)
DefaultAssay(R3D14.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
R3D14.integrated_PC50 <- ScaleData(R3D14.integrated_PC50, verbose = FALSE)
R3D14.integrated_PC50 <- RunPCA(R3D14.integrated_PC50, npcs = 50, verbose = FALSE)
R3D14.integrated_PC50 <- RunUMAP(R3D14.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "R3D14.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 14 PC 50
jpeg(file = "R3D14.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(R3D14.integrated_PC50, ndims = 50)
dev.off()
# We are using PC 20 for this sample?

########### Day 14 Cell Type Annotations ###############

R3D14.integrated_PC20 <- FindNeighbors(R3D14.integrated_PC20, dims = 1:20)
R3D14.integrated_PC20 <- FindClusters(R3D14.integrated_PC20, resolution = 0.6)

write.csv(R3D14.integrated_PC20@meta.data, file = "R3D14.integrated_PC20@meta.data.csv")
Idents(R3D14.integrated_PC20) <- "seurat_clusters"
#sample_ident (it said this command was not found)

jpeg(file = "R3D14.integrated_PC20_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC20, reduction = "umap", label = T)
dev.off()
###### ARZINA REMEMBER TO SAVE THIS, save image (line 669) ######
save(R3D14.integrated_PC20, file = "R3D14.integrated_PC20.Robj")

jpeg(file = "R3D14_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                                "CRABP1"))
dev.off()

jpeg(file = "R3D14_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("KRT14", "KRT17", "KRT15"))
dev.off()

jpeg(file = "R3D14_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "R3D14_selected_markers_Immune_Cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("AIF1", "NR4A3"))
dev.off()

jpeg(file = "R3D14_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("TH", "KIT"))
dev.off()

jpeg(file = "R3D14_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "R3D14_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("CD68", "FOLR2", "FOLR3"))
dev.off()

jpeg(file = "R3D14_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "R3D14_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "R3D14_selected_markers_Real_Endo.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("PECAM1", "FLT1", "LYVE1"))
dev.off()

jpeg(file = "R3D14_selected_markers_PTPRC.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("PTPRC"))
dev.off()

jpeg(file = "R3D14_selected_markers_neutrophil.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("S100A8", "S100A9"))
dev.off()

jpeg(file = "R3D14_selected_markers_IL17_CD74.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("IL17A", "CD74"))
dev.off()

jpeg(file = "R3D14_selected_markers_SRGN.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("SRGN", "SERPINA1", "SAT1"), min.cutoff = "q10")
dev.off()

jpeg(file = "R3D14_selected_markers_B_cell.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("VPREB1", "ENSBTAG00000017305", "JCHAIN"), min.cutoff = "q10")
dev.off()

jpeg(file = "R3D14_selected_markers_Macro_Myeloid.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("CD68", "NOS2", "DUSP2", "NCF1"))
dev.off()

jpeg(file = "R3D14_selected_markers_CSF1R.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("CSF1R"))
dev.off()


R3D14.integrated_PC20_9 <- FindMarkers(R3D14.integrated_PC20, ident.1 = 9, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
write.csv(R3D14.integrated_PC20_9, file = "R3D14.integrated_PC20_9.csv")

jpeg(file = "R3D14_selected_markers_Keratinocytes_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(R3D14.integrated_PC20, features = c("KRT14", "KRT17", "KRT15", "MT4"))
dev.off()

Idents(R3D14.integrated_PC20) <- "seurat_clusters"
new.cluster.ids <- c("Basal Keratinocyte", "Arterial Endothelial", "Basal Keratinocyte", "Basal Keratinocyte", "CD68+NOS+Macrophage", "T-cell", "Basal Keratinocyte", "Basal Keratinocyte", "Fibroblast", "Suprabasal Keratinocyte", "Fibroblast", "CD68+CSF1R+Macrophage", "CD68+CSF1R+Macrophage", "Melanocyte", "Fibroblast", "SRGN+UK", "VSM", "Lymphatic Endothelial", "Arterial Endothelial", "B-cell")
names(new.cluster.ids) <- levels(R3D14.integrated_PC20)
R3D14.integrated_PC20 <- RenameIdents(R3D14.integrated_PC20, new.cluster.ids)

jpeg(file = "R3D14.integrated_PC20_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC20, reduction = "umap", label = T)
dev.off()

jpeg(file = "R3D14.integrated_PC20_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC20, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

save(R3D14.integrated_PC20, file = "R3D14.integrated_PC20.Robj")

jpeg(file = "R3D14.integrated_PC20_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(R3D14.integrated_PC20, reduction = "umap", label = F)
dev.off()


### All DAYS Sample Redo

# Load the dataset
R3A.data <- Read10X(data.dir = "/home/ajaffer/Reindeer_All_Days_Agg/Round3_Agg_All_Days_Bovine/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data):
R3A <- CreateSeuratObject(counts = R3A.data, project = "reindeer_RNA", min.cells = 2, min.features = 100)
R3A


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
R3A[["percent.mt"]] <- PercentageFeatureSet(R3A, pattern = "^MT-")
# Visualize QC metrics as a violin plot
jpeg(file = "VlnPlot_R3A_Redo.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
### The bovine reference genome either is not annotated with mitochondrial DNA or the reindeer sample did not align to them
### The min cutoff was not included because we wanted to ganulocytes which have low transcriptional expression
R3A <- subset(R3A, subset = nFeature_RNA < 6000)
R3A <- subset(R3A, subset = nCount_RNA < 50000)
jpeg(file = "VlnPlot_R3A_Redo_After_Subset.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(R3A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

write.csv(R3A@meta.data, file = "R3A@meta.data.csv")
### Modify CSV on Excel (-1 = D0B, -2 = D0A, -3 = D3B, -4 = D3A, -5 = D7B, -6 = D7A, -7 = D14B, -8 = D14A)

######## BY DAY

### Modify CSV on Excel
meta_data = read.csv(file = "R3A_day@meta.data.csv") ####CHANGE THIS
new_col = data.frame(meta_data$sample_ident)
reindeer_all_days_by_day = R3A
reindeer_all_days_by_day$sample_ident = meta_data$sample_ident
save(reindeer_all_days_by_day, file = "reindeer_all_days_by_day.Robj")


reindeer_all_days_by_day.list <- SplitObject(reindeer_all_days_by_day, split.by = "sample_ident")
reindeer_all_days_by_day.list <- reindeer_all_days_by_day.list[c("Day0", "Day3", "Day7", "Day14")]

for (i in 1:length(reindeer_all_days_by_day.list)) {
  reindeer_all_days_by_day.list[[i]] <- NormalizeData(reindeer_all_days_by_day.list[[i]], verbose = FALSE)
  reindeer_all_days_by_day.list[[i]] <- FindVariableFeatures(reindeer_all_days_by_day.list[[i]], selection.method = "vst",
                                                             nfeatures = 2000, verbose = FALSE)
}
save(reindeer_all_days_by_day, file = "reindeer_all_days_by_day.Robj")
load(file = "reindeer_all_days_by_day.Robj")

########### Example with D0 - dims 1:10,  ###########
reindeer_all_days_by_day.list_PC10 <- reindeer_all_days_by_day.list[c("Day0", "Day3", "Day7", "Day14")]
reindeer_all_days_by_day.anchors_PC10 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_day.list_PC10, dims = 1:10)
reindeer_all_days_by_day.integrated_PC10 <- IntegrateData(anchorset = reindeer_all_days_by_day.anchors_PC10, dims = 1:10)
DefaultAssay(reindeer_all_days_by_day.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_day.integrated_PC10 <- ScaleData(reindeer_all_days_by_day.integrated_PC10, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC10 <- RunPCA(reindeer_all_days_by_day.integrated_PC10, npcs = 10, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC10 <- RunUMAP(reindeer_all_days_by_day.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "reindeer_all_days_by_day.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:20,  ###########
reindeer_all_days_by_day.list_PC20 <- reindeer_all_days_by_day.list[c("Day0", "Day3", "Day7", "Day14")]
reindeer_all_days_by_day.anchors_PC20 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_day.list_PC20, dims = 1:20)
reindeer_all_days_by_day.integrated_PC20 <- IntegrateData(anchorset = reindeer_all_days_by_day.anchors_PC20, dims = 1:20)
DefaultAssay(reindeer_all_days_by_day.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_day.integrated_PC20 <- ScaleData(reindeer_all_days_by_day.integrated_PC20, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC20 <- RunPCA(reindeer_all_days_by_day.integrated_PC20, npcs = 20, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC20 <- RunUMAP(reindeer_all_days_by_day.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "reindeer_all_days_by_day.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:30,  ###########
reindeer_all_days_by_day.list_PC30 <- reindeer_all_days_by_day.list[c("Day0", "Day3", "Day7", "Day14")]
reindeer_all_days_by_day.anchors_PC30 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_day.list_PC30, dims = 1:30)
reindeer_all_days_by_day.integrated_PC30 <- IntegrateData(anchorset = reindeer_all_days_by_day.anchors_PC30, dims = 1:30)
DefaultAssay(reindeer_all_days_by_day.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_day.integrated_PC30 <- ScaleData(reindeer_all_days_by_day.integrated_PC30, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC30 <- RunPCA(reindeer_all_days_by_day.integrated_PC30, npcs = 30, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC30 <- RunUMAP(reindeer_all_days_by_day.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "reindeer_all_days_by_day.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:40,  ###########
reindeer_all_days_by_day.list_PC40 <- reindeer_all_days_by_day.list[c("Day0", "Day3", "Day7", "Day14")]
reindeer_all_days_by_day.anchors_PC40 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_day.list_PC40, dims = 1:40)
reindeer_all_days_by_day.integrated_PC40 <- IntegrateData(anchorset = reindeer_all_days_by_day.anchors_PC40, dims = 1:40)
DefaultAssay(reindeer_all_days_by_day.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_day.integrated_PC40 <- ScaleData(reindeer_all_days_by_day.integrated_PC40, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC40 <- RunPCA(reindeer_all_days_by_day.integrated_PC40, npcs = 40, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC40 <- RunUMAP(reindeer_all_days_by_day.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "reindeer_all_days_by_day.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:50,  ###########
reindeer_all_days_by_day.list_PC50 <- reindeer_all_days_by_day.list[c("Day0", "Day3", "Day7", "Day14")]
reindeer_all_days_by_day.anchors_PC50 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_day.list_PC50, dims = 1:50)
reindeer_all_days_by_day.integrated_PC50 <- IntegrateData(anchorset = reindeer_all_days_by_day.anchors_PC50, dims = 1:50)
DefaultAssay(reindeer_all_days_by_day.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_day.integrated_PC50 <- ScaleData(reindeer_all_days_by_day.integrated_PC50, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC50 <- RunPCA(reindeer_all_days_by_day.integrated_PC50, npcs = 50, verbose = FALSE)
reindeer_all_days_by_day.integrated_PC50 <- RunUMAP(reindeer_all_days_by_day.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "reindeer_all_days_by_day.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 0 PC 50
jpeg(file = "reindeer_all_days_by_day.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(reindeer_all_days_by_day.integrated_PC50, ndims = 50)
dev.off()

########### All Days By Day Cell Type Annotations ###############

reindeer_all_days_by_day.integrated_PC30 <- FindNeighbors(reindeer_all_days_by_day.integrated_PC30, dims = 1:30)
reindeer_all_days_by_day.integrated_PC30 <- FindClusters(reindeer_all_days_by_day.integrated_PC30, resolution = 0.6)

write.csv(reindeer_all_days_by_day.integrated_PC30@meta.data, file = "reindeer_all_days_by_day.integrated_PC30@meta.data.csv")
Idents(reindeer_all_days_by_day.integrated_PC30) <- "seurat_clusters"
#sample_ident (it said this command was not found)

jpeg(file = "reindeer_all_days_by_day.integrated_PC30_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC30, reduction = "umap", label = T)
dev.off()
###### ARZINA REMEMBER TO SAVE THIS, save image (line 669) ######
save(reindeer_all_days_by_day.integrated_PC30, file = "reindeer_all_days_by_day.integrated_PC30.Robj")

jpeg(file = "R3A_day_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                                                   "CRABP1"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("KRT14", "KRT17", "KRT15"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Immune_Cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("AIF1", "NR4A3"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("TH", "KIT"))
dev.off()

jpeg(file = "R3A_day_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("CD68", "FOLR2", "FOLR3"))
dev.off()

jpeg(file = "R3A_day_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Real_Endo.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("PECAM1", "FLT1", "LYVE1"))
dev.off()

jpeg(file = "R3A_day_selected_markers_PTPRC.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("PTPRC"))
dev.off()

jpeg(file = "R3A_day_selected_markers_neutrophil.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("S100A8", "S100A9"))
dev.off()

jpeg(file = "R3A_day_selected_markers_IL17_CD74.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("IL17A", "CD74"))
dev.off()

jpeg(file = "R3A_day_selected_markers_SRGN.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("SRGN", "SERPINA1", "SAT1"), min.cutoff = "q10")
dev.off()

jpeg(file = "R3A_day_selected_markers_B_cell.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("VPREB1", "ENSBTAG00000017305", "JCHAIN"), min.cutoff = "q10")
dev.off()

jpeg(file = "R3A_day_selected_markers_Macro_Myeloid.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("CD68", "NOS2", "DUSP2", "NCF1"))
dev.off()

jpeg(file = "R3A_day_selected_markers_Keratinocytes_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("KRT14", "KRT17", "KRT15", "MT4"))
dev.off()

jpeg(file = "R3A_day_selected_markers_CSF1R.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("CSF1R"))
dev.off()

reindeer_all_days_by_day_17 <- FindMarkers(reindeer_all_days_by_day.integrated_PC30, ident.1 = 17, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
write.csv(reindeer_all_days_by_day_17, file = "reindeer_all_days_by_day_17.csv")

reindeer_all_days_by_day_22 <- FindMarkers(reindeer_all_days_by_day.integrated_PC30, ident.1 = 22, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
write.csv(reindeer_all_days_by_day_22, file = "reindeer_all_days_by_day_22.csv")

jpeg(file = "R3A_day_selected_markers_SERPINA1.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_day.integrated_PC30, features = c("SERPINA1", "SERPINB2", "PTPRC", "CD74", "CCR7"))
dev.off()


Idents(reindeer_all_days_by_day.integrated_PC30) <- "seurat_clusters"
new.cluster.ids <- c("DUSP2+ Macrophage", "CD68+NOS+Macrophage", "Arterial Endothelial", "Fibroblast", "Basal Keratinocyte", "Suprabasal Keratinocyte", "Basal Keratinocyte", "Basal Keratinocyte", "CD68+CSF1R+Macrophage", "Fibroblast", "T-cell", "Basal Keratinocyte", "Suprabasal Keratinocyte", "Melanocyte", "Basal Keratinocyte", "Fibroblast", "Fibroblast", "CD74+T-cell", "VSM", "Lymphatic Endothelial", "Suprabasal Keratinocyte", "Arterial Endothelial", "CD8+T-cell", "Schwann", "CD68+CSF1R+Macrophage", "Fibroblast", "Suprabasal Keratinocyte", "T-cell")
names(new.cluster.ids) <- levels(reindeer_all_days_by_day.integrated_PC30)
reindeer_all_days_by_day.integrated_PC30 <- RenameIdents(reindeer_all_days_by_day.integrated_PC30, new.cluster.ids)

jpeg(file = "reindeer_all_days_by_day.integrated_PC30_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC30, reduction = "umap", label = T)
dev.off()

jpeg(file = "reindeer_all_days_by_day.integrated_PC30_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC30, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

jpeg(file = "reindeer_all_days_by_day.integrated_PC30_annotated_nolabels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_day.integrated_PC30, reduction = "umap", label = F)
dev.off()

save(reindeer_all_days_by_day.integrated_PC30, file = "reindeer_all_days_by_day.integrated_PC30.Robj")


write.csv(reindeer_all_days_by_day.integrated_PC30@meta.data, file = "reindeer_all_days_by_day.integrated_PC30@meta.data.csv")
### Modify CSV on Excel (-1 = D0B, -2 = D0A, -3 = D3B, -4 = D3A, -5 = D7B, -6 = D7A, -7 = D14B, -8 = D14A)

### Modify CSV on Excel
meta_data = read.csv(file = "reindeer_all_days_by_day.integrated_PC30@meta.data.csv") ####CHANGE THIS
new_col = data.frame(meta_data$condition)
reindeer_all_days_by_day.integrated_PC30$condition = meta_data$condition
save(reindeer_all_days_by_day.integrated_PC30, file = "reindeer_all_days_by_day.integrated_PC30.Robj")


meta_data = read.csv(file = "reindeer_all_days_by_day.integrated_PC30@meta.data.csv") ####CHANGE THIS
new_col = data.frame(meta_data$day_condition)
reindeer_all_days_by_day.integrated_PC30$day_condition = meta_data$day_condition
save(reindeer_all_days_by_day.integrated_PC30, file = "reindeer_all_days_by_day.integrated_PC30.Robj")



###### BY CONDITION

### Modify CSV on Excel
meta_data = read.csv(file = "R3A_condition@meta.data.csv") ###CHANGE THIS
new_col = data.frame(meta_data$sample_ident)
reindeer_all_days_by_condition = R3A
reindeer_all_days_by_condition$sample_ident = meta_data$sample_ident
save(reindeer_all_days_by_condition, file = "reindeer_all_days_by_condition.Robj")


reindeer_all_days_by_condition.list <- SplitObject(reindeer_all_days_by_condition, split.by = "sample_ident")
reindeer_all_days_by_condition.list <- reindeer_all_days_by_condition.list[c("Back", "Antler")]

for (i in 1:length(reindeer_all_days_by_condition.list)) {
  reindeer_all_days_by_condition.list[[i]] <- NormalizeData(reindeer_all_days_by_condition.list[[i]], verbose = FALSE)
  reindeer_all_days_by_condition.list[[i]] <- FindVariableFeatures(reindeer_all_days_by_condition.list[[i]], selection.method = "vst",
                                                                   nfeatures = 2000, verbose = FALSE)
}


####START FROM HERE
save(reindeer_all_days_by_condition, file = "reindeer_all_days_by_condition.Robj")
load(file = "reindeer_all_days_by_condition.Robj")


########### Example with D0 - dims 1:10,  ###########
reindeer_all_days_by_condition.list_PC10 <- reindeer_all_days_by_condition.list[c("Back", "Antler")]
reindeer_all_days_by_condition.anchors_PC10 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_condition.list_PC10, dims = 1:10)
reindeer_all_days_by_condition.integrated_PC10 <- IntegrateData(anchorset = reindeer_all_days_by_condition.anchors_PC10, dims = 1:10)
DefaultAssay(reindeer_all_days_by_condition.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_condition.integrated_PC10 <- ScaleData(reindeer_all_days_by_condition.integrated_PC10, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC10 <- RunPCA(reindeer_all_days_by_condition.integrated_PC10, npcs = 10, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC10 <- RunUMAP(reindeer_all_days_by_condition.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:20,  ###########
reindeer_all_days_by_condition.list_PC20 <- reindeer_all_days_by_condition.list[c("Back", "Antler")]
reindeer_all_days_by_condition.anchors_PC20 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_condition.list_PC20, dims = 1:20)
reindeer_all_days_by_condition.integrated_PC20 <- IntegrateData(anchorset = reindeer_all_days_by_condition.anchors_PC20, dims = 1:20)
DefaultAssay(reindeer_all_days_by_condition.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_condition.integrated_PC20 <- ScaleData(reindeer_all_days_by_condition.integrated_PC20, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC20 <- RunPCA(reindeer_all_days_by_condition.integrated_PC20, npcs = 20, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC20 <- RunUMAP(reindeer_all_days_by_condition.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:30,  ###########
reindeer_all_days_by_condition.list_PC30 <- reindeer_all_days_by_condition.list[c("Back", "Antler")]
reindeer_all_days_by_condition.anchors_PC30 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_condition.list_PC30, dims = 1:30)
reindeer_all_days_by_condition.integrated_PC30 <- IntegrateData(anchorset = reindeer_all_days_by_condition.anchors_PC30, dims = 1:30)
DefaultAssay(reindeer_all_days_by_condition.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_condition.integrated_PC30 <- ScaleData(reindeer_all_days_by_condition.integrated_PC30, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC30 <- RunPCA(reindeer_all_days_by_condition.integrated_PC30, npcs = 30, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC30 <- RunUMAP(reindeer_all_days_by_condition.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:40,  ###########
reindeer_all_days_by_condition.list_PC40 <- reindeer_all_days_by_condition.list[c("Back", "Antler")]
reindeer_all_days_by_condition.anchors_PC40 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_condition.list_PC40, dims = 1:40)
reindeer_all_days_by_condition.integrated_PC40 <- IntegrateData(anchorset = reindeer_all_days_by_condition.anchors_PC40, dims = 1:40)
DefaultAssay(reindeer_all_days_by_condition.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_condition.integrated_PC40 <- ScaleData(reindeer_all_days_by_condition.integrated_PC40, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC40 <- RunPCA(reindeer_all_days_by_condition.integrated_PC40, npcs = 40, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC40 <- RunUMAP(reindeer_all_days_by_condition.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:50,  ###########
reindeer_all_days_by_condition.list_PC50 <- reindeer_all_days_by_condition.list[c("Back", "Antler")]
reindeer_all_days_by_condition.anchors_PC50 <- FindIntegrationAnchors(object.list = reindeer_all_days_by_condition.list_PC50, dims = 1:50)
reindeer_all_days_by_condition.integrated_PC50 <- IntegrateData(anchorset = reindeer_all_days_by_condition.anchors_PC50, dims = 1:50)
DefaultAssay(reindeer_all_days_by_condition.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_all_days_by_condition.integrated_PC50 <- ScaleData(reindeer_all_days_by_condition.integrated_PC50, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC50 <- RunPCA(reindeer_all_days_by_condition.integrated_PC50, npcs = 50, verbose = FALSE)
reindeer_all_days_by_condition.integrated_PC50 <- RunUMAP(reindeer_all_days_by_condition.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 0 PC 50
jpeg(file = "reindeer_all_days_by_condition.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(reindeer_all_days_by_condition.integrated_PC50, ndims = 50)
dev.off()


########### All Days By Day Cell Type Annotations ###############

reindeer_all_days_by_condition.integrated_PC30 <- FindNeighbors(reindeer_all_days_by_condition.integrated_PC30, dims = 1:30)
reindeer_all_days_by_condition.integrated_PC30 <- FindClusters(reindeer_all_days_by_condition.integrated_PC30, resolution = 0.6)

write.csv(reindeer_all_days_by_condition.integrated_PC30@meta.data, file = "reindeer_all_days_by_condition.integrated_PC30@meta.data.csv")
Idents(reindeer_all_days_by_condition.integrated_PC30) <- "seurat_clusters"
#sample_ident (it said this command was not found)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", label = T)
dev.off()
###### ARZINA REMEMBER TO SAVE THIS, save image (line 669) ######
save(reindeer_all_days_by_condition.integrated_PC30, file = "reindeer_all_days_by_condition.integrated_PC30.Robj")

jpeg(file = "R3A_condition_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                                                         "CRABP1"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Keratinocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("KRT14", "KRT17", "KRT15"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Endothelial.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("PECAM1", "ACKR1", "PTN", "APOA1"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Immune_Cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("AIF1", "NR4A3"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Melanocytes.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("TH", "KIT"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_T_cells.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("CD3E", "CD8A", "CD3D", "CD8B", "CD7", "CCR7", "IFNG", "CD52"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Macrophages.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("CD68", "FOLR2", "FOLR3"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_VSM.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("DES", "VIM4", "VIM", "PLN"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Schwaan.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("SOX10", "MBP"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Real_Endo.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("PECAM1", "FLT1", "LYVE1"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_PTPRC.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("PTPRC"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_neutrophil.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("S100A8", "S100A9"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_IL17_CD74.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("IL17A", "CD74"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_SRGN.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("SRGN", "SERPINA1", "SAT1"), min.cutoff = "q10")
dev.off()

jpeg(file = "R3A_condition_selected_markers_B_cell.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("VPREB1", "ENSBTAG00000017305", "JCHAIN"), min.cutoff = "q10")
dev.off()

jpeg(file = "R3A_condition_selected_markers_Macro_Myeloid.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("CD68", "NOS2", "DUSP2", "NCF1"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_Keratinocytes_MT4.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("KRT14", "KRT17", "KRT15", "MT4"))
dev.off()

jpeg(file = "R3A_condition_selected_markers_CSF1R.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_all_days_by_condition.integrated_PC30, features = c("CSF1R"))
dev.off()


reindeer_all_days_by_condition.integrated_PC30_26 <- FindMarkers(reindeer_all_days_by_condition.integrated_PC30, ident.1 = 26, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
write.csv(reindeer_all_days_by_condition.integrated_PC30_26, file = "reindeer_all_days_by_condition.integrated_PC30_26.csv")


jpeg(file = "R3A_condition_selected_markers_CCR7_Vln.jpeg", width = 30, height = 25, units = "cm", res = 500)
VlnPlot(reindeer_all_days_by_condition.integrated_PC30, features = c("CCR7"))
dev.off()

reindeer_all_days_by_condition.integrated_PC30_10 <- FindMarkers(reindeer_all_days_by_condition.integrated_PC30, ident.1 = 10, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
write.csv(reindeer_all_days_by_condition.integrated_PC30_10, file = "reindeer_all_days_by_condition.integrated_PC30_10.csv")

jpeg(file = "R3A_condition_selected_markers_CD45_Vln.jpeg", width = 30, height = 25, units = "cm", res = 500)
VlnPlot(reindeer_all_days_by_condition.integrated_PC30, features = c("CD45"))
dev.off()


Idents(reindeer_all_days_by_condition.integrated_PC30) <- "seurat_clusters"
new.cluster.ids <- c("DUSP2+ Macrophage", "Arterial Endothelial", "CD68+CSF1R+Macrophage", "CD68+NOS+Macrophage", "Basal Keratinocyte", "Fibroblast", "Suprabasal Keratinocyte", "Basal Keratinocyte", "Basal Keratinocyte", "T-cell", "CD45+Leukocyte", "Fibroblast", "Basal Keratinocyte", "Suprabasal Keratinocyte", "Fibroblast", "Melanocyte", "Fibroblast", "CD68+CSF1R+Macrophage", "VSM","Lymphatic Endothelial", "Suprabasal Keratinocyte", "Basal Keratinocyte", "B-cell", "Schwann", "CD8B+T-cell", "CD68+CSF1R+Macrophage", "CD68+CSF1R+Macrophage")
names(new.cluster.ids) <- levels(reindeer_all_days_by_condition.integrated_PC30)
reindeer_all_days_by_condition.integrated_PC30 <- RenameIdents(reindeer_all_days_by_condition.integrated_PC30, new.cluster.ids)

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", label = T)
dev.off()

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30_annotated_nolabels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", label = F)
dev.off()

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30_annotated_bycondition.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", label = F, group.by = "sample_ident")
dev.off()

jpeg(file = "reindeer_all_days_by_condition.integrated_PC30_annotated_byday.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_all_days_by_condition.integrated_PC30, reduction = "umap", label = F, group.by = "day")
dev.off()

save(reindeer_all_days_by_condition.integrated_PC30, file = "reindeer_all_days_by_condition.integrated_PC30.Robj")


save.image(file = "Reindeer_May8_2020_Analysis.RData")


write.csv(reindeer_all_days_by_condition.integrated_PC30@meta.data, file = "reindeer_all_days_by_condition.integrated_PC30@meta.data.csv")
### Modify CSV on Excel (-1 = D0B, -2 = D0A, -3 = D3B, -4 = D3A, -5 = D7B, -6 = D7A, -7 = D14B, -8 = D14A)

### Modify CSV on Excel
meta_data = read.csv(file = "reindeer_all_days_by_condition.integrated_PC30@meta.data.csv") ####CHANGE THIS
new_col = data.frame(meta_data$day)
reindeer_all_days_by_condition.integrated_PC30$day = meta_data$day
save(reindeer_all_days_by_condition.integrated_PC30, file = "reindeer_all_days_by_condition.integrated_PC30.Robj")

meta_data = read.csv(file = "reindeer_all_days_by_condition.integrated_PC30@meta.data.csv") ####CHANGE THIS
new_col = data.frame(meta_data$day_condition)
reindeer_all_days_by_condition.integrated_PC30$day_condition = meta_data$day_condition
save(reindeer_all_days_by_condition.integrated_PC30, file = "reindeer_all_days_by_condition.integrated_PC30.Robj")


save.image(file = "Reindeer_May16_2020_Analysis.RData")
