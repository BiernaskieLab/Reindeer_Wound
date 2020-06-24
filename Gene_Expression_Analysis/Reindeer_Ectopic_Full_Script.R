#### Reindeer Single-Cell Ectopic Analysis June 19

tmux
bsub -Is bash

export PATH="/home/ajaffer/anaconda3/bin:$PATH"
conda env list
source activate Seurat_v3.1
cd /home/ajaffer/Reindeer_Single_Cell_Analysis_Ectopic/
  R
getwd()
sessionInfo()

library(dplyr)
library(Seurat)
library(patchwork)

### Ectopic Analysis

# Load the dataset
E.data <- Read10X(data.dir = "/home/ajaffer/Reindeer_Ectopic_Back_scRNA/primary_alignment/Aggr/Round3_Agg_Ectopic_Bovine/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data):
E <- CreateSeuratObject(counts = E.data, project = "reindeer_RNA", min.cells = 2, min.features = 100)
E


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
E[["percent.mt"]] <- PercentageFeatureSet(E, pattern = "^MT-")
# Visualize QC metrics as a violin plot
jpeg(file = "VlnPlot_E_Filtering.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(E, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
### The bovine reference genome either is not annotated with mitochondrial DNA or the reindeer sample did not align to them
### The min cutoff was not included because we wanted to ganulocytes which have low transcriptional expression
E <- subset(E, subset = nFeature_RNA < 4500)
E <- subset(E, subset = nCount_RNA < 45000)
jpeg(file = "VlnPlot_E_Filtering_After_Subset.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(E, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

write.csv(E@meta.data, file = "E@meta.data.csv")
### Modify CSV on Excel (-1 = D0B, -2 = D0A, -3 = D0E)

### Modify CSV on Excel
meta_data = read.csv(file = "E@meta.data.csv") ####CHANGE THIS
new_col = data.frame(meta_data$sample_ident)
reindeer_ectopic = E
reindeer_ectopic$sample_ident = meta_data$sample_ident
save(reindeer_ectopic, file = "reindeer_ectopic.Robj")


reindeer_ectopic.list <- SplitObject(reindeer_ectopic, split.by = "sample_ident")
reindeer_ectopic.list <- reindeer_ectopic.list[c("Day0_Back", "Day0_Antler", "Day0_Ectopic")]

for (i in 1:length(reindeer_ectopic.list)) {
  reindeer_ectopic.list[[i]] <- NormalizeData(reindeer_ectopic.list[[i]], verbose = FALSE)
  reindeer_ectopic.list[[i]] <- FindVariableFeatures(reindeer_ectopic.list[[i]], selection.method = "vst",
                                                     nfeatures = 2000, verbose = FALSE)
}
save(reindeer_ectopic, file = "reindeer_ectopic.Robj")
load(file = "reindeer_ectopic.Robj")

########### Example with D0 - dims 1:10,  ###########
reindeer_ectopic.list_PC10 <- reindeer_ectopic.list[c("Day0_Back", "Day0_Antler", "Day0_Ectopic")]
reindeer_ectopic.anchors_PC10 <- FindIntegrationAnchors(object.list = reindeer_ectopic.list_PC10, dims = 1:10)
reindeer_ectopic.integrated_PC10 <- IntegrateData(anchorset = reindeer_ectopic.anchors_PC10, dims = 1:10)
DefaultAssay(reindeer_ectopic.integrated_PC10) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_ectopic.integrated_PC10 <- ScaleData(reindeer_ectopic.integrated_PC10, verbose = FALSE)
reindeer_ectopic.integrated_PC10 <- RunPCA(reindeer_ectopic.integrated_PC10, npcs = 10, verbose = FALSE)
reindeer_ectopic.integrated_PC10 <- RunUMAP(reindeer_ectopic.integrated_PC10, reduction = "pca", dims = 1:10)

jpeg(file = "reindeer_ectopic.integrated_PC10.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC10, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:20,  ###########
reindeer_ectopic.list_PC20 <- reindeer_ectopic.list[c("Day0_Back", "Day0_Antler", "Day0_Ectopic")]
reindeer_ectopic.anchors_PC20 <- FindIntegrationAnchors(object.list = reindeer_ectopic.list_PC20, dims = 1:20)
reindeer_ectopic.integrated_PC20 <- IntegrateData(anchorset = reindeer_ectopic.anchors_PC20, dims = 1:20)
DefaultAssay(reindeer_ectopic.integrated_PC20) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_ectopic.integrated_PC20 <- ScaleData(reindeer_ectopic.integrated_PC20, verbose = FALSE)
reindeer_ectopic.integrated_PC20 <- RunPCA(reindeer_ectopic.integrated_PC20, npcs = 20, verbose = FALSE)
reindeer_ectopic.integrated_PC20 <- RunUMAP(reindeer_ectopic.integrated_PC20, reduction = "pca", dims = 1:20)

jpeg(file = "reindeer_ectopic.integrated_PC20.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC20, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:30,  ###########
reindeer_ectopic.list_PC30 <- reindeer_ectopic.list[c("Day0_Back", "Day0_Antler", "Day0_Ectopic")]
reindeer_ectopic.anchors_PC30 <- FindIntegrationAnchors(object.list = reindeer_ectopic.list_PC30, dims = 1:30)
reindeer_ectopic.integrated_PC30 <- IntegrateData(anchorset = reindeer_ectopic.anchors_PC30, dims = 1:30)
DefaultAssay(reindeer_ectopic.integrated_PC30) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_ectopic.integrated_PC30 <- ScaleData(reindeer_ectopic.integrated_PC30, verbose = FALSE)
reindeer_ectopic.integrated_PC30 <- RunPCA(reindeer_ectopic.integrated_PC30, npcs = 30, verbose = FALSE)
reindeer_ectopic.integrated_PC30 <- RunUMAP(reindeer_ectopic.integrated_PC30, reduction = "pca", dims = 1:30)

jpeg(file = "reindeer_ectopic.integrated_PC30.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC30, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:40,  ###########
reindeer_ectopic.list_PC40 <- reindeer_ectopic.list[c("Day0_Back", "Day0_Antler", "Day0_Ectopic")]
reindeer_ectopic.anchors_PC40 <- FindIntegrationAnchors(object.list = reindeer_ectopic.list_PC40, dims = 1:40)
reindeer_ectopic.integrated_PC40 <- IntegrateData(anchorset = reindeer_ectopic.anchors_PC40, dims = 1:40)
DefaultAssay(reindeer_ectopic.integrated_PC40) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_ectopic.integrated_PC40 <- ScaleData(reindeer_ectopic.integrated_PC40, verbose = FALSE)
reindeer_ectopic.integrated_PC40 <- RunPCA(reindeer_ectopic.integrated_PC40, npcs = 40, verbose = FALSE)
reindeer_ectopic.integrated_PC40 <- RunUMAP(reindeer_ectopic.integrated_PC40, reduction = "pca", dims = 1:40)

jpeg(file = "reindeer_ectopic.integrated_PC40.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC40, reduction = "umap", group.by = "sample_ident")
dev.off()

########### Example with D0 - dims 1:50,  ########### I DIDNT RUN THE 50 CAUSE 30 LOOKS FINE
reindeer_ectopic.list_PC50 <- reindeer_ectopic.list[c("Day0_Back", "Day0_Antler", "Day0_Ectopic")]
reindeer_ectopic.anchors_PC50 <- FindIntegrationAnchors(object.list = reindeer_ectopic.list_PC50, dims = 1:50)
reindeer_ectopic.integrated_PC50 <- IntegrateData(anchorset = reindeer_ectopic.anchors_PC50, dims = 1:50)
DefaultAssay(reindeer_ectopic.integrated_PC50) <- "integrated"

# Run the standard workflow for visualization and clustering
reindeer_ectopic.integrated_PC50 <- ScaleData(reindeer_ectopic.integrated_PC50, verbose = FALSE)
reindeer_ectopic.integrated_PC50 <- RunPCA(reindeer_ectopic.integrated_PC50, npcs = 50, verbose = FALSE)
reindeer_ectopic.integrated_PC50 <- RunUMAP(reindeer_ectopic.integrated_PC50, reduction = "pca", dims = 1:50)

jpeg(file = "reindeer_ectopic.integrated_PC50.jpeg", width = 15, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC50, reduction = "umap", group.by = "sample_ident")
dev.off()

####ELBOW PLOT DAY 0 PC 50
jpeg(file = "reindeer_ectopic.integrated_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(reindeer_ectopic.integrated_PC40, ndims = 40)
dev.off()

########### All Days By Day Cell Type Annotations ###############

reindeer_ectopic.integrated_PC30 <- FindNeighbors(reindeer_ectopic.integrated_PC30, dims = 1:30)
reindeer_ectopic.integrated_PC30 <- FindClusters(reindeer_ectopic.integrated_PC30, resolution = 0.6)

write.csv(reindeer_ectopic.integrated_PC30@meta.data, file = "reindeer_ectopic.integrated_PC30@meta.data.csv")
Idents(reindeer_ectopic.integrated_PC30) <- "seurat_clusters"
#sample_ident (it said this command was not found)

jpeg(file = "reindeer_ectopic.integrated_PC30_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC30, reduction = "umap", label = T)
dev.off()
###### ARZINA REMEMBER TO SAVE THIS, save image (line 669) ######
save(reindeer_ectopic.integrated_PC30, file = "reindeer_ectopic.integrated_PC30.Robj")

jpeg(file = "RE_day_selected_markers_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(reindeer_ectopic.integrated_PC30, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                                           "CRABP1"))
dev.off()

reindeer_ectopic.integrated_PC30_subset_fibros = subset(reindeer_ectopic.integrated_PC30, idents = c("0","1", "13", "15", "16"))
jpeg(file = "reindeer_ectopic.integrated_PC30_subset_fibros_0.6_res_clustering.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC30_subset_fibros, reduction = "umap", label = T)
dev.off()
save(reindeer_ectopic.integrated_PC30_subset_fibros, file = "reindeer_ectopic.integrated_PC30_subset_fibros.Robj")

save.image(file = "Reindeer_Ectopic_June19_2020_Analysis.RData")

#### May 22 Full Annotations

###Fibroblast: "0","1", "13", "15", "16"

load(file = "reindeer_ectopic.integrated_PC30.Robj")

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


### IDENTITIES

Idents(reindeer_ectopic.integrated_PC30) <- "seurat_clusters"
new.cluster.ids <- c("Fibroblast", "Fibroblast", "Basal Keratinocyte", "Arterial Endothelial", "CD45+Il1B+Myeloid", "Basal Keratinocyte", "Suprabasal Keratinocyte", "Melanocyte", "Basal Keratinocyte", "T-cell", "Basal Keratinocyte", "Macrophage_Subset2", "Suprabasal Keratinocyte", "Fibroblast", "VSM", "Fibroblast", "Fibroblast", "Macrophage_Subset1", "Suprabasal Keratinocyte", "Lymphatic Endothelial", "Basal Keratinocyte", "Schwann", "Basal Keratinocyte")
names(new.cluster.ids) <- levels(reindeer_ectopic.integrated_PC30)
reindeer_ectopic.integrated_PC30 <- RenameIdents(reindeer_ectopic.integrated_PC30, new.cluster.ids)

jpeg(file = "reindeer_ectopic.integrated_PC30_annotated_labels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC30, reduction = "umap", label = T)
dev.off()

jpeg(file = "reindeer_ectopic.integrated_PC30_annotated_split_Back_Antler.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC30, reduction = "umap", label = F, split.by = "sample_ident")
dev.off()

jpeg(file = "reindeer_ectopic.integrated_PC30_annotated_nolabels.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(reindeer_ectopic.integrated_PC30, reduction = "umap", label = F)
dev.off()

write.csv(reindeer_ectopic.integrated_PC30@meta.data, "metadata_ectopic.csv")
write.csv(reindeer_ectopic.integrated_PC30@active.ident, "activeident_ectopic.csv")

meta_data = read.csv(file = "metadata_ectopic.csv") ####CHANGE THIS
new_col = data.frame(meta_data$cell_annotation)
reindeer_ectopic.integrated_PC30$cell_annotation = meta_data$cell_annotation

meta_data = read.csv(file = "metadata_ectopic.csv") ####CHANGE THIS
new_col = data.frame(meta_data$sample_celltype)
reindeer_ectopic.integrated_PC30$sample_celltype = meta_data$sample_celltype



save(reindeer_ectopic.integrated_PC30, file = "reindeer_ectopic.integrated_PC30.Robj")


