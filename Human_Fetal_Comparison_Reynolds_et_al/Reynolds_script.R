load(file = "NAME.Robj")
load(file = "fetalNAME.Robj")

library(dplyr)
library(Seurat)
library(patchwork)

NAME[["percent.mt"]] <- PercentageFeatureSet(NAME, pattern = "^MT-")
fetalNAME[["percent.mt"]] <- PercentageFeatureSet(fetalNAME, pattern = "^MT-")


jpeg(file = "VlnPlot_NAME_Filtering.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(NAME, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file = "VlnPlot_fetalNAME_Filtering.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(fetalNAME, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

NAME <- subset(NAME, subset = nFeature_RNA < 5200 & nCount_RNA < 50000 & percent.mt < 50)
fetalNAME <- subset(fetalNAME, subset = nFeature_RNA < 4700 & nCount_RNA < 25000 & percent.mt < 10)

jpeg(file = "VlnPlot_NAME_PostFiltering.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(NAME, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()




jpeg(file = "VlnPlot_fetalNAME_PostFiltering.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(fetalNAME, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

NAME <- NormalizeData(NAME, normalization.method = "LogNormalize", scale.factor = 10000)
fetalNAME <- NormalizeData(fetalNAME, normalization.method = "LogNormalize", scale.factor = 10000)

NAME <- FindVariableFeatures(NAME, selection.method = "vst", nfeatures = 2000)
fetalNAME <- FindVariableFeatures(fetalNAME, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(NAME)
NAME <- ScaleData(NAME, features = all.genes)

all.genes <- rownames(fetalNAME)
fetalNAME <- ScaleData(fetalNAME, features = all.genes)

NAME <- RunPCA(NAME, features = VariableFeatures(object = NAME))
fetalNAME <- RunPCA(fetalNAME, features = VariableFeatures(object = fetalNAME))

jpeg(file = "NAME_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(NAME, ndims = 50)
dev.off()

jpeg(file = "fetalNAME_elbow.jpeg", width = 40, height = 25, units = "cm", res = 500)
ElbowPlot(fetalNAME, ndims = 50)
dev.off()

NAME <- FindNeighbors(NAME, dims = 1:25)
NAME <- FindClusters(NAME, resolution = 0.6)

fetalNAME <- FindNeighbors(fetalNAME, dims = 1:25)
fetalNAME <- FindClusters(fetalNAME, resolution = 0.6)

NAME <- RunUMAP(NAME, dims = 1:25)
fetalNAME <- RunUMAP(fetalNAME, dims = 1:25)


jpeg(file = "NAME_dimplot.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(NAME, reduction = "umap", label = T)
dev.off()

jpeg(file = "fetalNAME_dimplot.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(fetalNAME, reduction = "umap", label = T)
dev.off()

jpeg(file = "NAME_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(NAME, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                               "CRABP1"))
dev.off()

jpeg(file = "fetalNAME_Fibro.jpeg", width = 30, height = 25, units = "cm", res = 500)
FeaturePlot(fetalNAME, features = c("COL1A1", "COL1A2", "PDGFRA", "DPT",
                                               "CRABP1"))
dev.off()

save(NAME, file = 'NAME_25.Robj')
save(fetalNAME, file = 'fetalNAME_25.Robj')


final_clustering

jpeg(file = "NAME_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(NAME, reduction = "umap", label = T, group.by = "final_clustering")
dev.off()

jpeg(file = "fetalNAME_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(fetalNAME, reduction = "umap", label = T, group.by = "anno_final")
dev.off()

load(file = 'fetalNAME_25.Robj')
load(file = 'NAME_25.Robj')

FibrofetalNAME = subset(fetalNAME, subset = anno_final == "D_fs_FB")

NAME@active.ident = NAME@meta.data$subset

FibroNAME = subset(NAME, subset = Status == "Healthy")
write.csv(FibroNAME@meta.data, file = "FibroNAME.csv")
FibroNAME = subset(FibroNAME, subset = subset == "Yes")
Yes = subset(NAME, idents = "Yes")
Yes = subset(NAME, idents = "Yes")


jpeg(file = "FibroNAME_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(FibroNAME, reduction = "umap", label = T)
dev.off()

write.csv(NAME@meta.data, file = "NAME_metadata.csv")
meta_data = read.csv(file = "NAME_metadata.csv") 
new_col = data.frame(meta_data$subset)
NAME$subset = meta_data$subset
write.csv(FibroNAME@meta.data, file = "FibronewNAME.csv")


save(FibroNAME, file = 'FibroNAME.Robj')
save(FibrofetalNAME, file = 'FibrofetalNAME.Robj')

NAME <- subset(NAME, subset = nFeature_RNA < 5200 & nCount_RNA < 50000 & percent.mt < 50)
write.csv(FibroNAME@meta.data, file = "FibronewNAME.csv")

jpeg(file = "VlnPlot_NAME_Filtering.jpeg", width = 15, height = 25, units = "cm", res = 500)
VlnPlot(NAME, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###### RESTART SUBSET

load(file = 'fetalNAME_25.Robj')

FibrofetalNAME = subset(fetalNAME, subset = anno_final == "D_fs_FB")
save(FibrofetalNAME, file = 'FibrofetalNAME.Robj')

load(file = 'NAME_25.Robj')

Healthy = subset(NAME, subset = Status == "Healthy")
write.csv(Healthy@meta.data, file = "Healthy_metadata.csv")
save(Healthy, file = 'Healthy.Robj')

meta_data = read.csv(file = "NAME_metadata.csv") 
new_col = data.frame(meta_data$subset)
NAME$subset = meta_data$subset

load(file = 'Healthy.Robj')
Healthy@active.ident = Healthy$final_clustering
Healthy_Fibros = subset(Healthy, idents = c("F1", "F2", "F3"))
write.csv(Healthy_Fibros@meta.data, file = "Healthy_Fibros_metadata.csv")
save(Healthy_Fibros, file = 'Healthy_Fibros.Robj')


load(file = 'FibrofetalNAME.Robj')
write.csv(FibrofetalNAME@meta.data, file = "FibrofetalNAME_metadata.csv")
meta_data = read.csv(file = "FibrofetalNAME_metadata.csv") 
new_col = data.frame(meta_data$stage)
FibrofetalNAME$stage = meta_data$stage
write.csv(FibrofetalNAME@meta.data, file = "FibrofetalNAME_metadata.csv")
save(FibrofetalNAME, file = 'FibrofetalNAME_new.Robj')

####Merge and Integrate the Adult and Fetal Fibroblasts

load(file = 'Healthy_Fibros.Robj')
load(file = 'FibrofetalNAME_new.Robj')

mergefibros = merge(x = Healthy_Fibros, y = FibrofetalNAME)
save(mergefibros, file = 'mergefibros.Robj')

### Integrate by adult and fetal
library(Seurat)
#library(SeuratData)
library(patchwork)

stage.list <- SplitObject(mergefibros, split.by = "stage")

stage.list <- lapply(X = stage.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = stage.list)

immune.anchors <- FindIntegrationAnchors(object.list = stage.list, anchor.features = features)

immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

jpeg(file = "merge_stage_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", group.by = "stage")
dev.off()

jpeg(file = "merge_clusters_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

jpeg(file = "merge_split_by_stage_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", split.by = "stage")
dev.off()

stage = immune.combined
save(stage, file = 'stage.Robj')

####Merge and Integrate the Adult and Fetal Fibroblasts

load(file = 'Healthy_Fibros.Robj')
load(file = 'FibrofetalNAME_new.Robj')

mergefibros = merge(x = Healthy_Fibros, y = FibrofetalNAME)
save(mergefibros, file = 'mergefibros.Robj')

### Integrate by donor id
library(Seurat)
#library(SeuratData)
library(patchwork)

donor_id.list <- SplitObject(mergefibros, split.by = "donor_id")

donor_id.list <- lapply(X = donor_id.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = donor_id.list)

immune.anchors <- FindIntegrationAnchors(object.list = donor_id.list, anchor.features = features)

immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

jpeg(file = "merge_donor_id_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", group.by = "stage")
dev.off()

jpeg(file = "merge_clusters_donor_id_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

jpeg(file = "merge_split_by_donor_id_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", split.by = "stage")
dev.off()

donor_id = immune.combined
save(donor_id, file = 'donor_id.Robj')

jpeg(file = "actual_merge_donor_id_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", group.by = "donor_id")
dev.off()

jpeg(file = "actual_merge_split_by_donor_id_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", split.by = "donor_id")
dev.off()

load(file = 'donor_id.Robj')
load(file = 'stage.Robj')

write.csv(donor_id@meta.data, file = "donor_id_metadata.csv")
write.csv(stage@meta.data, file = "stage_metadata.csv")

meta_data = read.csv(file = "donor_id_metadata.csv") 
new_col = data.frame(meta_data$enriched)
donor_id$enriched = meta_data$enriched
write.csv(donor_id@meta.data, file = "donor_id_metadata.csv")
save(donor_id, file = 'donor_id.Robj')

meta_data = read.csv(file = "stage_metadata.csv") 
new_col = data.frame(meta_data$enriched)
stage$enriched = meta_data$enriched
write.csv(stage@meta.data, file = "stage_metadata.csv")
save(stage, file = 'stage.Robj')

#### Subset Reindeer Fibros
load(file = 'R3D0.integrated_PC10.Robj')
Fibros = subset(R3D0.integrated_PC10, idents = "Fibroblast")
save(Fibros, file = 'Fibros.Robj')

####Merge
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

load(file = 'Fibros.Robj')
load(file = 'stage.Robj')
load(file = 'donor_id.Robj')

write.csv(Fibros@meta.data, file = "Fibros_metadata.csv")
meta_data = read.csv(file = "Fibros_metadata.csv") 
new_col = data.frame(meta_data$stage)
Fibros$stage = meta_data$stage
save(Fibros, file = 'Fibros.Robj')

mergestage = merge(x = Fibros, y = stage)
save(mergestage, file = 'mergestage.Robj')
mergedonor_id = merge(x = Fibros, y = donor_id)
save(mergedonor_id, file = 'mergedonor_id.Robj')



########START FRESHHHHH
load(file = 'Healthy_Fibros.Robj')
load(file = 'FibrofetalNAME_new.Robj')
load(file = 'Fibros.Robj')

AllFibrosmerge <- merge(Healthy_Fibros, y = c(FibrofetalNAME, Fibros))
save(AllFibrosmerge, file = 'AllFibrosmerge.Robj')
write.csv(AllFibrosmerge@meta.data, file = "AllFibrosmerge_metadata.csv")


library(Seurat)
library(patchwork)
load(file = 'AllFibrosmerge.Robj')
ifnb.list <- SplitObject(AllFibrosmerge, split.by = "stage")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

jpeg(file = "all_fibros_merge_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", group.by = "stage")
dev.off()

jpeg(file = "all_fibros_merge_clusters_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

jpeg(file = "all_fibros_merge_dimplot_anno_split_stage.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", split.by = "stage")
dev.off()

jpeg(file = "all_fibros_crabp1.jpeg", width = 20, height = 25, units = "cm", res = 500)
FeaturePlot(immune.combined, features = "CRABP1", split.by = "stage")
dev.off()

allfibrosintegrate = immune.combined

write.csv(allfibrosintegrate@meta.data, 'allfibrosintegrate_metadata.csv')
save(allfibrosintegrate, file = 'allfibrosintegrate.Robj')


library(Seurat)
library(patchwork)
load(file = 'AllFibrosmerge.Robj')

############### Subsampling equal number of cells across both conditions ###############
load(file = "AllFibrosmerge.Robj") #### Specific to object you want to subsample
obj.split <- SplitObject(AllFibrosmerge, split.by = 'stage') #### Specific to object you want to subsample
Back = obj.split$Back
Antler = obj.split$Antler
adult = obj.split$adult
fetal = obj.split$fetal
View(Back@meta.data)
length(rownames(Back@meta.data)) #3354
View(Antler@meta.data)
length(rownames(Antler@meta.data)) #806
View(adult@meta.data)
length(rownames(adult@meta.data)) #17418
View(fetal@meta.data)
length(rownames(fetal@meta.data)) #33365
## Pick the lowest of four

Back.downsample = subset(Back, cells = sample(Cells(Back), 806)) #lowest of the two
Antler.downsample = subset(Antler, cells = sample(Cells(Antler), 806))#lowest of the two
adult.downsample = subset(adult, cells = sample(Cells(adult), 806)) #lowest of the two
fetal.downsample = subset(fetal, cells = sample(Cells(fetal), 806))

# umap_Back = Embeddings(Back.downsample, reduction = "umap")
# umap_Antler = Embeddings(Antler.downsample, reduction = "umap")
# umap_adult = Embeddings(adult.downsample, reduction = "umap")
# umap_fetal = Embeddings(fetal.downsample, reduction = "umap")
# umap = rbind(umap_Back, umap_Antler, umap_adult, umap_fetal)
# View(umap)

AllFibrossubset <- merge(Back.downsample, y = c(Antler.downsample, adult.downsample, fetal.downsample))

# Back_Antler <- merge(Back.downsample, y = Antler.downsample, 
#                      project = "D3", merge.data = TRUE)
# Back_Antler[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(Back_Antler))

write.csv(AllFibrossubset@meta.data, file = "AllFibrossubset_metadata.csv")
save(AllFibrossubset, file = "AllFibrossubset.Robj")

# DimPlot(Back_Antler, reduction = "umap", split.by = "sample_ident")
# FeaturePlot(Back_Antler, features = "CCL2", split.by = "sample_ident", reduction = "umap", min.cutoff = "q2")
# VlnPlot(Back_Antler, features = "CCL2", split.by = "sample_ident")


save(AllFibrossubset, file = "AllFibrossubset.Robj")
########Integrate Subsetted Object
library(Seurat)
library(patchwork)
load(file = 'AllFibrossubset.Robj')
ifnb.list <- SplitObject(AllFibrossubset, split.by = "stage")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

jpeg(file = "all_fibros_subset_merge_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", group.by = "stage")
dev.off()

jpeg(file = "all_fibros_subset_merge_clusters_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

jpeg(file = "all_fibros_subset_merge_dimplot_anno_split_stage.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", split.by = "stage")
dev.off()

jpeg(file = "all_fibros_subset_crabp1.jpeg", width = 40, height = 25, units = "cm", res = 500)
FeaturePlot(immune.combined, features = "CRABP1", split.by = "stage")
dev.off()

allfibrossubsetintegrate = immune.combined

write.csv(allfibrossubsetintegrate@meta.data, 'allfibrossubsetintegrate_metadata.csv')
save(allfibrossubsetintegrate, file = 'allfibrossubsetintegrate.Robj')



####Attempt 2 Reference Based Integration
library(ggplot2)
library(cowplot)
library(patchwork)

load(file = 'Fibros.Robj')
load(file = 'stage.Robj')
load(file = 'donor_id.Robj')

pancreas.integrated <- stage
pancreas.query <- Fibros
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$stage, 
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query, 
                           refdata = list(celltype = "stage"), reference.reduction = "pca", reduction.model = "umap")

jpeg(file = "reference_integrated_stage.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(pancreas.integrated, reduction = "umap", group.by = "stage", label = TRUE, label.size = 3, 
        repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
dev.off()

jpeg(file = "query_integrated_stage.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, 
        label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
dev.off()

jpeg(file = "query_integrated_stage_split.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", split.by = "stage", label = TRUE, 
        label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
dev.off()

save(pancreas.query, file = "pancreas.query.Robj")
save(pancreas.integrated, file = "pancreas.integrated.Robj")

pancreas.query$umap = pancreas.query$ref.umap


integratedmerge <- merge(pancreas.integrated, y = pancreas.query)


jpeg(file = "integratedmerge_integrated_stage.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(integratedmerge, reduction = "umap", group.by = "stage", label = TRUE, 
        label.size = 3, repel = TRUE) + NoLegend() + ggtitle("integratedmerge")
dev.off()

write.csv(pancreas.query@meta.data, file = "pancreas.query_metadata.csv")


##### Integrate with Subset Object

load(file = "allfibrossubsetintegrate.Robj")

pancreas.list <- SplitObject(allfibrossubsetintegrate, split.by = "stage")
pancreas.list <- pancreas.list[c("adult", "fetal", "Back", "Antler")]

adult <- NormalizeData(pancreas.list[[i]], verbose = FALSE)


for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 1000, 
                                             verbose = FALSE)
}

####Merge Big Object


load(file = 'Healthy_Fibros.Robj')
load(file = 'FibrofetalNAME_new.Robj')
load(file = 'reindeerfibros.Robj')

write.csv(reindeerfibros@meta.data, file = "reindeerfibros_metadata.csv")
meta_data = read.csv(file = "reindeer_all_days_by_condition.integrated_PC30_metadata.csv") 
new_col = data.frame(meta_data$stage)
reindeer_all_days_by_condition.integrated_PC30$stage = meta_data$stage
save(reindeer_all_days_by_condition.integrated_PC30, file = 'reindeer_all_days_by_condition.integrated_PC30.Robj')

AllFibrosDays <- merge(Healthy_Fibros, y = c(FibrofetalNAME, reindeerfibros))
save(AllFibrosDays, file = 'AllFibrosDays.Robj')
write.csv(AllFibrosDays@meta.data, file = "AllFibrosDays_metadata.csv")


load(file = 'AllFibrosDays.Robj')
# split the dataset into a list of two seurat objects (stim and CTRL)
AllFibroAllDays.list <- SplitObject(AllFibrosDays, split.by = "stage")

# normalize and identify variable features for each dataset independently
AllFibroAllDays.list <- lapply(X = AllFibroAllDays.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = AllFibroAllDays.list)
immune.anchors <- FindIntegrationAnchors(object.list = AllFibroAllDays.list, anchor.features = features)

immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

jpeg(file = "all_fibros_subset_merge_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", group.by = "stage")
dev.off()

jpeg(file = "all_fibros_subset_merge_clusters_dimplot_anno.jpeg", width = 20, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

jpeg(file = "all_fibros_subset_merge_dimplot_anno_split_stage.jpeg", width = 40, height = 25, units = "cm", res = 500)
DimPlot(immune.combined, reduction = "umap", split.by = "stage")
dev.off()

jpeg(file = "all_fibros_subset_crabp1.jpeg", width = 40, height = 25, units = "cm", res = 500)
FeaturePlot(immune.combined, features = "CRABP1", split.by = "stage")
dev.off()

AllFibroAllDays.list = immune.combined

save(AllFibroAllDays.list, file = 'AllFibroAllDays_list.Robj')
write.csv(AllFibroAllDays.list@meta.data, file = "AllFibroAllDays_list_metadata.csv")


load(file = 'reindeer_all_days_by_condition.integrated_PC30.Robj')
reindeerfibros = subset(reindeer_all_days_by_condition.integrated_PC30, idents = "Fibroblast")
save(reindeerfibros, file = 'reindeerfibros.Robj')




