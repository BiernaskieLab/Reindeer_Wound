library(ggplot2)
library(garnett)
library(Seurat)
library(dittoSeq)

#Convert Training Seurat object into CellDataSet object 
load("back_antler_fibro.Robj")
R3D0.integrated_PC10 <- back_antler_fibro
#add umap coordinates to metadata
UMAP = as.data.frame(R3D0.integrated_PC10[['umap']]@cell.embeddings)
R3D0.integrated_PC10@meta.data$UMAP_1 = UMAP$UMAP_1
R3D0.integrated_PC10@meta.data$UMAP_2 = UMAP$UMAP_2

Idents(R3D0.integrated_PC10) <- "fibro_state"
R3D0.integrated_PC10@meta.data$celltype <- Idents(R3D0.integrated_PC10)

data <- as(as.matrix(R3D0.integrated_PC10@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = R3D0.integrated_PC10@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
d0_fibro_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

# generate size factors for normalization later
d0_fibro_cds <- estimateSizeFactors(d0_fibro_cds)
saveRDS(d0_fibro_cds, file = "d0_fibro_cds.RDS")

#Check Markers 
marker_check <- check_markers(d0_cds, "fibro_markers.txt",
                              db='none')

svg(filename = "test_markers.svg")
print(plot_markers(marker_check)) #32 genes
dev.off()

#Verifying the classifier 
d0_fibro_classifier <- train_cell_classifier(cds = d0_fibro_cds,
                                         marker_file = "fibro_markers.txt",
                                         db='none',
                                         num_unknown = 50)

feature_genes <- get_feature_genes(d0_fibro_classifier,
                                   node = "root",
                                   db = 'none')
head(feature_genes)

#Testing the classifier  
d0_fibro_cds <- classify_cells(d0_fibro_cds, d0_fibro_classifier,
                           db = "none",
                           cluster_extend = TRUE,
                          cluster_extend_max_frac_incorrect	= 0.05)


svg(filename = "d0_fibro_garnettPrediction_nonextended.svg")
print(ggplot(data = pData(d0_fibro_cds), aes(x = UMAP_1, y = UMAP_2, colour = cell_type)) + geom_point())
dev.off()

svg(filename = "d0_fibro_garnettPrediction_extended.svg")
print(ggplot(data = pData(d0_fibro_cds), aes(x = UMAP_1, y = UMAP_2, colour = cluster_ext_type)) + geom_point())
dev.off()

#Using classifier on human fibroblasts 
human_fibro_cds <- classify_cells(human_fibro_cds, d0_fibro_classifier,
                            db = "none",
                            cluster_extend = TRUE,
                            cluster_extend_max_frac_incorrect	= 0.05)

svg(filename = "human_fibro_prediction_extended.svg")
print(ggplot(data = pData(human_fibro_cds), aes(x = UMAP_1, y = UMAP_2, colour = cluster_ext_type)) + geom_point())
dev.off()

svg(filename = "human_fibro_prediction_nonextended.svg")
print(ggplot(data = pData(human_fibro_cds), aes(x = UMAP_1, y = UMAP_2, colour = cell_type)) + geom_point())
dev.off()

#Convert back to Seurat objects and make percentage classified bar charts 
#Reindeer fibroblasts
count.mat <- Biobase::exprs(d0_fibro_cds)
meta.df <- Biobase::pData(d0_fibro_cds)
d0_reindeer_seurat <- CreateSeuratObject(counts = count.mat,
                                      project = "my.project",
                                      assay = "RNA",
                                      meta.data = meta.df)
save(d0_reindeer_seurat, file = "d0_reindeer_seurat.Robj")

p1 <- dittoBarPlot(d0_fibro_seurat, "cell_type", group.by = "fibro_state")
p2 <- dittoBarPlot(d0_fibro_seurat, "cell_type", group.by = "sample_ident")
p3 <- dittoBarPlot(d0_fibro_seurat, "cluster_ext_type", group.by = "fibro_state")
p4 <- dittoBarPlot(d0_fibro_seurat, "cluster_ext_type", group.by = "sample_ident")

svg(filename = "reindeer_percentClassified_nonextended.svg")
print(p1 + p2)
dev.off()

svg(filename = "reindeer_percentClassified_extended.svg")
print(p3 + p4)
dev.off()

#Human fiborblasts
count.mat <- Biobase::exprs(human_fibro_cds)
meta.df <- Biobase::pData(human_fibro_cds)
d0_human_seurat <- CreateSeuratObject(counts = count.mat,
                                         project = "my.project",
                                         assay = "RNA",
                                         meta.data = meta.df)
save(d0_human_seurat, file = "d0_human_seurat.Robj")

svg(filename = "human_percentClassified_nonextended.svg")
print(dittoBarPlot(d0_fibro_seurat, "cell_type", group.by = "stage"))
dev.off()

svg(filename = "human_percentClassified_extended.svg")
print(dittoBarPlot(d0_fibro_seurat, "cluster_ext_type", group.by = "stage"))
dev.off()
