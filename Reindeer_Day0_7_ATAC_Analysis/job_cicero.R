setwd("/work/biernaskie_lab/sarthak_sinha/D0_7_scATAC_Seq_Reindeer/signac/")
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(patchwork)
set.seed(1234)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
set.seed(1234)
#BiocManager::install("BSgenome.Btaurus.UCSC.bosTau9")
library(BSgenome.Btaurus.UCSC.bosTau9)
library(hdf5r)
library(ensembldb)
library(AnnotationHub)
library(devtools)
library(monocle3)
library(cicero)
library(SeuratWrappers)


load(file = "reindeer_d0_7_back_velvet.fibros.Robj")
DefaultAssay(reindeer_d0_7_back_velvet.fibros) = "peaks"
reindeer_d0_7_back_velvet.fibros.cds <- as.cell_data_set(x = reindeer_d0_7_back_velvet.fibros)
reindeer_d0_7_back_velvet.fibros.cds.cicero <- make_cicero_cds(reindeer_d0_7_back_velvet.fibros.cds, reduced_coordinates = reducedDims(reindeer_d0_7_back_velvet.fibros.cds)$UMAP)
genome <- seqlengths(reindeer_d0_7_back_velvet.fibros)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns <- run_cicero(reindeer_d0_7_back_velvet.fibros.cds.cicero, genomic_coords = genome.df, sample_num = 100)
head(conns)
save(conns, file = "conns.Robj")
ccans <- generate_ccans(conns)
head(ccans)
save(ccans, file = "ccans.Robj")
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(reindeer_d0_7_back_velvet.fibros) <- links

save(reindeer_d0_7_back_velvet.fibros, file = "reindeer_d0_7_back_velvet.fibros.Robj")

jpeg(file = "CoveragePlot_1_reindeer_d0_7_back_velvet.fibros_NFRKB_2_co-assisible.jpeg", width = 40, height = 15, units = "cm", res = 500)
CoveragePlot(
  object = reindeer_d0_7_back_velvet.fibros,
  region = "29-36670323-36675093",
  extend.upstream = 20,
  extend.downstream = 20)
dev.off()










