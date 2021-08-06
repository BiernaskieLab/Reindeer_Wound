setwd("/home/ajaffer/SCENIC_Runs/SCENIC_Reindeer_all_fibros")
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)
scenicOptions <- readRDS("int/scenicOptions.Rds")

load("x.Robj")
exprMat <- x@assays$RNA@counts
exprMat <- as.matrix(exprMat)
genesKept <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- exprMat[genesKept,]
exprMat_filtered <- log2(exprMat_filtered+1) 

runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
#runSCENIC_4_aucell_binarize(scenicOptions)

