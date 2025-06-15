loadPackage('sceasy')
loadPackage(SingleCellExperiment)
rds_filename <- 'scATAC_GeneScoreMatrix.rds'
h5ad_filename <- 'scATAC_GeneScoreMatrix.norm.h5ad'
#rds_filename <- 'scATAC_PeakMatrix.rds'
#h5ad_filename <- 'scATAC_PeakMatrix.norm.h5ad'

rds <- readRDS(rds_filename)
sce <- as(assays(rds)[[1]], 'SingleCellExperiment')
names(assays(sce)) <- c('normcounts')
sceasy:::sce2anndata(sce, outFile=h5ad_filename, main_layer='normcounts')
