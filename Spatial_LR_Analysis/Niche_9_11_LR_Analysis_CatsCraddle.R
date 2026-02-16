####1. Load Required Libraries####
library(CatsCradle)
library(Seurat)

#### 2. Load & Process Xenium Data ####

#Load data
path_Xenium_RDS <- "Path_to_Xenium_data_as_RDS"
Xenium <- readRDS(file = path_Xenium_RDS)

#Prepare Seurat object for L-R
Xenium <- UpdateSeuratObject(Xenium)

#Add FOV label metadata
Xenium$FOV <- sub("_[^_]+$", "", Cells(Xenium))

#Filter only samples that come from affected regions
affected <- unique(Xenium$FOV)[c(2,3,5,6,8,11,14,16)]
Xenium <- Xenium[ ,Xenium$FOV %in% affected]

#Update annotations
Xenium$annot_new[Xenium$annot_new == ""] <- Xenium$annot[Xenium$annot_new == ""] 
Xenium$cell_type = factor(Xenium$annot_new)

#### 3. Calculate DelaunayNeighbours per cell per FOV (cut off at 100um distance) ####

#Prepare required objects
clusters = Xenium$cell_type
graphList = list()

#Run neighbor computation
for (sample in names(Xenium@images)){
  print(sample)
  centroids = GetTissueCoordinates(Xenium, image = sample)
  rownames(centroids) = centroids$cell
  print("computing graph")
  spatialGraph = computeNeighboursDelaunay(centroids) # delaunayNeighbours! 
  print("adding distances")  
  spatialGraph = edgeLengthsAndCellTypePairs(spatialGraph, clusters, centroids)
  print("filtering by distance")
  spatialGraph = spatialGraph[spatialGraph$length < 100,] # 100 micron cut-off
  print("extending to 4th degree neighbours")
  extendedNeighboursList = getExtendedNBHDs(spatialGraph[,c(1,2)], 4)
  spatialGraph = collapseExtendedNBHDs(extendedNeighboursList, 4)
  print("storing processed spatial graph")
  graphList[[sample]] = spatialGraph
}

#Combine results across FOVs
delaunayNeighbours_merged = do.call(rbind, graphList)

#### 4. Run L-R analysis in Niche 11 ####

#Filter for cells in niche 11
sample_cells <- colnames(Xenium)[Xenium$spatial_cluster15 == "11"]

#Run analysis
ligandReceptorResults = performLigandReceptorAnalysis(Xenium[, sample_cells], 
                                                      delaunayNeighbours_merged[
                                                        delaunayNeighbours_merged$nodeA %in% sample_cells & 
                                                          delaunayNeighbours_merged$nodeB %in% sample_cells, 1:2], 
                                                      "human", 
                                                      clusters, 
                                                      Xenium$cell_type[sample_cells], 
                                                      method = "analytical", 
                                                      conditional = FALSE,
                                                      minEdgesPos = 10)

#Path to output directory
outs <- "Path_to_Outputs"
saveRDS(ligandReceptorResults, paste0(outs,"/Niche11_LR_Affected_FOVs.RDS"))

#### 5. Run L-R analysis in Niche 9 ####

#Filter for cells in niche 9
sample_cells <- colnames(Xenium)[Xenium$spatial_cluster15 == "9"]

#Run analysis
ligandReceptorResults = performLigandReceptorAnalysis(Xenium[, sample_cells], 
                                                      delaunayNeighbours_merged[
                                                        delaunayNeighbours_merged$nodeA %in% sample_cells & 
                                                          delaunayNeighbours_merged$nodeB %in% sample_cells, 1:2], 
                                                      "human", 
                                                      clusters, 
                                                      Xenium$cell_type[sample_cells], 
                                                      method = "analytical", 
                                                      conditional = FALSE,
                                                      minEdgesPos = 10)
saveRDS(ligandReceptorResults, paste0(outs,"/Niche9_LR_Affected_FOVs.RDS"))

#### 6. Session Information ####

sessionInfo()

# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Rocky Linux 8.10 (Green Obsidian)
# 
# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-0           CatsCradle_1.4.2  