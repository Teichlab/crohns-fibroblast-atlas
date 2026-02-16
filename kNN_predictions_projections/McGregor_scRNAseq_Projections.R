####1. Load Required Libraries####
library(Seurat)
library(SCP)
library(ggplot2)

####2. Set Visualization options and variables####

options(bitmapType='cairo-png')

#Paths to data
path_McGregor <- "Path_to_McGregor2025_scRNAseq_RDS_object"
path_Simon <- "Path_to_stroma_scRNAseq_as_RDS"

#Path to output directory
outs <- "Path_to_Outputs"

####3. Load and process McGregor scRNA-seq stromal data####

#Load data
Seurat_McGregor <- readRDS(path_McGregor)

#Filter only mesenchymal stromal cells
Seurat_McGregor$McGregor_2025 <- Idents(Seurat_McGregor)
Seurat_McGregor <- Seurat_McGregor[,Seurat_McGregor$McGregor_2025 %in% c("Fistula Stroma",
                                                                     "Stromal 1","Stromal 2",
                                                                     "Stromal 3","Stromal 4",
                                                                     "Muscle","Pericytes", "Myofibroblasts")]

####4. Load and process Simon scRNA-seq stromal data####

#Load the processed Simon's scRNA-seq fibroblast data
Seurat_Simon <- readRDS(path_Simon)

#Process data
Seurat_Simon <- NormalizeData(Seurat_Simon)

####5. Apply k-NN based predictions on scRNA-seq stromal data####

srt_query <- RunKNNMap(srt_query = Seurat_McGregor,
                       srt_ref = Seurat_Simon,
                       ref_umap = "umap",
                       nfeatures = length(intersect(rownames(Seurat_McGregor),rownames(Seurat_Simon)))
                       )

####6. Make Projection Plots####

#UMAP McGregor Projection onto Simon's dataset
p1 <- ProjectionPlot(
  srt_query = srt_query, 
  srt_ref = Seurat_Simon,
  query_group = "McGregor_2025", 
  ref_group = "cell_subtype",
  query_param = list(palcolor = c("Stromal 3" = "gold",
                                  "Fistula Stroma" = "grey",
                                  "Muscle" = "red",
                                  "Myofibroblasts" = "#FF7F7F",
                                  "Pericytes" =  "brown",
                                  "Stromal 1" = "purple",
                                  "Stromal 2" = "orange", 
                                  "Stromal 4" = "darkgreen"),
                     cells.highlight = TRUE),
  ref_param = list(palcolor = c("Lamina propria fibroblasts" = "purple",
                                "Pericryptal fibroblasts" = "orange", 
                                "SM fibroblasts" = "gold",
                                "Universal fibroblasts" = "skyblue",
                                "FRC-like fibroblasts" = "darkgreen",
                                "Inflammatory fibroblasts" = "grey",
                                "Contractile muscle 1" = "#FF7F7F" ,
                                "Contractile muscle 2" = "red",
                                "Pericyte" =  "brown")))

ggsave(filename = file.path(outs, "Projections_McGregor.pdf"),
       plot = p1,
       width = 12, height = 7, units = "in", dpi = "retina")


####7. Session Information####

sessionInfo()
 
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.10 (Green Obsidian)
# 
# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SCP_0.5.6          Seurat_5.3.1       SeuratObject_5.2.0 sp_2.1-2           ggplot2_4.0.1     
