####1. LOAD LIBRARIES####
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

####2. Set Visualization options and variables####

options(bitmapType='cairo-png')

#Paths to data
path_McGregor <- "Path_to_McGregor2025_scRNAseq_RDS_object"
path_TAURUS <- "Path_to_TAURUS_scRNAseq_as_RDS_object"
path_Simon <- "Path_to_stroma_scRNAseq_as_RDS"

#Path to output directory
outs <- "Path_to_Outputs"

####3. Load Ligand and Receptor Genes of Interest####

#Load data
LR_Pairs <- read.csv("Path_to_LRPairs.csv")

#Get Genes
Ligands <- unique(sub(" .*", "", LR_Pairs$Ligand.Symbols))
Receptors <- unique(sub(" .*", "", LR_Pairs$Receptor.Symbols))
Receptors <- Receptors[!Receptors %in% c("SDC2", "ITGB6", "ITGB8", "ITGAV", "ACVRL1", "CD48")] #Remove non-canonical

####4. Load, process and plot L-Rs in McGregor scRNA-seq####

#Load data
Seurat_McGregor <- readRDS(path_McGregor)

#Set Idents based metadata column
Seurat_McGregor$McGregor_2025 <- Idents(Seurat_McGregor)

#Filter out Healthy Samples
Seurat_McGregor <- Seurat_McGregor[,Seurat_McGregor$Type2 != "Healthy"]

#Re-order for plot aesthetics
Seurat_McGregor$McGregor_2025 <- factor(Seurat_McGregor$McGregor_2025, 
                                        levels = c("Enterocytes", "BEST4", "Goblets", "Tuft", "Stem/Undiff", "EECs", "TAs",
                                                   "Myeloid", "Mast", "B-cells", "Plasma", "T-cells",
                                                   "Stromal 1", "Stromal 2", "Stromal 3", "Stromal 4", "Fistula Stroma",
                                                   "Myofibroblasts", "Muscle", "Pericytes", "Endothelium", "Lymphatic", "Glia"))
                                        
#Make Dotplot
p1 <- DotPlot(Seurat_McGregor, 
              Ligands,
              group.by = "McGregor_2025",
              scale = TRUE,
              col.min = -3,
              col.max = 3) +
  coord_flip() +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  scale_size_continuous(range=c(0.1,8)) + ggtitle("Ligands McGregor_2025") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        title = element_text(size = 12, face = "bold", hjust = 10), 
        legend.position = "right") +
  labs(x = "", y = "")

#Save merged plot
ggsave(filename = file.path(outs, "McGregor_2025_Ligands_by_Cell.pdf"),
       plot = p1,
       width = 20, height = 8, units = "in", dpi = "retina")


####5. Load, process and plot L-Rs in TAURUS scRNA-seq####

#Load data
Seurat_TAURUS <- readRDS(path_TAURUS)

#Process data
Seurat_TAURUS <- NormalizeData(Seurat_TAURUS)

#Filter only Crohn's
Seurat_TAURUS <- Seurat_TAURUS[ ,Seurat_TAURUS$Disease == "CD"]
                               
#Make Dotplot
p1 <- DotPlot(Seurat_TAURUS, 
              Ligands,
              group.by = "major",
              scale = TRUE,
              col.min = -3,
              col.max = 3) +
  coord_flip() +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  scale_size_continuous(range=c(0.1,8)) + ggtitle("Ligands TAURUS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        title = element_text(size = 12, face = "bold", hjust = 10), 
        legend.position = "right") +
  labs(x = "", y = "")

#Save merged plot
ggsave(filename = file.path(outs, "TAURUS_Ligands_by_Cell.pdf"),
       plot = mixed,
       width = 20, height = 8, units = "in", dpi = "retina")

#### 6. Session details ####

sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Rocky Linux 8.10 (Green Obsidian)
# 
# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS;  LAPACK version 3.11.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4        patchwork_1.3.2    ggplot2_4.0.1      Seurat_5.3.1       SeuratObject_5.2.0 sp_2.1-2          