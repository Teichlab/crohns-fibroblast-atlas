####1. Load Required Libraries####
library(Matrix)
library(ggplot2)
library(Seurat)

####2. Set Visualization options and variables####

options(bitmapType='cairo-png')

#Paths to data
path_Xenium <- "Path_to_Xenium_data_as_RDS"

#Path to output directory
outs <- "Path_to_Outputs"

####3. Load and process McGregor scRNA-seq stromal data####

#Load data
Xenium <- readRDS(file = path_Xenium)

#Normalize data
Xenium <- NormalizeData(Xenium)

####4. Make Ligand Expression per Niche Plots####

#Load Ligand & Receptor data
LR_Pairs <- read.csv("Path_to_LRPairs.csv")
Ligands <- unique(sub(" .*", "", LR_Pairs$Ligand.Symbols))
Receptors <- unique(sub(" .*", "", LR_Pairs$Receptor.Symbols))
Receptors <- Receptors[!Receptors %in% c("SDC2", "ITGB6", "ITGB8", "ITGAV", "ACVRL1", "CD48")] #Remove non-canonical

#Make Ligand Plot by niche
dp <- DotPlot(Xenium, 
              Ligands,
              group.by = "spatial_cluster15",
              scale = TRUE,
              col.min = -3,
              col.max = 3) +
  coord_flip() +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  scale_size_continuous(range=c(0.1,8)) + ggtitle("Ligand Expression by Niche") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        title = element_text(size = 12, face = "bold", hjust = 10),
        plot.title =  element_text(hjust = 0.5),
        legend.position = "right") +
  labs(x = "", y = "")

#Save merged plot
ggsave(filename = file.path(outs, "Ligands_by_Niches.pdf"),
       plot = dp,
       width = 8, height = 6, units = "in", dpi = "retina")

####5. Session Information####

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
#   [1] Seurat_5.3.1       SeuratObject_5.2.0 sp_2.1-2           ggplot2_3.5.2      Matrix_1.6-4   