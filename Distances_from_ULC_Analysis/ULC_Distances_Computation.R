#### 1. Load Libraries ####
library(Matrix)
library(ggplot2)
library(Seurat)
library(phenoptr)
library(dplyr)

#### 2. Load & Process Xenium Data ####

#Load data
path_Xenium_RDS <- "Path_to_Xenium_data_as_RDS"
Xenium <- readRDS(file = path_Xenium_RDS)

#Set cells of interest for distance analysis
fibro <- c("SM PI16hi universal fibroblast", "SM CA12hi universal fibroblast", 
           "SM GREM1hi fibroblast", "MP DPThi fibroblast","M COL7A1hi inflammatory fibroblast",
           "S SOCS3hi fibroblast", "M ADAMDEChi fibroblast", 
           "CMP SYNMhi muscle", "CMP GREM1hi muscle", "LMP NR4A1hi muscle", 
           "SM&M PDGFRBhi pericyte", "M COL4A1hi pericryptal fibroblast", 
           "SM CCL19hi FRC-like fibroblast", "SM&M CCL19hi pericyte", 
           "SM GREM1hi fibrotic muscle", "MM NR4A1hi muscle", 
           "SM MCAMhi perivascular muscle","SM TNChi fibrotic myofibroblast")

#Create annotations for distance analysis
Xenium$clusters <- Xenium$annot
Xenium$clusters[Xenium$spatial_cluster15 == "11"] <- "ULC"

### 3. Calculate distances to closest Niche_11/ULC cell ####

#Start metadata column to store results
Xenium$Dist_to_ULC <- NA

for (i in 1:length(names(Xenium@images))) {
  
  #Get FOV data 
  df <- GetTissueCoordinates(Xenium[[ names(Xenium@images)[i]]], which = "centroids")
  cells  <- df$cell #extract cell IDs
  rownames(df) <- cells #name rows with cell IDs
  df <- df[,-3] #remove cell IDs column
  
  #Add Cell/Niche Annotations
  df$clusters <- Xenium$clusters[match(rownames(df), Cells(Xenium))]
  
  #Prepare for distance calculation
  cds <- df %>%
    mutate(
      `Cell X Position` = x,
      `Cell Y Position` = y,
      Phenotype = as.character(clusters),
      `Cell ID` = row_number()
    ) %>%
    select(`Cell X Position`, `Cell Y Position`, Phenotype, `Cell ID`) %>%
    as_tibble()
  rownames(cds) <- cells #name rows with cell IDs
  
  #Filter relevant cells
  Filtered.cells <- cds %>% 
    filter(Phenotype %in% c(fibro,"ULC"))
  
  #If no relevant cells, skip FOV
  if (nrow(Filtered.cells) == 0) {
    message("FOV ", i, "/", length(names(Xenium@images)), " (", names(Xenium@images)[i], "): no fibro/ULC cells -> skipping")
    next
  }
  
  #Compute per-cell nearest neighbor distances
  dist <- find_nearest_distance(Filtered.cells)
  gc()
  
  #Join data
  Cell_with_distance <- bind_cols(Filtered.cells, dist)
  
  #Add distances
  Xenium$Dist_to_ULC[match(rownames(Cell_with_distance), Cells(Xenium))] <- Cell_with_distance$`Distance to ULC`
  
  #End loop
  if (i == length(names(Xenium@images))) {
    
    saveRDS(Xenium, file = paste0(dirname(path_Xenium_RDS),
                                  "/All_Xenium_with_dist_RDSobj.RDS"))
    message("All done! Object with distances saved in: ", dirname(path_Xenium_RDS))
    
  } else {
    
    print(i)
    
  }

}

#### 4. Session Information ####

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
#   [1] dplyr_1.1.4        phenoptr_0.3.2     Seurat_5.3.1       SeuratObject_5.2.0 sp_2.1-2           ggplot2_3.5.2     
# [7] Matrix_1.6-4      