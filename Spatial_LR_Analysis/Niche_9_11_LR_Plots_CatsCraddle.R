#### 1. Load Libraries ####
library(CatsCradle)
library(Seurat)
library(pheatmap)
library(tidyverse)

####2. Set Visualization options and variables####

options(bitmapType='cairo-png')

#Paths to data
path_Niche11_LR <- "Path_to_Niche11_LR_Affected_FOVs.RDS"
path_Niche9_LR <- "Path_to_Niche9_LR_Affected_FOVs.RDS"
Path_to_LR_Pairs <- "Path_to_LRPairs.csv"

#Path to output directory
outs <- "Path_to_Outputs"

#### 3. Load Data ####

#Load Niche 11 results
ligandReceptorResults_N11 <- readRDS(path_Niche11_LR)

#Load Niche 9 results
ligandReceptorResults_N9 <- readRDS(path_Niche9_LR)

#### 4. Plot L-R Total Number of Interactions per niche ####

#Filter all immune cells and relevant stroma
fibro <- c("SM PI16hi universal fibroblast",
           "SM CA12hi universal fibroblast",
           "SM GREM1hi fibroblast",
           "SM TNChi fibrotic myofibroblast",
           "M COL7A1hi inflammatory fibroblast")
immune <- c("Mast","T_CD8_TE_Cycling","Monocyte","TH_17",
            "T_CD8_TE","Subepithelial macrophage","CLEC9Ahi DC",
            "CD14hi macrophage","B cell-adjacent macrophage",
            "B_plasma_IgA1","T_G/D","T_CD4_EM", "B_memory",
            "Macrophage","TH-17","T_CD4_CM","T_CD8_EM",
            "T_Reg","T_mixed","T_CD4_Cycling","CCL19hi DC",
            "T_CD8_RM","DC_pDC","CXCL12 macrophage",
            "SMC-associated macrophage","Neuron-associated macrophage",
            "T_CD8_CM","CXCL9 macrophage","Phagocytic macrophage",
            "TH_17_Cycling","T_FH","B_germinal center")
cells <- c(fibro, immune)

#Compute total L-R interaction between cells of interest in Niche 9
clusters <- levels(factor(ligandReceptorResults_N9$interactionsOnEdgesMeta$cellTypeA))
clusters <- factor(clusters)
cellTypePerCellTypeLigRecMatrix = makeSummedLRInteractionHeatmap(ligandReceptorResults_N9,
                                                                 clusters, 
                                                                 "total")
dev.off()
cellTypePerCellTypeLigRecMatrix <- cellTypePerCellTypeLigRecMatrix[rownames(cellTypePerCellTypeLigRecMatrix) %in% cells,
                                                                   colnames(cellTypePerCellTypeLigRecMatrix) %in% cells]

#Plot total L-R interaction between cells of interest in Niche 9
pdf(paste0(outs,"/Total_LR_interactions_Fibrotic.pdf"),
    width = 9,   # adjust if needed
    height = 7)
ht1 <- pheatmap(cellTypePerCellTypeLigRecMatrix,
               cluster_rows = T,
               cluster_cols = T,
               cellwidth = 12,
               cellheight = 9,  
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               border_color = NA,
               scale = "none",
               cutree_rows = 4,
               show_rownames = T,
               fontsize_row = 7,
               fontsize_col = 7,
               breaks = seq(0,13, by = 1),
               color = colorRampPalette(c("blue","lightblue","white","orange","red"))(length(seq(0,13, by =1))),
               angle_col = 315,
               main = "Total L-R interactions in Fibrotic Niche")
dev.off()

#Compute total L-R interaction between cells of interest in Niche 11
clusters <- levels(factor(ligandReceptorResults_N11$interactionsOnEdgesMeta$cellTypeA))
clusters <- factor(clusters)
cellTypePerCellTypeLigRecMatrix = makeSummedLRInteractionHeatmap(ligandReceptorResults_N11,
                                                                 clusters, 
                                                                 "total")
dev.off()
cellTypePerCellTypeLigRecMatrix <- cellTypePerCellTypeLigRecMatrix[rownames(cellTypePerCellTypeLigRecMatrix) %in% cells,
                                                                   colnames(cellTypePerCellTypeLigRecMatrix) %in% cells]

#Plot total L-R interaction between cells of interest in Niche 11
pdf(paste0(outs,"/Total_LR_interactions_ULC.pdf"),
    width = 9,
    height = 7)
ht2 <- pheatmap(cellTypePerCellTypeLigRecMatrix,
               cluster_rows = T,
               cluster_cols = T,
               cellwidth = 12,
               cellheight = 9,  
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               border_color = NA,
               scale = "none",
               cutree_rows = 3,
               show_rownames = T,
               fontsize_row = 7,
               fontsize_col = 7,
               breaks = seq(min(cellTypePerCellTypeLigRecMatrix),13, by = 1),
               color = colorRampPalette(c("blue","lightblue","white","orange","red"))(length(seq(min(cellTypePerCellTypeLigRecMatrix),max(cellTypePerCellTypeLigRecMatrix), by =1))),
               angle_col = 315,
               main = "Total L-R interactions in ULC")
dev.off()

#### 5. Plot L-R Significant Interactions per niche ####

#Load L-R Paris to plot
LR_to_filter <- read.csv(Path_to_LR_Pairs)
Pairs <- gsub(pattern = " ",replacement = "_",LR_to_filter$LR.Pair)

#Plot significant L-R interaction between cells of interest in Niche 9

#Extract results
pvals <- ligandReceptorResults_N9$pValues
means <- ligandReceptorResults_N9$meanInteractionsByCluster

#Clean and Pivot P-Values and Mean Interactions
p_long <- pvals %>%
  rownames_to_column("ClusterPair") %>%
  pivot_longer(cols = -ClusterPair, 
               names_to = "LRPair", 
               values_to = "pVal")
m_long <- means %>%
  rownames_to_column("ClusterPair") %>%
  pivot_longer(cols = -ClusterPair, 
               names_to = "LRPair", 
               values_to = "MeanScore")

#Merge data by cluster pairs and by L-R pairs
plot_data <- left_join(p_long, 
                       m_long, 
                       by = c("ClusterPair", "LRPair"))

#Adjust p.values
plot_data$adjpval <- p.adjust(plot_data$pVal,method = "fdr")

#Calculate L-R Score
plot_data$score <- -log10(plot_data$adjpval + as.numeric(.Machine$double.eps)) * plot_data$MeanScore

#Filter significant L-R interactions on cells of interest (immune -> fibro)
top_hits <- plot_data %>% filter(adjpval < 0.05, 
                                 plot_data$MeanScore > 0.01, 
                                 plot_data$ClusterPair %in% 
                                   grep(paste(paste0("-",fibro), collapse="|"), 
                                        plot_data$ClusterPair, value = TRUE),
                                 plot_data$ClusterPair %in% 
                                   grep(paste(paste0(immune,"-"), collapse="|"), 
                                        plot_data$ClusterPair, value = TRUE))

#Filter significant L-R interactions on cells of interest (immune -> fibro)
top_hits <- top_hits[top_hits$LRPair %in% Pairs,]

#Make plot for Fibrotic Niche
colnames(top_hits)[6] <- "Interaction Score"
p1 <- ggplot(top_hits, aes(x = LRPair, y = ClusterPair)) +
  geom_point(aes(size = MeanScore, 
                 fill = `Interaction Score`,
                 alpha = `Interaction Score`),
             shape = 21,color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 12)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "L-R Interactions in Fibrotic Niche",
       size = "Interaction Freq",
       x = "",
       y = "") +
  theme(
    axis.text.x = element_text(angle = 45,size = 10, 
                               hjust = 0.5,
                               vjust = 0.5),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


#Plot significant L-R interaction between cells of interest in Niche 11

#Extract results
pvals <- ligandReceptorResults_N11$pValues
means <- ligandReceptorResults_N11$meanInteractionsByCluster

#Clean and Pivot P-Values and Mean Interactions
p_long <- pvals %>%
  rownames_to_column("ClusterPair") %>%
  pivot_longer(cols = -ClusterPair, 
               names_to = "LRPair", 
               values_to = "pVal")
m_long <- means %>%
  rownames_to_column("ClusterPair") %>%
  pivot_longer(cols = -ClusterPair, 
               names_to = "LRPair", 
               values_to = "MeanScore")

#Merge data by cluster pairs and by L-R pairs
plot_data <- left_join(p_long, 
                       m_long, 
                       by = c("ClusterPair", "LRPair"))

#Adjust p.values
plot_data$adjpval <- p.adjust(plot_data$pVal,method = "fdr")

#Calculate L-R Score
plot_data$score <- -log10(plot_data$adjpval + as.numeric(.Machine$double.eps)) * plot_data$MeanScore

#Filter significant L-R interactions on cells of interest (immune -> fibro)
top_hits <- plot_data %>% filter(adjpval < 0.05, 
                                 plot_data$MeanScore > 0.01, 
                                 plot_data$ClusterPair %in% 
                                   grep(paste(paste0("-",fibro), collapse="|"), 
                                        plot_data$ClusterPair, value = TRUE),
                                 plot_data$ClusterPair %in% 
                                   grep(paste(paste0(immune,"-"), collapse="|"), 
                                        plot_data$ClusterPair, value = TRUE))

#Filter significant L-R interactions on cells of interest (immune -> fibro)
top_hits <- top_hits[top_hits$LRPair %in% Pairs,]

#Make plot for ULC Niche
colnames(top_hits)[6] <- "Interaction Score"
p2 <- ggplot(top_hits, aes(x = LRPair, y = ClusterPair)) +
  geom_point(aes(size = MeanScore, 
                 fill = `Interaction Score`,
                 alpha = `Interaction Score`),
             shape = 21,color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 12)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "L-R Interactions in ULC Niche",
       size = "Interaction Freq",
       x = "",
       y = "") +
  theme(
    axis.text.x = element_text(angle = 45,size = 10, 
                               hjust = 0.5,
                               vjust = 0.5),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#Merge plots
mixed <- p1 | p2

#Save merged plot
ggsave(filename = file.path(outs, "LR_Fibrotic_vs_ULC.pdf"),
       plot = mixed,
       width = 13, height = 7, units = "in", dpi = "retina")

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.4       
# [5] purrr_1.0.2        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1      
# [9] ggplot2_4.0.1      tidyverse_2.0.0    pheatmap_1.0.12    Seurat_5.3.1      
# [13] SeuratObject_5.2.0 sp_2.1-2           CatsCradle_1.1.0  