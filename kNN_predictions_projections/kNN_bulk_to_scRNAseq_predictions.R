####1. Load Required Libraries####
library(Seurat)
library(SCP)
library(ggplot2)
library(Matrix)
library(scatterpie)
library(dplyr)
library(tidyr)

####2. Set Visualization options and variables####

options(bitmapType='cairo-png')

#Paths to data
path_fibro_tx <- "Path_to_invitro_treated_stroma_kallisto_output"
path_scRNAseq <- "Path_to_stroma_scRNAseq"

#Path to output directory
outs <- "Path_to_Outputs"

#Function Variables
Non_studied <- "Combo_none|17uT_none|Combo_Belinostat|Combo_GSKJ4"
freq_filter <- 20

####3. Load and process in vitro bulk RNAseq fibro data####

#Load data
fibro_treat_bulk <- readRDS(path_fibro_tx)

#Filter TPMs
fibro_treat_bulk <- log1p(fibro_treat_bulk$abundance)

#Filter out unnecessary columns
fibro_treat_bulk <- fibro_treat_bulk[, !colnames(fibro_treat_bulk) %in% grep(Non_studied,colnames(fibro_treat_bulk),value = T)]

#Treatments for scRNA-seq predictions:
table(gsub("^[0-9]+_","",colnames(fibro_treat_bulk)))

#Collapse technical replicates
treatments <- unique(gsub("^[0-9]+_", "", colnames(fibro_treat_bulk)))
fibro_treat_merged <- do.call(
  cbind,
  lapply(treatments, function(trt) {
    cols <- which(gsub("^[0-9]+_","",colnames(fibro_treat_bulk)) == trt)
    rowMeans(fibro_treat_bulk[, cols, drop = FALSE])
  })
)
colnames(fibro_treat_merged) <- treatments

#Transform into dgCMatrix
fibro_treat_merged <- Matrix(fibro_treat_merged, sparse = T)

####4. Load and process scRNA-seq stromal data####

#Load scRNA-seq stromal data
scRNAseq <- readRDS(path_scRNAseq)

#Normalize data
scRNAseq <- NormalizeData(scRNAseq)

####5. Apply k-NN based predictions on scRNA-seq stromal data####

#Apply kNN based prediction
scRNAseq_2 <- RunKNNPredict(srt_query = scRNAseq, 
                            bulk_ref = fibro_treat_merged, 
                            filter_lowfreq = freq_filter,
                            nfeatures = length(intersect(rownames(scRNAseq),rownames(fibro_treat_bulk)))
                            )

####6. Plot Prediction Results ####

#Modify function's code for compatibility (SCP)
new_get_legend <- function(plot) {
  plot <- cowplot::as_gtable(plot)
  grob_names <- plot$layout$name
  grobs <- plot$grobs
  grobIndex <- which(grepl("guide-box-bottom", grob_names))
  grobIndex <- grobIndex[1]
  matched_grobs <- grobs[[grobIndex]]
  return(matched_grobs)
}
assignInNamespace("get_legend", new_get_legend, ns = "SCP")

#Modify labels for plotting
colnames(scRNAseq_2@meta.data)[30] <- "Treatment Prediction"
scRNAseq_2$`Treatment Prediction`[scRNAseq_2$KNNPredict_simil < 0.3] <- NA #Remove uncertain cells
scRNAseq_2@meta.data$`Treatment Prediction`[scRNAseq_2@meta.data$`Treatment Prediction` == "unstim_none"] <- "Unstimulated"
scRNAseq_2@meta.data$`Treatment Prediction`[scRNAseq_2@meta.data$`Treatment Prediction` == "Combo_DMSO"] <- "Cytokine mix"


#Plot predictions on top of UMAP
p1 <- CellDimPlot(srt = scRNAseq_2,
                  group.by = "Treatment Prediction",
                  palcolor = c("Unstimulated" = "blue",
                               "Cytokine mix" = "red"),
                  pt.alpha = 0.3,
                  pt.size = 0.6,
                  reduction = "UMAP", 
                  theme_use = "theme_blank") +
  theme_void() +
  theme(
    legend.title = element_text(face =  "bold", size = 13),
    legend.text = element_text(size = 13),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = unit(c(0,0,-10,-10), "pt")
  )

#Prepare cluster centroids for Pie charts plotting
umap <- Embeddings(scRNAseq_2, reduction = "umap") %>% as.data.frame()
colnames(umap) <- c("UMAP_1", "UMAP_2")
df <- cbind(scRNAseq_2@meta.data[c("cell_subtype","Treatment Prediction")], umap)
centroids <- df %>%
  group_by(cell_subtype) %>%
  summarise(UMAP_1 = mean(UMAP_1, na.rm = TRUE),
            UMAP_2 = mean(UMAP_2, na.rm = TRUE),
            n_cells = n()) %>%
  ungroup()

#Compute cell counts for each cell subtype and KNNPredictions
comp <- df %>%
  group_by(cell_subtype, `Treatment Prediction`) %>%
  tally(name = "count") %>%
  pivot_wider(names_from = `Treatment Prediction`, values_from = count, values_fill = 0)

#Join  cell counts composition to cluster centroids
centroid_pies <- centroids %>%
  left_join(comp, by = "cell_subtype")

#Plot predictions on top of UMAP with Pie Charts
p_with_pies <- p1 + geom_scatterpie(data = centroid_pies,
                                    aes(x = UMAP_1, y = UMAP_2, r = 0.7),
                                    cols = c("Cytokine mix", "Unstimulated"),
                                    color = "black",    
                                    linewidth = 0.3,
                                    alpha = 1,
                                    show.legend = FALSE) +
  scale_fill_manual(values = c("Unstimulated" = "blue", "Cytokine mix" = "red"),
                    na.value = "grey70", name = "KNNPredict") +
  theme(
    plot.margin = unit(c(0,0,-15,-10), "pt")
  )
print(p_with_pies)

#Save plot
ggsave(filename = file.path(outs, "kNN_Prediction_with_pies.pdf"),
       plot = p_with_pies,
       width = 6, height = 6, units = "in", dpi = "retina")

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
#   [1] tidyr_1.3.0        dplyr_1.1.4        scatterpie_0.2.1   Matrix_1.6-5       SCP_0.5.6          Seurat_5.3.1      
# [7] SeuratObject_5.2.0 sp_2.1-2           ggplot2_4.0.1  