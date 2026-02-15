####1. LOAD LIBRARIES####
library(Matrix)
library(ggplot2)
library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)

####2. Set Visualization options and variables####

options(bitmapType='cairo-png')

#Paths to data
path_Xenium <- "Path_to_Xenium_data_with_computed_ULC_distances"
path_DEGs <- "Path_to_scRNAseq_stroma_DEGS.csv"

#Path to output directory
outs <- "Path_to_Outputs"

####3. Load and process Xenium data####

#Load data
Xenium <- readRDS(file = path_Xenium)

#Filter data according to calculated distances
Xenium <- Xenium[,!is.na(Xenium$Dist_to_ULC)]

#Normalize data
Xenium <- NormalizeData(Xenium)

#Add FOV label metadata
Xenium$FOV <- sub("_[^_]+$", "", Cells(Xenium))

#Add "IF_Fibroblast score based on gene set scRNAseq DEGs"
Markers <- readRDS(path_DEGs)
Markers <- Markers[Markers$cluster == "Inflammatory fibroblasts", ]
Genes <- Markers$gene[Markers$p_val_adj < 0.05 & Markers$avg_log2FC > 1.5]
#Xenium <- Xenium[,Xenium$transcript_counts >= 100] #same results using high quality cells only
Xenium <- AddModuleScore(object = Xenium,
                         slot = "data",
                         features = list(Genes),
                         ctrl = 100,
                         seed = 13,
                         name = "IF_SCORE_scRNAseq")

####4. Plot fibroblast frequency and IF fibroblast score across distance to ULC niche ####

###FIBROBLAST FREQUENCY
#Create dataframe for plot data
df <- data.frame("Distance_to_ULC" = Xenium$Dist_to_ULC,
                 "Fibroblast_type" = Xenium$clusters,
                 "FOV" = Xenium$FOV,
                 "IF_Score" = Xenium$IF_SCORE_scRNAseq1)

#Filter unwanted fibroblast types out
df2 <- df %>%
  filter(Fibroblast_type %in% c(
    "SM PI16hi universal fibroblast",
    "SM CA12hi universal fibroblast",
    "SM GREM1hi fibroblast",
    "SM TNChi fibrotic myofibroblast"))

#Filter only FOVs from affected samples
df2 <- df2[df2$FOV %in% unique(df$FOV)[c(2,3,5,6,8,11,14,16)],]

#Calculate cell density between the min and max distance from ULC niche
#Create distance grid
x_grid <- seq(min(df2$Distance_to_ULC, na.rm = TRUE),
              max(df2$Distance_to_ULC, na.rm = TRUE),
              length.out = 1000)# ~fibroblast density per um
#Helper function to calculate cell density at a certain point of the grid
get_density_on_grid <- function(x, x_grid) {
  if(length(unique(x)) < 2) return(tibble(x = x_grid, y = 0))
  d <- density(x, from = min(x_grid), to = max(x_grid), n = length(x_grid))
  tibble(x = d$x, y = d$y)
}
#Calculate per fibroblast_type/FOV density
dens_by_fov <- df2 %>%
  group_by(Fibroblast_type, FOV) %>%
  summarise(dist_vec = list(Distance_to_ULC), .groups = "drop") %>%
  mutate(dens = map(dist_vec, ~ get_density_on_grid(.x, x_grid))) %>%
  unnest(dens)
#Get density mean and Standard Error(SE) across FOVs
summary_by_type <- dens_by_fov %>%
  group_by(Fibroblast_type, x) %>%
  summarise(
    mean_y = mean(y, na.rm = TRUE),
    sd_y = sd(y, na.rm = TRUE),
    n_fov = n(),
    se_y = sd_y / sqrt(n_fov),
    .groups = "drop"
  )

#Compute Standard Error per Fibroblast type
summary_by_type <- summary_by_type %>%
  mutate(
    se_low = pmax(mean_y - se_y, 0),
    se_high = mean_y + se_y
  )

#Plot densities across distance
fibro_cols <- c(
  "SM CA12hi universal fibroblast" = "darkgreen",
  "SM GREM1hi fibroblast" = "red",
  "SM PI16hi universal fibroblast" = "magenta",
  "SM TNChi fibrotic myofibroblast" =  "purple")

p1 <- ggplot(summary_by_type, 
             aes(x = x, y = mean_y, 
                 color = Fibroblast_type, 
                 fill = Fibroblast_type)) +
  geom_ribbon(aes(ymin = se_low, ymax = se_high), 
              alpha = 0.1, colour = NA) +
  geom_line(size = 1) +
  labs(x = "Distance to ULC (µm)", y = "Fibroblast Density",
       color = "Fibroblast type", fill = "Fibroblast type") +
  ggtitle("Cell Frequency") +
  scale_color_manual(values = fibro_cols) +
  scale_fill_manual(values = fibro_cols) +
  xlim(0, quantile(df$Distance_to_ULC, probs = 0.75)) +
  theme_classic(base_size = 14) +
  theme(title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5)) + NoLegend()

###FIBROBLAST IF GENESET EXPRESSION ACROSS ULC DISTANCE
p2 <- ggplot(df2,
       aes(x = Distance_to_ULC,
           y = IF_Score,
           color = Fibroblast_type,
           fill = Fibroblast_type)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.1) +
  labs(x = "Distance to ULC (µm)", y = "Relative IF Score",
    title = "Inflammatory Fibroblast Score") + 
  xlim(0,quantile(df$Distance_to_ULC,probs = 0.75)) +
  scale_color_manual(values = fibro_cols) +
  scale_fill_manual(values = fibro_cols) + 
  theme_classic(base_size = 14) +
  theme(title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))

#Merge plots
mixed <- p1 | p2

#Save merged plot
ggsave(filename = file.path(outs, "Fibroblast_Density_IFScore_from_ULC.pdf"),
       plot = mixed,
       width = 13, height = 5, units = "in", dpi = "retina")

#### 5. Session details ####

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
#   [1] purrr_1.0.2        tidyr_1.3.0        dplyr_1.1.4        patchwork_1.3.2    Seurat_5.3.1       SeuratObject_5.2.0 sp_2.1-2          
# [8] ggplot2_3.5.2      Matrix_1.6-4 