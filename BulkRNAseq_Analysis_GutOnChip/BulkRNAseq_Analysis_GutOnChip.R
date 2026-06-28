library(tximport)
library(rhdf5)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(rstatix)

`%nin%` = Negate(`%in%`) 


txi.kallisto<-readRDS(file="txi.kallisto.proj010.vascular.samples.rds")
countdata<-txi.kallisto[["counts"]]
countdata<-countdata [-1,]
dim(countdata) #12 total samples

#metadata proj010
meta<-read.csv("proj010_metadata.csv", header=T, na.string=c(""," ","NA","NA "," NA","#DIV/0!","NA  ","none"),stringsAsFactors = T)
(colnames(countdata) %in% meta$sample_name) #all sequenced samples have metadata
meta[duplicated(meta[,]) | duplicated(meta[,], fromLast=TRUE) ,] #fromLast to show all rows, not just duplicated ones

#check library sizes 
librarySizes<-colSums(countdata)
librarySizes<-as.data.frame(librarySizes)
row.names(librarySizes)->librarySizes$sample
librarySizes %>% ggplot (aes(x=sample, y=librarySizes)) + 
  geom_bar(stat="identity",colour="grey",fill="black")+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#check # of genes per sample
iterations = ncol(countdata)
output<-matrix(ncol=iterations,nrow=1)
for(i in 1:iterations){
  output[,i]<-length(which(countdata[,i]!=0))
}
output<-data.frame(output)
ngenes<-data.frame(t(output))
ngenes$sample<-colnames(countdata)
ngenes %>% ggplot (aes(x=sample, y=t.output.)) + 
  geom_bar(stat="identity",colour="grey",fill="black")+
  xlab("sample")+
  ylab("detected genes")+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#plot both
merged<-merge(librarySizes, ngenes, by = "sample", all = TRUE)
colnames(merged)[3]<-"ngenes"
merged %>% ggplot (aes(x=librarySizes, y=ngenes)) + 
  geom_point()+
  geom_text_repel(aes(label=sample),max.overlaps=Inf)+
  #scale_color_manual(values=mycolors)+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"))
#no strong outliers by #genes (all > 18000) or #reads

# detect outlier samples by hierarchical clustering
htree <- hclust(dist(t(countdata)), method = "average")
plot(htree)

#create DESEQ object
#match counts and metadata
colnames(txi.kallisto$counts)
sample_info<-meta
sample_info$sample_name
reorder_idx<-match(colnames(txi.kallisto$counts),sample_info$sample_name) #match returns positions to match onto first vector
sample_info_reordered<-sample_info[reorder_idx,] #using the indices to reorder second vector
sample_info<-sample_info_reordered
dds <- DESeqDataSetFromTximport(txi.kallisto, sample_info, design = ~1) 
dds<-dds[-1,]
metadata<-as.data.frame(colData(dds))

# filter for at least 2 samples with a count of 15 or higher
dim(dds)
keep <- rowSums(counts(dds) > 15) >= 2
dds <- dds[keep,]
dim(dds)

#removing qc outliers before normalisation
#dds <- dds[,!colnames(dds) %in% raw_qc_outliers]
#dim(dds)

#normalise using vst 
vsd<-vst(dds, blind=TRUE) 
head(assay(vsd),3)

#check distributions
as.data.frame(assay(vsd)) %>% ggplot(aes(x=rowSums(assay(vsd)))) + 
  geom_density() + 
  scale_x_continuous(trans='log2')

# Screeplot of PCs
## calculate the variance for each gene
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
##plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:24)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_point(stat="identity")+theme_bw()
# => first 2-3 components cover most variance

#PCA(using DESEQ package this time) from scratch with ggplot
pcaData<-plotPCA(vsd, pcsToUse=1:2 , intgroup=c("sample_id","sample_name","treatment_1","treatment_2", "treatment_3", "treatment_combined","replicate","chip_stim_side"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData[["data"]], aes(x = PC1, y = PC2, color = chip_stim_side)) +
  geom_point(size =1) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  #geom_text_repel(aes(label=name),max.overlaps = Inf)+
  ggtitle("PCA with VST data") + theme_bw()

# as expected, PC1 separates vascular from epithelial stimulation of chip
# PC2 possibly separating unstim vs stim
# some separation by replicate on PC2 on vascular stimulations: may need to remove R2 (vascular)

#Testing for differential genes that change across +/- Belinostat
dds <- DESeqDataSetFromTximport(txi.kallisto, sample_info, design = ~treatment_combined) 
dds<-dds[-1,]
dds_subset<-dds[,c(7,8,15,16,23,24)] #subsetting to vascular only, treatment 3 and 4 only
metadata<-as.data.frame(colData(dds_subset))

dim(dds_subset)
keep <- rowSums(counts(dds_subset) > 15) >= 2
dds_subset <- dds_subset[keep,]
dim(dds_subset)
dds_subset$treatment_combined<-factor(dds_subset$treatment_combined, levels = c("stim_fibroblast_SN","stim_BEL_fibroblast_SN"))
deseq<-DESeq(dds_subset) 
results<-as.data.frame(results(deseq)) %>% arrange(padj)

# first 503 genes are p < 0.05)
counts<-as.data.frame(vsd@assays@data@listData)
colnames(counts) <- gsub(".", "-", colnames(counts), fixed = TRUE)
counts$gene<-row.names(counts)
counts.tidy<-gather(counts, sample_name,vst, -gene)
counts.tidy<-merge(counts.tidy, sample_info_reordered, by="sample_name")
counts.tidy$treatment_combined<-factor(counts.tidy$treatment_combined, 
                                       levels = c("control","unstim_fibroblast_SN","stim_fibroblast_SN","stim_BEL_fibroblast_SN"))
selection<-counts.tidy %>% filter(chip_stim_side=="vascular") %>% pull (sample_name)
selection<-which(colnames(counts)%in%selection)
heatmap <- counts[,selection] #subsetting to vascular only
heatmap<-heatmap[c(1,5,9,2,6,10,3,7,11,4,8,12)]
View(sample_info_reordered)
row.names(sample_info_reordered)<-sample_info_reordered$sample_name
goi<-rownames(results)[1:503]
heatmap <- heatmap [goi,]
pheatmap (heatmap,scale='row',
         clustering_method = "average",
         border_color = NA,
         cellwidth = 10,
         cellheight = 10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #breaks = seq(-4, 4, 0.05),
         #color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = sample_info_reordered[,c(5,6,7)]
)

interesting_hits <- c("RELB","SOX6","TFF3","CDX2","CA12","COL6A1","MUC3A","CXCL3","CXCL2","KRT20","VCAM1","MUC17","S100A8","S100A9")

##validation on normalised vst - tidyverse

#removing replicate (=Chip1) of control group, as clear outlier and noticed leakiness of chip at experiment#
counts.tidy %>% filter (gene %in% interesting_hits) %>%
  filter (sample_name != "Replicate-1---vascular---Treatment-1") %>%
  filter (chip_stim_side == "vascular") %>%
  ggplot(aes(x=treatment_combined, y=vst,group=replicate))+ 
  geom_violin(aes(group=treatment_combined,color=treatment_combined))+
  geom_jitter(aes(group=treatment_combined),width=0.1,size=1)+
  #geom_path(linewidth=0.1)+
  stat_summary(aes(group=treatment_combined),fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1)+
  xlab("")+ylab("Raw counts (kallisto)")+labs(fill="stimulation")+
  facet_wrap(.~gene, nrow=2, scales="free")+
  #scale_y_continuous(trans='log10')+
  #stat_compare_means(method= "wilcox.test")+
  #scale_fill_manual(values=c("#FF1500","#F7A41D","#363B8E","#FFFFFF"))+
  #geom_text_repel(aes(label=sample_name),max.overlaps = Inf)+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4,colour="black"))

stat.test<-counts.tidy %>% 
  filter (gene %in% interesting_hits) %>% 
  filter (chip_stim_side == "vascular") %>%
  group_by(gene) %>%
  pairwise_t_test(vst~treatment_combined, ref.group = "stim_fibroblast_SN", p.adjust.method = "fdr")
View(stat.test)

#significant: CDX2, COL6A1, CXCL2, MUC17, MUC3A, RELB, S100A9, SOX6, TFF3
