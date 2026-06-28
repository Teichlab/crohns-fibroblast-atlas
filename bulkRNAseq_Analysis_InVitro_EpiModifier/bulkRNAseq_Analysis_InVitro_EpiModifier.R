library(ggplot2)
library(Rtsne)
library(tidyverse)
library(DESeq2)
library(tximport)
library(rhdf5)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(GOplot)
library(VennDiagram)
library(WGCNA)
library(CorLevelPlot)
library(clusterProfiler)
#library("org.Hs.eg.db")


### importing data from Kallisto quantifications ###

dir<-"/Users/mfriedrich/NDM/Data/proj003//kallisto.dir"
samples<-list.files("/Users/mfriedrich/NDM/Data/proj003/kallisto.dir")
samples<-samples[c(1,5:19)]
samples<-as.data.frame(samples)
colnames(samples)<-"sample"
files <- file.path(dir, samples$sample , "abundance.h5")
names(files)<-samples$sample

#listMarts(host="https://www.ensembl.org")
#ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL")
#listDatasets(ensembl_us_west)
#human_ensembl = useEnsembl("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#listAttributes(human_ensembl)
#tx2genes <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), mart = human_ensembl)
#saveRDS(tx2genes, file="txgenes.rds")
#txi.kallisto <- tximport(files, type = "kallisto", txOut = F, tx2gene=tx2genes)
#saveRDS(txi.kallisto, file="txi.kallisto")

#importing kallisto output (Robject) and metadata 

txi.kallisto<-readRDS(file="txi.kallisto")
countdata<-txi.kallisto[["counts"]]
countdata<-countdata [-1,]
countdata<-countdata[,order(colnames(countdata))]

#metadata 
sample_info<-as.data.frame(samples)
sample_info$combined_label<-sample_info$sample
sample_info<-sample_info %>% separate(sample, into = c("sample_id", "cytokine_treatment",'epimodulator_treatment'), sep="_") %>%
  arrange(sample_id)
row.names(sample_info)<-sample_info$combined_label

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
        panel.background=element_rect(fill="white"))
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
        panel.background=element_rect(fill="white"))
#plot both
merged<-merge(librarySizes, ngenes, by = "sample", all = TRUE)
colnames(merged)[3]<-"ngenes"
merged$stimulation<-sample_info$stimulation[match(merged$sample, sample_info$id)]
merged$combined_label<-merged$sample
merged<-merged %>% separate(sample, into = c("sample_id", "cytokine_treatment",'epimodulator_treatment'), sep="_") %>%
  arrange(sample_id)
merged %>% ggplot (aes(x=librarySizes, y=ngenes, shape=cytokine_treatment,color=epimodulator_treatment)) + 
  geom_point()+
  geom_text_repel(aes(label=combined_label))+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"))

#no poor quality samples

###create DESEQ object###

#make sure rows exactly correspond to columns of the txi.kallisto!
colnames(countdata) %in% sample_info$combined_label
sample_info<-sample_info[(which(colnames(countdata) %in% sample_info$combined_label)),]
sample_info<-as.data.frame(sample_info)
row.names(sample_info)<-sample_info$combined_label

dds <- DESeqDataSetFromTximport(txi.kallisto, sample_info, design = ~epimodulator_treatment)
dds<-dds[-1,]
# filter for at least 4 samples (= all replicates of one treatment group) with a count of 10 or higher
keep <- rowSums(counts(dds) > 10) >= 4 
dds <- dds[keep,]

#normalise using vst (performs better than rlog for bigger datasets)
vsd<-vst(dds, blind=TRUE) #TRUE makes it fully unsupervised, assuming no knowledge on confounding factors
head(assay(vsd),3)

htree <- hclust(dist(t(vsd@assays@data@listData[[1]])), method = "average")
plot(htree)
# all clustering as expected

# Screeplot of PCs
# calculate the variance for each gene
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
#plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:16)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_point(stat="identity")+theme_bw()
# => first 8 components explain most of the variation

#PCA
pcaData<-plotPCA(vsd, intgroup=c("cytokine_treatment", "epimodulator_treatment","sample_id","combined_label"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, shape = cytokine_treatment,color=epimodulator_treatment)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text_repel(aes(label=combined_label),max.overlaps = Inf)+
    ggtitle("PCA with VST data") + theme_classic()
#appears that each treatment group of n=4 has n=1 outlier; possible due to well effect on culture plate
outliers<-c("17_Combo_Belinostat","16_Combo_DMSO","22_Combo_GSKJ4","1_unstim_none")

##exploration on normalised vst - tidyverse
counts<-as.data.frame(vsd@assays@data@listData)
counts$gene<-row.names(counts)
counts.tidy<-gather(counts, sample,vst, -gene)
counts.tidy$sample<-gsub("X","",counts.tidy$sample)
counts.tidy$combined_label<-counts.tidy$sample
counts.tidy<-counts.tidy %>% separate(sample, into = c("sample_id", "cytokine_treatment",'epimodulator_treatment'), sep="_")
counts.tidy$cytokine_treatment<-factor(counts.tidy$cytokine_treatment,
                                       levels=c("unstim","17uT","Combo"))
counts.tidy$epimodulator_treatment<-factor(counts.tidy$epimodulator_treatment,
                                           levels=c("none","DMSO","Belinostat","GSKJ4"))

counts.tidy %>% filter (gene %in% c("CXCL5","CXCL6","SOCS3","ACTA2","DES","MYH11")) %>% 
  ggplot(aes(x=cytokine_treatment, y=vst, fill=epimodulator_treatment))+ 
  geom_boxplot(outlier.color = NA)+
  geom_point(position=position_jitterdodge(),size=1)+
  xlab("")+ylab("raw kallisto counts")+labs(shape="epimodulator_treatment")+
  facet_wrap(.~gene,ncol=1,scales = "free")+
  #stat_compare_means(method= "wilcox.test")+
  #scale_fill_manual(values=c("#FF1500","#F7A41D","#363B8E","#FFFFFF"))+
  #geom_text_repel(aes(label=sample))+
  theme_classic()
#outliers persist, removing now 1 from each treatment group
outliers
dds_no_outliers <- dds[,-which((colnames(dds) %in% outliers)==TRUE)]
dim(dds)
dim(dds_no_outliers)

#re-normalise
vsd<-vst(dds_no_outliers, blind=TRUE) #TRUE makes it fully unsupervised, assuming no knowledge on confounding factors
head(assay(vsd),3)
htree <- hclust(dist(t(vsd@assays@data@listData[[1]])), method = "average")
plot(htree)

#re-plot final PCA
pcaData<-plotPCA(vsd, intgroup=c("cytokine_treatment", "epimodulator_treatment","sample_id","combined_label"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2,color=epimodulator_treatment)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #geom_text_repel(aes(label=combined_label),max.overlaps = Inf)+
  ggtitle("PCA with VST data") + theme_classic()


##differential gene expression; above we defined ~ epimodulator_treatment;
design(dds_no_outliers)
dds_no_outliers$epimodulator_treatment<-factor(dds_no_outliers$epimodulator_treatment, levels = c("none","DMSO","Belinostat","GSKJ4") )
dds_no_outliers$epimodulator_treatment <- relevel(dds_no_outliers$epimodulator_treatment ,ref="DMSO")
deseq<-DESeq(dds_no_outliers)
resultsNames(deseq)

res_Combo<- as.data.frame(results(deseq, name = c("epimodulator_treatment_none_vs_DMSO"))) #note this is inverse: negative FC means up in Combo
res_Bel <- as.data.frame(results(deseq, name = c("epimodulator_treatment_Belinostat_vs_DMSO")))
res_GSK <- as.data.frame(results(deseq, name = c("epimodulator_treatment_GSKJ4_vs_DMSO")))

write.csv(res_Combo, file="DESeq2_[inverse]_unstim_over_DMSO.csv")
write.csv(res_Bel, file="DESeq2_Belinostat_over_DMSO.csv")
write.csv(res_GSK, file="DESeq2_GSKJ4_over_DMSO.csv")

#Exploration by volcano plots
Bel<- res_Bel %>% mutate(gene=row.names(res_Bel)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 1 & padj < 0.05)
Bel %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=Bel[Bel$flag==TRUE,],aes())+
  geom_label_repel(data=Bel %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=10))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
GSK<- res_GSK %>% mutate(gene=row.names(res_GSK)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 1 & padj < 0.05)
GSK %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=GSK[GSK$flag==TRUE,],aes())+
  geom_label_repel(data=GSK %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=20))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
Combo<- res_Combo %>% mutate(gene=row.names(res_Combo)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 1 & padj < 0.05)
Combo %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=Combo[Combo$flag==TRUE,],aes())+
  geom_label_repel(data=Combo %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=20))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()

##plotting comparative heatmaps with all genes regulated by any cytokine##
Bel_reg<-Bel %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame() 
GSK_reg<-GSK %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()

list<-list(Bel_reg,GSK_reg)
all_reg_genes<-Reduce(function(x,y) merge (x,y,all=TRUE), list)
reg_genes<-unique(all_reg_genes$gene) %>% as.data.frame()
colnames(reg_genes)<-"gene"

row.names<-vsd@colData@rownames
norm_matrix<-as.data.frame(vsd@assays@data@listData)
colnames(norm_matrix)<-row.names
colnames(norm_matrix)
row.names(sample_info) #check if order is matching

reg_genes_matrix<-merge(norm_matrix,reg_genes, by.x=0, by.y=1)
row.names(reg_genes_matrix)<-reg_genes_matrix[,1]
reg_genes_matrix<-as.matrix(reg_genes_matrix[,-1])
reg_genes_matrix<-reg_genes_matrix[ , order(colnames(reg_genes_matrix))]
pheatmap(reg_genes_matrix,
         scale="row",
         legend_breaks = seq(0,14, by=1),
         border_color=NA,
         #cutree_rows = 5,
         #gaps_col = seq(4,40,by=4),
         cellwidth = 2.5,
         cellheight = 1.3,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         fontsize_row = 1,
         fontsize_col = 2.5,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         )

