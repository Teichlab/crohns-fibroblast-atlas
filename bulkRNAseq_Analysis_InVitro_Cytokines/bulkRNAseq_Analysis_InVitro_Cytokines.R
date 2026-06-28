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
#library("org.Hs.eg.db")

### importing data from Kallisto quantifications ###

dir<-"/Users/mfriedrich/NDM/Data/proj002/kallisto_individual_cytokines/kallisto.dir"
samples<-list.files("/Users/mfriedrich/NDM/Data/proj002/kallisto_individual_cytokines/kallisto.dir")
samples<-as.data.frame(samples)
files <- file.path(dir, samples$samples, "abundance.h5")
names(files)<-samples$samples

#listMarts(host="https://www.ensembl.org")
#ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL")
#listDatasets(ensembl_us_west)
#human_ensembl = useEnsembl("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#listAttributes(human_ensembl)
#tx2genes <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), mart = human_ensembl)
#saveRDS(tx2genes, file="txgenes.rds")
#txi.kallisto <- tximport(files, type = "kallisto", txOut = F, tx2gene=tx2genes)
#saveRDS(txi.kallisto, file="txi.kallisto_individual_cytokines")


#importing kallisto output (Robject) and metadata 
txi.kallisto<-readRDS(file="txi.kallisto_individual_cytokines")
countdata<-txi.kallisto[["counts"]]
countdata<-countdata [-1,]

#metadata 
sample_info<-read.csv("metadata_individual_cytokines.csv", stringsAsFactors = TRUE)
names(sample_info)<-c("id","stimulation")

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
nb.cols<-length(levels(factor((merged$stimulation))))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
merged %>% ggplot (aes(x=librarySizes, y=ngenes, color=stimulation)) + 
  geom_point()+
  geom_text_repel(aes(label=sample))+
  scale_color_manual(values=mycolors)+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"))

#38,27,9,21 are slight outliers, but still well within range of expected #of genes detected; 

###create DESEQ object###

#make sure rows exactly correspond to columns of the txi.kallisto!
colnames(countdata) %in% sample_info$id
sample_info<-sample_info[(which(sample_info$id %in% colnames(countdata))),]
sample_info<-as.data.frame(sample_info)
row.names(sample_info)<-sample_info$id
sample_info<-sample_info[sort(colnames(countdata)),] 

dds <- DESeqDataSetFromTximport(txi.kallisto, sample_info, design = ~stimulation)
dds<-dds[-1,]
# filter for at least 4 samples (= all replicates of one treatment group) with a count of 10 or higher
keep <- rowSums(counts(dds) > 10) >= 4 
dds <- dds[keep,]

#normalise using vst (performs better than rlog for bigger datasets)
vsd<-vst(dds, blind=TRUE) #TRUE makes it fully unsupervised, assuming no knowledge on confounding factors
head(assay(vsd),3)

htree <- hclust(dist(t(vsd@assays@data@listData[[1]])), method = "average")
plot(htree)
# not striking outliers

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
scree_plot[,2]<- c(1:44)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_point(stat="identity")+theme_bw()
# => first two components explain most of the variation

#PCA
pcaData<-plotPCA(vsd, intgroup=c("id", "stimulation"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = stimulation)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text_repel(aes(label=id),max.overlaps = Inf)+
    ggtitle("PCA with VST data") + theme_classic()
#all well within expected

##exploration on normalised vst - tidyverse
counts<-as.data.frame(vsd@assays@data@listData)
counts$gene<-row.names(counts)
counts.tidy<-gather(counts, sample,vst, -gene)
counts.tidy$sample<-gsub("X","",counts.tidy$sample)
counts.tidy$stimulation[counts.tidy$sample %in% c(1:4)] <- "unstimulated"  
counts.tidy$stimulation[counts.tidy$sample %in% c(9:12)] <- "IL1B"  
counts.tidy$stimulation[counts.tidy$sample %in% c(13:16)] <- "IL6"
counts.tidy$stimulation[counts.tidy$sample %in% c(17:20)] <-"IL13"
counts.tidy$stimulation[counts.tidy$sample %in% c(21:24)] <-"IL17A"
counts.tidy$stimulation[counts.tidy$sample %in% c(25:28)] <-  "IL22"
counts.tidy$stimulation[counts.tidy$sample %in% c(29:32)] <-  "IFNA"
counts.tidy$stimulation[counts.tidy$sample %in% c(33:36)] <-  "IFNG"
counts.tidy$stimulation[counts.tidy$sample %in% c(37:40)] <-  "TGFB1"
counts.tidy$stimulation[counts.tidy$sample %in% c(41:44)] <-  "TNF"
counts.tidy$stimulation[counts.tidy$sample %in% c(45:48)] <-  "OSM"

counts.tidy$stimulation<-factor(counts.tidy$stimulation, levels=c("unstimulated","IL1B","TNF","IL17A","OSM","IL6","IL13","IL22","IFNA","IFNG","TGFB1"))

counts.tidy %>% filter (gene %in% c("PDPN","THY1","COL7A1","MMP1","MMP3","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","IL11")) %>% 
  ggplot(aes(x=stimulation, y=vst, fill=stimulation))+ 
  geom_boxplot()+
  geom_jitter(width=0.1,size=1)+
  xlab("")+ylab("vst-normalised expression")+labs(fill="stimulation")+
  facet_wrap(.~gene,ncol=1,scales="free")+
  #scale_y_continuous(breaks=seq(0,15,by=1))+
  #stat_compare_means(method= "wilcox.test")+
  #scale_fill_manual(values=c("#FF1500","#F7A41D","#363B8E","#FFFFFF"))+
  #geom_text_repel(aes(label=sample))+
  theme_classic()
#no obvious outliers by checking of replicates, so continuing with all samples

#re-plot final PCA (FIGURE)
pcaData<-plotPCA(vsd, intgroup=c("id", "stimulation"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData %>% 
  ggplot(aes(x = PC1, y = PC2, color = stimulation)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  #geom_text_repel(aes(label=id),max.overlaps = Inf)+
  scale_colour_manual(values=c("purple3", "yellow2","red","green","tan","turquoise",'brown',"magenta","blue",'cyan',"black"))+
  ggtitle("PCA with VST data") + theme_classic()


##differential gene expression; above we defined ~ stimulation; so by default it will compare stimulation by highest level against lowest
#setting unstimulated as reference level:
dds$stimulation <- relevel(dds$stimulation, ref = "unstimulated")
deseq<-DESeq(dds)
resultsNames(deseq)
res_IFNA <- as.data.frame(results(deseq, name = c("stimulation_IFNa_vs_unstimulated")))
res_IFNG <- as.data.frame(results(deseq, name = c("stimulation_IFNg_vs_unstimulated")))
res_IL13 <- as.data.frame(results(deseq, name = c("stimulation_IL13_vs_unstimulated")))
res_IL17A <- as.data.frame(results(deseq, name = c("stimulation_IL17A_vs_unstimulated")))
res_IL22 <- as.data.frame(results(deseq, name = c("stimulation_IL22_vs_unstimulated")))
res_IL6 <- as.data.frame(results(deseq, name = c("stimulation_IL6_vs_unstimulated")))
res_IL1B <- as.data.frame(results(deseq, name = c("stimulation_IL1B_vs_unstimulated")))
res_OSM <- as.data.frame(results(deseq, name = c("stimulation_OSM_vs_unstimulated" )))
res_TNF <- as.data.frame(results(deseq, name = c("stimulation_TNFa_vs_unstimulated")))
res_TGFB1 <- as.data.frame(results(deseq, name = c("stimulation_TGFb1_vs_unstimulated")))

res_IFNA$gene<-row.names(res_IFNA) 
res_IFNG$gene<-row.names(res_IFNG)
res_IL13$gene <- row.names(res_IL13)
res_IL17A$gene <- row.names(res_IL17A)
res_IL22$gene <- row.names(res_IL22)
res_IL1B$gene <- row.names(res_IL1B)
res_IL6$gene <- row.names(res_IL6)
res_OSM$gene <- row.names(res_OSM)
res_TNF$gene <- row.names(res_TNF)
res_TGFB1$gene <- row.names(res_TGFB1)

write.csv(res_IFNA, file="DESeq2_IFNA.csv")
write.csv(res_IFNG, file="DESeq2_IFNG.csv")
write.csv(res_IL13, file="DESeq2_IL13.csv")
write.csv(res_IL17A, file="DESeq2_IL17A.csv")
write.csv(res_IL22, file="DESeq2_IL22.csv")
write.csv(res_IL1B, file="DESeq2_IL1B.csv")
write.csv(res_IL6, file="DESeq2_IL6.csv")
write.csv(res_OSM, file="DESeq2_OSM.csv")
write.csv(res_TNF, file="DESeq2_TNF.csv")
write.csv(res_TGFB1, file="DESeq2_TGFB1.csv")

#Exploration by volcano plots
IL1B<- res_IL1B %>% mutate(gene=row.names(res_IL1B)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 4 & padj < 0.05)
IL1B %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IL1B[IL1B$flag==TRUE,],aes())+
  geom_label_repel(data=IL1B %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,120,by=10))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
TNF<- res_TNF %>% mutate(gene=row.names(res_TNF)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 4 & padj < 0.05)
TNF %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=TNF[TNF$flag==TRUE,],aes())+
  geom_label_repel(data=TNF %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=20))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
OSM<- res_OSM %>% mutate(gene=row.names(res_OSM)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 2 & padj < 0.05)
OSM %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=OSM[OSM$flag==TRUE,],aes())+
  geom_label_repel(data=OSM %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=5))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
IL17A<- res_IL17A %>% mutate(gene=row.names(res_IL17A)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 1 & padj < 0.05)
IL17A %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IL17A[IL17A$flag==TRUE,],aes())+
  geom_label_repel(data=IL17A %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=10))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
IL22<- res_IL22 %>% mutate(gene=row.names(res_IL22)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 1 & padj < 0.05)
IL22 %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IL22[IL22$flag==TRUE,],aes())+
  geom_label_repel(data=IL22 %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=1))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
IL6<- res_IL6 %>% mutate(gene=row.names(res_IL6)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 1 & padj < 0.05)
IL6 %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IL6[IL6$flag==TRUE,],aes())+
  geom_label_repel(data=IL6 %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=1))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
IL13<- res_IL13 %>% mutate(gene=row.names(res_IL13)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 2 & padj < 0.05)
IL13 %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IL13[IL13$flag==TRUE,],aes())+
  geom_label_repel(data=IL13 %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-20,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=10))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
TGFB1<- res_TGFB1 %>% mutate(gene=row.names(res_TGFB1)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 2 & padj < 0.05)
TGFB1 %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=TGFB1[TGFB1$flag==TRUE,],aes())+
  geom_label_repel(data=TGFB1 %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-20,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=5))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
IFNA<- res_IFNA %>% mutate(gene=row.names(res_IFNA)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 4 & padj < 0.05)
IFNA %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IFNA[IFNA$flag==TRUE,],aes())+
  geom_label_repel(data=IFNA %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-20,10,by=1))+
  scale_y_continuous(breaks=seq(0,100,by=5))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()
IFNG<- res_IFNG %>% mutate(gene=row.names(res_IFNG)) %>% arrange (padj,abs(log2FoldChange)) %>% mutate(flag = abs(log2FoldChange) >= 4 & padj < 0.05)
IFNG %>% ggplot(aes(x=log2FoldChange, y=-log10(padj),color=flag)) +
  geom_point(alpha=0.3)+
  geom_point(data=IFNG[IFNG$flag==TRUE,],aes())+
  geom_label_repel(data=IFNG %>% filter(flag), aes(label=gene),fill="white",alpha=0.7, size=2, max.overlaps = Inf,fontface="bold")+
  geom_hline(yintercept=-log10(0.05), col="blue",linetype="dashed")+
  geom_vline(xintercept=1, col="blue",linetype="dashed")+
  geom_vline(xintercept=-1, col="blue",linetype="dashed")+
  scale_x_continuous(breaks=seq(-20,10,by=1))+
  scale_y_continuous(breaks=seq(0,200,by=20))+
  scale_colour_manual(values=c("red","black"))+
  theme_classic()

##look at number of significantly (Padjust<0.05, abs(FC)>2) regulated genes per cytokine

list<-list(IL1B=IL1B,TNF=TNF,IL6=IL6,IL13=IL13,IL17A=IL17A,IL22=IL22,IFNA=IFNA,IFNG=IFNG,OSM=OSM,TGFB1=TGFB1)
iterations = length(list)
output_up<-matrix(ncol=iterations,nrow=1)
for(i in 1:iterations){
  output_up[,i]<-list[[i]] %>% filter (padj<0.05 & log2FoldChange>1) %>% pull(gene) %>% length()
}
output_up<-data.frame(output_up)
colnames(output_up)<-c("IL1B","TNF","IL6","IL13","IL17A","IL22","IFNA","IFNG","OSM","TGFB1")

list<-list(IL1B=IL1B,TNF=TNF,IL6=IL6,IL13=IL13,IL17A=IL17A,IL22=IL22,IFNA=IFNA,IFNG=IFNG,OSM=OSM,TGFB1=TGFB1)
iterations = length(list)
output_down<-matrix(ncol=iterations,nrow=1)
for(i in 1:iterations){
  output_down[,i]<-list[[i]] %>% filter (padj<0.05 & log2FoldChange<(-1)) %>% pull(gene) %>% length()
}
output_down<-data.frame(output_down)
colnames(output_down)<-c("IL1B","TNF","IL6","IL13","IL17A","IL22","IFNA","IFNG","OSM","TGFB1")

ngenes_up_regulated<-as.data.frame(t(output_up))
colnames(ngenes_up_regulated)<-"ngenes_up_regulated"
ngenes_up_regulated$sample<-row.names(ngenes_up_regulated)
ngenes_up_regulated$sample <- factor(ngenes_up_regulated$sample, levels = c("IL1B","IL6","IL13","IL17A","IL22","IFNA","IFNG","OSM","TGFB1","TNF"))
ngenes_up_regulated %>% ggplot (aes(x=sample, y=ngenes_up_regulated)) + 
  geom_bar(stat="identity")+
  theme_classic()

ngenes_down_regulated<-as.data.frame(t(output_down))
colnames(ngenes_down_regulated)<-"ngenes_down_regulated"
ngenes_down_regulated$sample<-row.names(ngenes_down_regulated)
ngenes_down_regulated$sample <- factor(ngenes_down_regulated$sample, levels = c("IL1B","IL6","IL13","IL17A","IL22","IFNA","IFNG","OSM","TGFB1","TNF"))
ngenes_down_regulated %>% ggplot (aes(x=sample, y=ngenes_down_regulated)) + 
  geom_bar(stat="identity")+
  theme_classic()

#IL-6 and IL-22 have very little effect, # of regulated genes are in the range of false positives

##plotting comparative heatmaps with all genes regulated by any cytokine##
IL1B_reg<-IL1B %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame() 
IL6_reg<-IL6 %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()
IL13_reg<-IL13 %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()
IL17A_reg<-IL17A %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()
IL22_reg<-IL22 %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame() 
IFNA_reg<-IFNA %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame() 
IFNG_reg<-IFNG %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()
OSM_reg<-OSM %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame() 
TGFB1_reg<-TGFB1 %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()  
TNF_reg<-TNF %>% filter(padj<0.05 & (abs(log2FoldChange))>1) %>% as.data.frame()

list<-list(IL1B_reg,IL6_reg,IL13_reg,IL17A_reg,IL22_reg,IFNA_reg,IFNG_reg,OSM_reg,TGFB1_reg,TNF_reg)
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
         legend_breaks = seq(-5,5, by=1),
         border_color=NA,
         #cutree_rows = 5,
         gaps_col = seq(4,40,by=4),
         cellwidth = 2.5,
         cellheight = 0.1,
         #cluster_rows=FALSE,
         cluster_cols=FALSE,
         fontsize_row = 1,
         fontsize_col = 2.5,
         clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         )
###> MAKE SUPPLEMENT FIGURE? ###

## assessing impact of cytokines on top 'top m' markers of each cluster ##
## using eigengene scores - this is redundant to what we have already, but a sanity check ##
topm=100
clus_markers<-read.csv("scRNAseq_markers.csv", header=T, stringsAsFactors = F) %>% 
  filter(pvals_adj<0.05) %>% arrange(group,desc(logfoldchanges)) %>% group_by(group) %>% slice_head(n=topm)
clusters<-group_split(clus_markers)

MF1top<-clusters[[1]] %>% slice_head(n=topm) %>% pull (group, names)  %>% as.data.frame() %>% rename(.="color")
MF2top<-clusters[[2]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
PCtop<-clusters[[3]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
S1top<-clusters[[4]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
S2top<-clusters[[5]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
S3top<-clusters[[6]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
S3xtop<-clusters[[7]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
S4top<-clusters[[8]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")
S5top<-clusters[[9]] %>% slice_head(n=topm) %>% pull (group, names) %>% as.data.frame() %>% rename(.="color")

MF1<-merge(MF1top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
MF1col<- MF1[,1:2] %>% column_to_rownames(var="Row.names")
MF2<-merge(MF2top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
MF2col<- MF2[,1:2] %>% column_to_rownames(var="Row.names")
PC<-merge(PCtop,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
PCcol<- PC[,1:2] %>% column_to_rownames(var="Row.names")
S1<-merge(S1top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
S1col<- S1[,1:2] %>% column_to_rownames(var="Row.names")
S2<-merge(S2top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
S2col<- S2[,1:2] %>% column_to_rownames(var="Row.names")
S3<-merge(S3top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
S3col<- S3[,1:2] %>% column_to_rownames(var="Row.names")
S3x<-merge(S3xtop,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
S3xcol<- S3x[,1:2] %>% column_to_rownames(var="Row.names")
S4<-merge(S4top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
S4col<- S4[,1:2] %>% column_to_rownames(var="Row.names")
S5<-merge(S5top,reg_genes_matrix,by.x=0,by.y=0,all.y=TRUE)
S5col<- S5[,1:2] %>% column_to_rownames(var="Row.names")

ME_expr<-t(reg_genes_matrix)

ME_MF1 = moduleEigengenes(ME_expr, colors = MF1col$color)
ME_MF2 = moduleEigengenes(ME_expr, colors = MF2col$color)
ME_S1 = moduleEigengenes(ME_expr, colors = S1col$color)
ME_S2 = moduleEigengenes(ME_expr, colors = S2col$color)
ME_S3 = moduleEigengenes(ME_expr, colors = S3col$color)
ME_S3x = moduleEigengenes(ME_expr, colors = S3xcol$color)
ME_S4 = moduleEigengenes(ME_expr, colors = S4col$color)
ME_S5 = moduleEigengenes(ME_expr, colors = S5col$color)

MEscore_MF1<-ME_MF1[["eigengenes"]] %>% as.data.frame()
MEscore_MF2<-ME_MF2[["eigengenes"]] %>% as.data.frame()
MEscore_S1<-ME_S1[["eigengenes"]] %>% as.data.frame()
MEscore_S2<-ME_S2[["eigengenes"]] %>% as.data.frame()
MEscore_S3<-ME_S3[["eigengenes"]] %>% as.data.frame()
MEscore_S3x<-ME_S3x[["eigengenes"]] %>% as.data.frame()
MEscore_S4<-ME_S4[["eigengenes"]] %>% as.data.frame()
MEscore_S5<-ME_S5[["eigengenes"]] %>% as.data.frame()

#list<-list(MEscore_MF1,MEscore_MF2,MEscore_S1,MEscore_S2,MEscore_S3,MEscore_S3x,MEscore_S4,MEscore_S5)
#test<-Reduce(function(x,y) merge (x,y), list)

one<-merge(MEscore_MF1,MEscore_MF2,by=0)
two<-merge(one, MEscore_S1,by.x=1,by.y=0)
three<-merge(two, MEscore_S2,by.x=1,by.y=0)
four<-merge(three, MEscore_S3,by.x=1,by.y=0)
five<-merge(four, MEscore_S3x,by.x=1,by.y=0)
six<-merge(five, MEscore_S4,by.x=1,by.y=0)
seven<-merge(six, MEscore_S5,by.x=1,by.y=0)

MEscores_cyt_clusters<-t(seven)
colnames(MEscores_cyt_clusters)<-MEscores_cyt_clusters[1,]
MEscores_cyt_clusters<-as.data.frame(MEscores_cyt_clusters[-1,])
rn<-row.names(MEscores_cyt_clusters)
MEscores_cyt_clusters<-sapply(MEscores_cyt_clusters[,], as.numeric)
row.names(MEscores_cyt_clusters)<-rn
colnames(MEscores_cyt_clusters)<-sample_info$stimulation
pheatmap(MEscores_cyt_clusters,
         scale="row",
         legend_breaks = seq(4,-4),
         border_color=NA,
         cutree_rows = 3,
         #gaps_col = seq(4,40,by=4),
         cellwidth = 10,
         cellheight = 10,
         #cluster_rows=FALSE,
         cluster_cols=FALSE,
         fontsize_row = 10,
         fontsize_col = 10,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
)

##projecting cytokine responses onto SCENIC+ analysis##
scenic_TFs<-read.csv("SCENIC_nodes_TF.csv",header=T,stringsAsFactors = T)
colnames(reg_genes_matrix)
sample_info$id
colnames(reg_genes_matrix)<-sample_info$stimulation
reg_TF_matrix<-reg_genes_matrix[which((row.names(reg_genes_matrix)) %in% scenic_TFs$TF),]

pheatmap(reg_TF_matrix,
         scale="row",
         legend_breaks = seq(4,-4),
         border_color=NA,
         cutree_rows = 3,
         #gaps_col = seq(4,40,by=4),
         cellwidth = 10,
         cellheight = 10,
         #cluster_rows=FALSE,
         cluster_cols=FALSE,
         fontsize_row = 10,
         fontsize_col = 10,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
)
#none of the cluster top TFs stand out as being significantly regulated, except IFN-induced ones
#moving on with regulated target genes instead of TFs
scenic_targets<-read.csv("SCENIC_nodes_gene.csv",header=T,stringsAsFactors = T)
scenic_targets_S5<-scenic_targets %>% filter(GEX_celltype_Log2FC_S5>1) %>% pull(Gene)
reg_target_matrix_S5<-reg_genes_matrix[which((row.names(reg_genes_matrix)) %in% scenic_targets_S5),]
reg_target_matrix_S5<-reg_target_matrix_S5[,order(colnames(reg_target_matrix_S5))]

pheatmap(reg_target_matrix_S5,
         scale="row",
         legend_breaks = seq(4,-4),
         border_color=NA,
         cutree_rows = 3,
         #gaps_col = seq(4,40,by=4),
         cellwidth = 10,
         cellheight = 10,
         #cluster_rows=FALSE,
         cluster_cols=FALSE,
         fontsize_row = 10,
         fontsize_col = 10,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
)

