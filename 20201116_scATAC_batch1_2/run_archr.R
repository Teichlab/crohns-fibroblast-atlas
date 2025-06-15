## scATAC analysis with ArchR ##
# Load packages ####
loadPackage(tidyverse)
loadPackage(RColorBrewer)
loadPackage(ArchR)
loadPackage(Seurat)
loadPackage(SingleCellExperiment)
set.seed(1)
addArchRThreads(8)
addArchRGenome('hg38')

# Define constants ####
sampleNames <- sprintf('A%02d', seq(1, 13))
ArrowFiles <- dir('scatac/ArrowFiles', full.names=T)
projName <- 'all'
seRNA_h5ad <- '../20200715_scRNA_batch1_2/scRNA_batch1_2.downsampled.raw.h5ad'

# Load data ####
## Load ArchR project
proj <- loadArchRProject(path=projName)
## Load sample metadata
meta_df <- read_csv('../data/clinical_data.csv')
## Merge metadata
proj@cellColData$Sample@values <- sprintf('A%02d', 1:13)
orig_cellColData <- rownames(proj@cellColData)
obs_df <- proj@cellColData %>% merge(meta_df, by.x='Sample', by.y='ATAC_SID')
rownames(obs_df) <- rownames(orig_cellColData)
proj@cellColData <- obs_df
## Load scRNA
seRNA <- ReadH5AD(seRNA_h5ad)
Seurat::NormalizeData(seRNA)
Seurat::FindVariableFeatures(seRNA)
seRNA@meta.data$annot1 <- factor(
  seRNA@meta.data$annot1,
  levels = c(
    'S1', 'S2', 'S3', 'S3x', 'S4', 'S5', 'MF1', 'MF2', 'PC', 'GL', 'IM', 'LE',
    'VE arterial', 'VE arterial tip', 'VE capillary', 'VE postcapillary venlue'
  )
)
seRNA <- Seurat::SetIdent(seRNA, value='annot1')

# Filter, Dimensionality-reduction, Batch-correction, Clustering ####
## find doublets
doubleScores <- addDoubletScores(
  input=ArrowFiles,
  k=10,
  knnMethod='UMAP',
  LSIMethod=1
)
## create project
proj <- ArchRProject(
  ArrowFiles=ArrowFiles,
  outputDirectory=projName,
  copyArrows=T,
  showLogo=F
)
## filter doublets
proj <- filterDoublets(ArchRProj=proj)
## dimensionality reduction
proj <- addIterativeLSI(
  ArchRProj=proj,
  useMatrix="TileMatrix",
  name="IterativeLSI"
)
## UMAP
proj <- addUMAP(ArchRProj=proj, reducedDims="IterativeLSI")
## clustering
proj <- addClusters(
  input=proj,
  reducedDims="IterativeLSI",
  method='Seurat',
  name='Clusters',
  resolution=0.8
)
## alternative DR, UMAP and clustering with batch correction
proj <- addHarmony(
  ArchRProj=proj,
  reducedDims="IterativeLSI",
  name="Harmony",
  groupBy="Sample"
)
proj <- addUMAP(ArchRProj=proj, reducedDims="Harmony", name="UMAP_hm")
proj <- addClusters(
  input=proj,
  reducedDims="Harmony",
  method='Seurat',
  name='Clusters_hm',
  resolution=0.7
)
## imputation
proj <- addImputeWeights(proj, k=5)
## Save project
proj <- saveArchRProject(ArchRProj=proj)

# Export gene score matrix ####
gsm <- getMatrixFromProject(proj, 'GeneScoreMatrix')
saveRDS(gsm, file='scATAC_GeneScoreMatrix.rds')
rm(gsm)

# Label transfer from scRNA ####
## without harmony
### Unconstrained label transfer
proj <- addGeneIntegrationMatrix(
  proj,
  seRNA=seRNA,
  addToArrow=F,
  force=T,
  useMatrix='GeneScoreMatrix',
  matrixName='GeneIntegrationMatrix',
  reducedDims='IterativeLSI',
  groupRNA='annot1',
  nameCell='predictedCell_Un',
  nameGroup='predictedGroup_Un',
  nameScore='predictedScore_Un'
)
### Constrained label transfer
cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
cM1 <- t(t(cM) / colSums(cM)) / rowSums(t(t(cM) / colSums(cM)))
preClust <- colnames(cM1)[apply(cM1, 1, which.max)]
clustFB <- rownames(cM)[grepl('^S[0-9]', preClust)]
clustVE <- rownames(cM)[grepl('^VE', preClust)]
clustOther <- rownames(cM)[!grepl('^(S[0-9]|VE)', preClust)]
rnaFB <- colnames(seRNA)[grepl('^S[0-9]', seRNA@meta.data$annot1)]
rnaVE <- colnames(seRNA)[grepl('^VE', seRNA@meta.data$annot1)]
rnaOther <- colnames(seRNA)[!grepl('^(S[0-9]|VE)', seRNA@meta.data$annot1)]
groupList <- SimpleList(
  FB=SimpleList(
    ATAC=proj$cellNames[proj$Clusters %in% clustFB],
    RNA=rnaFB
  ),
  VE=SimpleList(
    ATAC=proj$cellNames[proj$Clusters %in% clustVE],
    RNA=rnaVE
  ),
  Other=SimpleList(
    ATAC=proj$cellNames[proj$Clusters %in% clustOther],
    RNA=rnaOther
  )
)
proj <- addGeneIntegrationMatrix(
  proj,
  seRNA=seRNA,
  groupList=groupList,
  force=T,
  useMatrix='GeneScoreMatrix',
  matrixName='GeneIntegrationMatrix',
  reducedDims='IterativeLSI',
  addToArrow=F,
  groupRNA='annot1',
  groupATAC='Clusters',
  nameCell='predictedCell_Co',
  nameGroup='predictedGroup_Co',
  nameScore='predictedScore_Co'
)

## with harmony
### Unconstrained label transfer
proj <- addGeneIntegrationMatrix(
  proj,
  seRNA=seRNA,
  addToArrow=F,
  force=T,
  useMatrix='GeneScoreMatrix',
  matrixName='GeneIntegrationMatrix',
  reducedDims='Harmony',
  groupRNA='annot1',
  nameCell='predictedCell_Un_hm',
  nameGroup='predictedGroup_Un_hm',
  nameScore='predictedScore_Un_hm'
)
### Constrained label transfer
cM_hm <- as.matrix(confusionMatrix(proj$Clusters_hm, proj$predictedGroup_Un_hm))
cM1_hm <- t(t(cM_hm) / colSums(cM_hm)) / rowSums(t(t(cM_hm) / colSums(cM_hm)))
preClust_hm <- colnames(cM1_hm)[apply(cM1_hm, 1, which.max)]
clustFB_hm <- rownames(cM_hm)[grepl('^S[0-9]', preClust_hm)]
clustVE_hm <- rownames(cM_hm)[grepl('^VE', preClust_hm)]
clustOther_hm <- rownames(cM_hm)[!grepl('^(S[0-9]|VE)', preClust_hm)]
groupList_hm <- SimpleList(
  FB=SimpleList(
    ATAC=proj$cellNames[proj$Clusters_hm %in% clustFB_hm],
    RNA=rnaFB
  ),
  VE=SimpleList(
    ATAC=proj$cellNames[proj$Clusters_hm %in% clustVE_hm],
    RNA=rnaVE
  ),
  Other=SimpleList(
    ATAC=proj$cellNames[proj$Clusters_hm %in% clustOther_hm],
    RNA=rnaOther
  )
)
proj <- addGeneIntegrationMatrix(
  proj,
  seRNA=seRNA,
  groupList=groupList_hm,
  useMatrix='GeneScoreMatrix',
  matrixName='GeneIntegrationMatrix',
  reducedDims='Harmony',
  addToArrow=F,
  force=T,
  groupRNA='annot1',
  nameCell='predictedCell_Co_hm',
  nameGroup='predictedGroup_Co_hm',
  nameScore='predictedScore_Co_hm'
)

## Choose this one to add
proj <- addGeneIntegrationMatrix(
  proj,
  seRNA=seRNA,
  groupList=groupList_hm,
  useMatrix='GeneScoreMatrix',
  matrixName='GeneIntegrationMatrix',
  reducedDims='Harmony',
  addToArrow=T,
  force=T,
  groupRNA='annot1',
  nameCell='predictedCell',
  nameGroup='predictedGroup',
  nameScore='predictedScore'
)

# Annotate clusters ####
cM <- confusionMatrix(proj$Clusters, proj$predictedGroup)
cM1 <- t(t(cM) / colSums(cM)) / rowSums(t(t(cM) / colSums(cM)))
labelOld <- rownames(cM1)
labelNew <- colnames(cM)[apply(cM1, 1, which.max)]
proj$Annot <- mapLabels(proj$Clusters, newLabels=labelNew, oldLabels=labelOld)

# Call peaks by cluster ####
## Make pseudobulk
## Rename "Annot" to "Clusters", as addReproduciblePeakSet seems to only work
## with "Clusters"
proj$Clusters1 <- proj$Clusters
proj$Clusters <- proj$Annot
proj2 <- addGroupCoverages(ArchRProj=proj, groupBy="Clusters")

## Call peaks using macs2
pathToMacs2 <- findMacs2()
proj2 <- addReproduciblePeakSet(
  proj2, groupby='Clusters', pathToMacs2=pathToMacs2
)
proj2 <- loadArchRProject(path=projName)
proj3 <- addPeakMatrix(proj2)
saveArchRProject(ArchRProj=proj3, load=F)

# Export peak matrix ####
pm <- getMatrixFromProject(proj3, useMatrix='PeakMatrix')
saveRDS(pm, file='scATAC_PeakMatrix.rds')
rm(pm)
peakSetDT <- as.data.table(proj3@peakSet)
peakSetDT$peakName <- paste0(
  peakSetDT$seqnames, ':', peakSetDT$start, '-', peakSetDT$end, ',',
  peakSetDT$nearestGene, ',', peakSetDT$peakType
)
write_tsv(peakSetDT, 'scATAC_PeakSet.tsv')

# Link peak-to-gene linkage ####
proj3 <- addPeak2GeneLinks(ArchRProj=proj3, reducedDims='Harmony')
p2g <- getPeak2GeneLinks(
  ArchRProj=proj3,
  corCutOff=0.45,
  resolution=1,
  returnLoops=F
)

# Find marker peaks for each cluster ####
markersPeaks <- getMarkerFeatures(
  ArchRProj=proj3,
  useMatrix="PeakMatrix",
  groupBy="Clusters",
  bias=c("TSSEnrichment", "log10(nFrags)"),
  testMethod="wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff="FDR <= 0.01 & Log2FC >= 1", returnGR=T)

celltypes <- names(table(proj3@cellColData$Annot))
for (ct in celltypes) {
  if (ct %in% names(markerList)) {
    export.bed(
      markerList[[ct]],
      con=paste0(
        '../20201116_scATAC_batch1_2/all/markerPeaks/',
        ct,'.fdr_0_01.lfc_1.bed'
      )
    )
  }
}

# DE peak between inflamed non-inflamed ####
dePeak <- getMarkerFeatures(
  proj3,
  useMatrix='PeakMatrix',
  groupBy='state',
  bias=c("TSSEnrichment", "log10(nFrags)"),
  testMethod="wilcoxon"
)
dePeakList <- getMarkers(
  dePeak, cutOff="FDR <= 0.05 & Log2FC >= 0.53", returnGR=T
)

# DE peak between inflamed non-inflamed within each cluster ####
proj3@cellColData$DEgroups <- paste(
  proj3@cellColData$Annot,
  proj3@cellColData$state,
  sep='-'
)
dePeakByCelltype <- list()
dePeakListByCelltype <- list()
for (ct in celltypes) {
  print(ct)
  tstGroup <- paste(ct, 'inflamed', sep='-')
  bgdGroup <- paste(ct, 'non-inflamed', sep='-')
  if (sum(proj3@cellColData$DEgroups == tstGroup) < 10 ||
      sum(proj3@cellColData$DEgroups == bgdGroup) < 10) {
    message('too few cells in group')
    next
  }
  dePeakByCelltype[[ct]] <- getMarkerFeatures(
    proj3,
    useMatrix='PeakMatrix',
    groupBy='DEgroups',
    useGroups=tstGroup,
    bgdGroups=bgdGroup,
    bias=c("TSSEnrichment", "log10(nFrags)"),
    testMethod="wilcoxon"
  )
}
dePeakListByCelltype <- lapply(
  dePeakByCelltype,
  getMarkers,
  cutOff='Pval <= 0.001 & Log2FC >= 0.53',
  returnGR=T
)

for (ct in celltypes) {
  if (ct %in% names(dePeakListByCelltype)) {
    export.bed(
      dePeakListByCelltype[[ct]][[1]],
      con=paste0(
        '../20201116_scATAC_batch1_2/all/dePeaks/',
        ct,'.inflamed.p_0_001.lfc_0_53.bed'
      )
    )
  }
}

# DE gene between inflamed non-inflamed within each cluster ####
deGeneByCelltype <- list()
deGeneListByCelltype <- list()
for (ct in celltypes) {
  print(ct)
  tstGroup <- paste(ct, 'inflamed', sep='-')
  bgdGroup <- paste(ct, 'non-inflamed', sep='-')
  if (sum(proj3@cellColData$DEgroups == tstGroup) < 10 ||
      sum(proj3@cellColData$DEgroups == bgdGroup) < 10) {
    message('too few cells in group')
    next
  }
  deGeneByCelltype[[ct]] <- getMarkerFeatures(
    proj3,
    useMatrix='GeneScoreMatrix',
    groupBy='DEgroups',
    useGroups=tstGroup,
    bgdGroups=bgdGroup,
    bias=c("TSSEnrichment", "log10(nFrags)"),
    testMethod="wilcoxon"
  )
  deGeneListByCelltype[[ct]] <-
    getMarkers(deGeneByCelltype[[ct]], cutOff="FDR <= 0.05 & Log2FC >= 0.53")
}

# Make coverage tracks ####
proj3@cellColData$AnnotBySample <- paste(
  proj3@cellColData$Annot,
  proj3@cellColData$Sample,
  sep='-'
)
getGroupBW(ArchRProj=proj3, groupBy='Sample')
getGroupBW(ArchRProj=proj3, groupBy='Annot')
getGroupBW(ArchRProj=proj3, groupBy='DEgroups')
getGroupBW(ArchRProj=proj3, groupBy='AnnotBySample')


# Add motif annotation ####
proj3 <- addMotifAnnotations(proj3, motifSet='cisbp', name='Motif')

motifsUp <- peakAnnoEnrichment(
  seMarker=dePeakByCelltype[['S4']],
  ArchRProj=proj3,
  peakAnnotation="Motif",
  cutOff="FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF=rownames(motifsUp), mlog10Padj=assay(motifsUp)[, 1])
df <- df[order(df$mlog10Padj, decreasing=TRUE), ]
df$rank <- seq_len(nrow(df))




# Trajectory ####


# Plot ####
p1 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="Sample",
  embedding="UMAP"
)
p2 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="Clusters",
  embedding="UMAP"
)
p3 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="predictedGroup_Un",
  embedding="UMAP"
)
p4 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="predictedGroup_Co",
  embedding="UMAP"
)
plotPDF(
  plotList=list(p1, p2, p3, p4),
  name="Plot-UMAP.pdf",
  ArchRProj=proj,
  addDOC=FALSE,
  width=5,
  height=5
)

p5 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="Sample",
  embedding="UMAP_hm"
)
p6 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="Clusters_hm",
  embedding="UMAP_hm"
)
p7 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="predictedGroup_Un_hm",
  embedding="UMAP_hm"
)
p8 <- plotEmbedding(
  ArchRProj=proj,
  colorBy="cellColData",
  name="predictedGroup_Co_hm",
  embedding="UMAP_hm"
)
plotPDF(
  plotList=list(p5, p6, p7, p8),
  name="Plot-UMAP-W-Harmony.pdf",
  ArchRProj=proj,
  addDOC=FALSE,
  width=5,
  height=5
)

markerGenes <- c(
  'CDH19',
  'MZB1',
  'JCHAIN',
  'PROX1',
  'HHIP',
  'DES',
  'COX4I2',
  'SFTA1P',
  'PCSK6',
  'PI16',
  'IGF1',
  'NPB',
  'PCDH11X',
  'CXCL5',
  'GJA5',
  'CA4',
  'CHST1',
  'MADCAM1'
)

p9s <- plotEmbedding(
  ArchRProj=proj,
  colorBy="GeneScoreMatrix",
  name=markerGenes,
  embedding="UMAP",
  imputeWeights=getImputeWeights(proj)
)
p9 <- lapply(p9s, function(x) {
  x + guides(color=FALSE, fill=FALSE) +
    theme_ArchR(baseSize=6.5) +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
# do.call(cowplot::plot_grid, c(list(ncol=3), p9))
plotPDF(
  plotList=p9,
  name="Plot-UMAP-Marker-Genes-W-Imputation.pdf",
  ArchRProj=proj,
  HiaddDOC=FALSE,
  width=5,
  height=5
)

p10s <- plotEmbedding(
  ArchRProj=proj,
  colorBy="GeneIntegrationMatrix",
  name=markerGenes,
  continuousSet="horizonExtra",
  embedding="UMAP",
  imputeWeights=NULL
)
p10 <- lapply(p10s, function(x) {
  x + guides(color=FALSE, fill=FALSE) +
    theme_ArchR(baseSize=6.5) +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
plotPDF(
  plotList=p10,
  name="Plot-UMAP-Marker-Genes-W-Integration.pdf",
  ArchRProj=proj,
  addDOC=FALSE,
  width=5,
  height=5
)

plotEmbedding(proj, colorBy="cellColData", name="DoubletScore")
plotEmbedding(proj3, colorBy="cellColData", name="site_taken_fine")
p11 <- plotEmbedding(proj, colorBy="cellColData", name="Annot")
plotPDF(
  p11,
  name='Plot-UMAP-Annot.pdf',
  width=5,
  height=5,
  ArchRProj=proj3,
  addDOC=F
)

# p12 ####
p12 <- plotBrowserTrack(
  proj3,
  groupBy='Annot',
  geneSymbol='CXCL5',
  features=GRanges(markerList[['S5']]),
  upstream=20000,
  downstream=20000
)
plotPDF(
  p12,
  name='Plot-Tracks-With-Features-S5-CXCL5.pdf',
  width=5,
  height=5,
  ArchRProj=proj3,
  addDOC=F
)

p13 <- plotBrowserTrack(
  proj3,
  groupBy='DEgroups',
  #useGroups=c('S5-inflamed', 'S5-non-inflamed'),
  region=flank(dePeakListByCelltype$S5$`S5-inflamed`[5:6], 50000, both=T),
  features=dePeakListByCelltype$S5$`S5-inflamed`
)
plotPDF(
  p13,
  name='Plot-Tracks-With-Features-S5-inflamed.pdf',
  width=5,
  height=5,
  ArchRProj=proj3,
  addDOC=F
)

heatmapPeaks <- markerHeatmap(
  seMarker=markersPeaks,
  cutOff="FDR <= 0.1 & Log2FC >= 0.5",
  transpose=TRUE
)
draw(
  heatmapPeaks,
  heatmap_legend_side="bot",
  annotation_legend_side="bot"
)
pma <- plotMarkers(
  seMarker=markersPeaks,
  name="MF2",
  cutOff="FDR <= 0.1 & Log2FC >= 1",
  plotAs="MA"
)
pma






# p10s <- plotEmbedding(
#   ArchRProj=proj, colorBy="GeneScoreMatrix", name=markerGenes,
#   embedding="UMAP_hm", imputeWeights=getImputeWeights(proj))
# p10 <- lapply(p10s, function(x) {
#   x + guides(color=FALSE, fill=FALSE) +
#     theme_ArchR(baseSize=6.5) +
#     theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
#     theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())
# })
# plotPDF(
#   plotList=p10, name="Plot-UMAP-Marker-Genes-W-Harmony-Imputation.pdf", ArchRProj=proj,
#   addDOC=FALSE, width=5, height=5)
#
# p11s <- plotEmbedding(
#   ArchRProj=proj,
#   colorBy="GeneIntegrationMatrix",
#   name=markerGenes,
#   continuousSet="horizonExtra",
#   embedding="UMAP_hm",
#   imputeWeights=NULL
#   # imputeWeights=getImputeWeights(proj)
# )
# p11 <- lapply(p11s, function(x){
#   x + guides(color=FALSE, fill=FALSE) +
#     theme_ArchR(baseSize=6.5) +
#     theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
#     theme(
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank(),
#       axis.text.y=element_blank(),
#       axis.ticks.y=element_blank()
#     )
# })
# plotPDF(
#   plotList=p11, name="Plot-UMAP-Marker-Genes-W-Harmony-Integration.pdf",
#   ArchRProj=proj, addDOC=FALSE, width=5, height=5)


# df <- data.frame(table(as.vector(proj@cellColData$Sample), proj@cellColData$predictedGroup_Un))
# colnames(df) <- c('Sample', 'CellType', 'Count')
# pdf(file=paste0(projName, '/scATAC_projected_celltype_composition_by_sample.pdf'))
# print(df %>% ggplot(aes(x=Sample, y=Count, fill=CellType)) + geom_bar(position='fill', stat='identity') +
#         scale_fill_manual(values=c(brewer.pal(15, 'Paired'), grey.colors(1))))
# dev.off()
