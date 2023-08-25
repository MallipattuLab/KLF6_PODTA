
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(biovizBase)
set.seed(1234)


counts <- Read10X_h5(filename = "SC393/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "SC393/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'SC393/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

WT_SHAM <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
WT_SHAM$OE <- "WTSHAM"


counts <- Read10X_h5(filename = "SC391/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "SC391/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'SC391/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

OE_SHAM <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
OE_SHAM$OE <- "OESHAM"

counts <- Read10X_h5(filename = "SC394/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "SC394/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'SC394/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

WT_STZ <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
WT_STZ$OE <- "WTSTZ"

counts <- Read10X_h5(filename = "SC392/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "SC392/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'SC392/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

OE_STZ <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
OE_STZ$OE <- "OESTZ"

#######Merging the groups######
Data.Combined <- merge(WT_SHAM, y = c(OE_SHAM, WT_STZ, OE_STZ), add.cell.ids = c("WT_SHAM", "OE_SHAM", "WT_STZ", "OE_STZ"), project = "Combined")
###############################

Data.Combined
Data.Combined[['peaks']]
granges(Data.Combined)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(Data.Combined) <- annotations

saveRDS(Data.Combined, file = "preQC_ATAC.rds")

##########################################Perform QC##################
# compute nucleosome signal score per cell
Data.Combined <- NucleosomeSignal(object = Data.Combined)

# compute TSS enrichment score per cell
Data.Combined <- TSSEnrichment(object = Data.Combined, fast = FALSE)
head(Data.Combined@meta.data)

# add blacklist ratio and fraction of reads in peaks
Data.Combined$pct_reads_in_peaks <- Data.Combined$peak_region_fragments / Data.Combined$passed_filters * 100
Data.Combined$blacklist_ratio <- Data.Combined$blacklist_region_fragments / Data.Combined$peak_region_fragments

Data.Combined$high.tss <- ifelse(Data.Combined$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(Data.Combined, group.by = 'high.tss') + NoLegend()

Data.Combined$nucleosome_group <- ifelse(Data.Combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Data.Combined, group.by = 'nucleosome_group')
FragmentHistogram(object = Data.Combined, region = "chr1-1-5000000")

head(Data.Combined@meta.data)

VlnPlot(
  object = Data.Combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
saveRDS(Data.Combined, file = "preSubset_ATAC.rds")

Data.Combined <- subset(
  x = Data.Combined,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 5 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)
###########################
Data.Combined.list <- SplitObject(Data.Combined, split.by = "OE")
Data.Combined.list <- lapply(X = Data.Combined.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = Data.Combined.list)
Data.Combined.list <- lapply(X = Data.Combined.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunSVD(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = Data.Combined.list, reduction = "rlsi",
                                  dims = 1:50)
Data.Combined <- IntegrateData(anchorset = anchors, dims = 1:30)
Data.Combined <- RunSVD(Data.Combined)
Data.Combined <- RunUMAP(object = Data.Combined, reduction = 'lsi', dims = 2:30)
Data.Combined <- FindNeighbors(object = Data.Combined, reduction = 'lsi', dims = 2:30)
Data.Combined <- FindClusters(object = Data.Combined, verbose = FALSE, algorithm = 3)
DimPlot(object = Data.Combined, label = TRUE) + NoLegend()
DimPlot(object = Data.Combined, group.by = "OE", label = F)
###########################

Data.Combined <- RunTFIDF(Data.Combined)
Data.Combined <- FindTopFeatures(Data.Combined, min.cutoff = 'q0')
Data.Combined <- RunSVD(Data.Combined)

memory.limit(size = 100000)
options(future.globals.maxSize= 15912896000)
####

DepthCor(Data.Combined)

gene.activities <- GeneActivity(Data.Combined)
Data.Combined[['GA']] <- CreateAssayObject(counts = gene.activities)
head(Data.Combined@meta.data)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
Data.Combined <- NormalizeData(
  object = Data.Combined,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(Data.Combined$nCount_GA)
)

DefaultAssay(Data.Combined) <- 'GA'

saveRDS(Data.Combined, file = "Data.Combined_integrated.rds")

FeaturePlot(
  object = Data.Combined,
  features = c('Nphs2'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
markers.to.plot<-c("Nphs1", "Flt1", "Cfh", "Cd44", "Top2a", "Ptprc", "Insrr", "Clnk", "Frmpd4", "Egfem1", "Trpm6", "Slc12a3", "Slc12a1", "Cyp7b1", "Keg1", "Cntnap5a", "Slc5a12", "Tmem27")
DotPlot(Data.Combined, features = rev(markers.to.plot), cols = c("yellow", "red"), 
        dot.scale = 8) + RotatedAxis()
DotPlot(Data.Combined, features = rev(markers.to.plot), split.by = "OE",  cols = c("yellow", "red", "blue", "green"),
        dot.scale = 8) + RotatedAxis()

##############################################################################################################################
# Load the pre-processed snRNA-seq data 
K6_rna <- Data.Combined1
DefaultAssay(K6_rna) <- "integrated"

transfer.anchors <- FindTransferAnchors(
  reference = K6_rna,
  query = Data.Combined,
  reduction = 'cca'
)

K6_rna$celltype <- Idents(K6_rna)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = K6_rna$celltype,
  weight.reduction = Data.Combined[['lsi']],
  dims = 2:30
)

Data.Combined <- AddMetaData(object = Data.Combined, metadata = predicted.labels)


plot1 <- DimPlot(
  object = K6_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('snRNA-seq')

plot2 <- DimPlot(
  object = Data.Combined,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle('snATAC-seq')

plot2
plot1 + plot2

table(Data.Combined@meta.data[["predicted.id"]], Data.Combined@meta.data$OE)
table(Data.Combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$OE)
plot2 <- DimPlot(object = Data.Combined,  group.by = 'predicted.id', split.by = "OE",  label = TRUE,
  repel = TRUE) + ggtitle('snATAC-seq')


install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
remotes::install_github('satijalab/seurat-wrappers')

library(cicero)
library(SeuratWrappers)
# convert to CellDataSet format and make the cicero object
bone.cds <- as.cell_data_set(x = Data.Combined)
bone.cicero <- make_cicero_cds(bone.cds, reduced_coordinates = reducedDims(bone.cds)$UMAP)

DefaultAssay(Data.Combined) <- "peaks"
genome <- seqlengths(Data.Combined)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(bone.cicero, genomic_coords = genome.df, sample_num = 100)
head(genome)
head(conns)
ccans <- generate_ccans(conns)
head(ccans)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(Data.Combined) <- links



DefaultAssay(Data.Combined) <- 'peaks'
Idents(Data.Combined) <- Data.Combined$predicted.id 
CoveragePlot(Data.Combined, region = 'Acss2', features = 'Acss2', assay = 'peaks', expression.assay = 'RNA', peaks = FALSE)
immune.combined <- subset(Data.Combined, idents = c("Novel(PT)", "PT(S3)", "PT(S1-S2)", "PT(S3)/LH(DL)", "PT(S3)-2"))
head(immune.combined@meta.data)
Idents(immune.combined) <- "celltype.OE"
saveRDS(immune.combined, "PTs_subset.rds")
Idents(immune.combined) <- "predicted.id"
Novel <- subset(immune.combined, idents = "Novel(PT)")
PTS3 <- subset(immune.combined, idents = "PT(S3)")



CoveragePlot(Data.Combined, region = 'Nphs1', features = 'Nphs1', assay = 'peaks', expression.assay = 'GA', peaks = FALSE)
CoveragePlot(Data.Combined, region = 'Clu', extend.upstream = 40000, extend.downstream = 20000)
CoveragePlot(PodSHAM, region = 'Clu', extend.upstream = 20000, extend.downstream = 10000)


immune.combined <- subset(Data.Combined, idents = c("PT(S3)", "Novel(PT)", "PT(S1-S2)", "PT(S3)-2", "PT(S3)/LH(DL)"))
immune.combined <- RenameIdents(object = immune.combined, "PT(S3)" = "PT(S3)", "Novel(PT)" = "Novel(PT)", "PT(S1-S2)" = "PT(S1-S2)", "PT(S3)-2" = "PT(S3)", "PT(S3)/LH(DL)" = "PT(S3)-LH(DL)")
CoveragePlot(immune.combined, region = 'Klf15', extend.upstream = 40000, extend.downstream = 20000)
head(Data.Combined@meta.data)
Data.Combined$celltype.OE <- paste(Idents(Data.Combined), Data.Combined$celltype, sep = '_')
saveRDS(Data.Combined, "AfterCicero_9_2_21.rds")

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
Data.Combined <- AddMotifs(
  object = Data.Combined,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
#Make a subset of the podocyte and find uperegulated peaks
#Set idents to comparison feature, ie OESHAM vs WTSHAM
da_peaks <- FindMarkers(
  object = Data.Combined,
  ident.1 = 'OESHAM',
  ident.2 = 'WTSHAM',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])


#Plot podocyte cluster Differentially accessible regions
# find cell type specific TF and make a unique list
#make a subset of the podocyte cluster
tf_celltype.df <- FindMarkers(PodClus, ident.1 = 'OESHAM', ident.2 = 'WTSHAM', only.pos = T, test.use = 'LR', min.pct = 0.05, 
                              latent.vars = 'nCount_peaks')
top.da.peak <- rownames(tf_celltype.df[tf_celltype.df$p_val < 0.005, ])
enriched.motifs <- FindMotifs(
  object = Data.Combined_ATAC,
  features = top.da.peak
)
#save enriched.motifs file
write.csv(enriched.motifs, "PodMotifs_OESHAMvWTSHAM.csv")

top.chromvar <- rownames(enriched.motifs[enriched.motifs$pvalue < 0.05, ])
aver_chromvar <- AverageExpression(PodClus, assays = "chromvar", features = top.chromvar) %>%
  as.data.frame()

# visualize results #change scaling as necessary for better visualization
fig <- pheatmap::pheatmap(aver_chromvar,scale = "row",
                          cluster_cols=F,cluster_rows = F,
                          color = jdb_palette("brewer_yes"),
                          show_rownames=F) #450x640


saveRDS(fig1, "Motif_Pod_OESHAMvsWTSHAM.rds")


head(Data.Combined_ATAC@meta.data)
scoring <- as.data.frame(Data.Combined_ATAC@meta.data)
write.csv(scoring, "metadata_ATAC_for_scoring.csv")
scoring <- subset(scoring, select = -c(x,y) )
#Get the columns needed to plot the heatmap (predicted.id scores and cluster names)

fig2 <- pheatmap::pheatmap(scoring1,scale = "row",
                           cluster_cols=F,cluster_rows = F,
                           color = jdb_palette("brewer_yes"),
                           show_rownames=F) #450x640
saveRDS(fig2, "scoring_ATAC_labeltransfer.rds")


