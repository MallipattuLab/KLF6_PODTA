if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
ArchR::installExtraPackages()
devtools::install_github("immunogenomics/presto")
library(ArchR)
library(presto)
addArchRThreads(threads = 2) 
addArchRGenome("mm10")

inputFiles_391<-c("/gpfs/projects/MallipattuGroup/Neha/SC39/SC391/fragments.tsv.gz")
inputFiles_391

inputFiles_392<-c("/gpfs/projects/MallipattuGroup/Neha/SC39/SC392/fragments.tsv.gz")
inputFiles_392

inputFiles_393<-c("/gpfs/projects/MallipattuGroup/Neha/SC39/SC393/fragments.tsv.gz")
inputFiles_393

inputFiles_394<-c("/gpfs/projects/MallipattuGroup/Neha/SC39/SC394/fragments.tsv.gz")
inputFiles_394

ArrowFiles_391<-c("E:/Neha/OESHAM.arrow")
ArrowFiles_392<-c("E:/Neha/OE.arrow")
ArrowFiles_393<-c("E:/Neha/OE.arrow")
ArrowFiles_394<-c("E:/Neha/WT.arrow")

ArrowFiles_391 <- createArrowFiles(
  inputFiles = inputFiles_391,
  sampleNames = c("OESHAM"),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE
)


ArrowFiles_392 <- createArrowFiles(
  inputFiles = inputFiles_392,
  sampleNames = c("OE"),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE
)
ArrowFiles_393 <- createArrowFiles(
  inputFiles = inputFiles_393,
  sampleNames = c("WTSHAM"),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE
)

ArrowFiles_394 <- createArrowFiles(
  inputFiles = inputFiles_394,
  sampleNames = c("WT"),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE
)


doubScores_391 <- addDoubletScores(
  input = ArrowFiles_391,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

doubScores_392 <- addDoubletScores(
  input = ArrowFiles_392,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

doubScores_393 <- addDoubletScores(
  input = ArrowFiles_393,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

doubScores_394 <- addDoubletScores(
  input = ArrowFiles_394,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#creating an ArchRProjec
Klf6_Proj1_OESHAM <- ArchRProject(
  ArrowFiles = c(ArrowFiles_391), 
  outputDirectory = "Klf6_Proj1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

#creating an ArchRProjec
Klf6_Proj1_WTSHAM <- ArchRProject(
  ArrowFiles = c(ArrowFiles_392), 
  outputDirectory = "Klf6_Proj1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

Klf6_Proj1_WTSHAM <- ArchRProject(
  ArrowFiles = c(ArrowFiles_393), 
  outputDirectory = "Klf6_Proj1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

#creating an ArchRProjec
Klf6_Proj1_WT <- ArchRProject(
  ArrowFiles = c(ArrowFiles_394), 
  outputDirectory = "Klf6_Proj1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)



#creating an ArchRProjec
Klf6_Proj1 <- ArchRProject(
  ArrowFiles = c(ArrowFiles_391, ArrowFiles_392, ArrowFiles_393, ArrowFiles_394), 
  outputDirectory = "Klf6_Proj1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

Klf6_Proj1<-addGeneScoreMatrix(
  input = Klf6_Proj1,
  genes = getGenes(Klf6_Proj1),
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "GeneScoreMatrix",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),
  geneUpstream = 5000,
  geneDownstream = 0,
  useGeneBoundaries = TRUE,
  useTSS = FALSE,
  extendTSS = FALSE,
  tileSize = 500,
  ceiling = 4,
  geneScaleFactor = 5,
  scaleTo = 5000,
  excludeChr = c("chrY", "chrM"),
  blacklist = getBlacklist(Klf6_Proj1),
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  force = TRUE,
  logFile = createLogFile("addGeneScoreMatrix")
)

Klf6_Proj1

saveArchRProject(ArchRProj = Klf6_Proj1, outputDirectory = "Klf6_Proj1", load = FALSE)

#check how much memory it takes to store project in R
paste0("Memory Size = ", round(object.size(Klf6_Proj1) / 10^6, 3), " MB")

loadArchRProject(path = "/gpfs/projects/MallipattuGroup/Neha/01.Raw_data/Archr/Archr_Predicted/Klf6_Proj1_predictedgroup/Klf6_Proj2", force = FALSE, showLogo = TRUE)



#which matrices are available in project
getAvailableMatrices(Klf6_Proj1)

#access names associated with each cell
head(Klf6_Proj1$cellNames)

#access the sample names associated with each cell
head(Klf6_Proj1$Sample)

#access the TSS Enrichment Scores for each cell
quantile(Klf6_Proj1$TSSEnrichment)

#Make a ridge plot for each sample for the TSS enrichment scores
p1 <- plotGroups(
  ArchRProj = Klf6_Proj1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

#Make a violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
  ArchRProj = Klf6_Proj1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2

#Make a ridge plot for each sample for the log10(unique nuclear fragments)
p3 <- plotGroups(
  ArchRProj = Klf6_Proj1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p3

#Make a violin plot for each sample for the log10(unique nuclear fragments)
p4 <- plotGroups(
  ArchRProj = Klf6_Proj1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4

#plot fragment size distributions
p5 <- plotFragmentSizes(Klf6_Proj1)
p5

#plot TSS enrichment profiles
p6<-plotTSSEnrichment(ArchRProj = Klf6_Proj1, flank = 1000, smooth = 1)
p6

#Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(Klf6_Proj1_OE, select = c("log10(nFrags)", "TSSEnrichment"))
df

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p

saveRDS(Klf6_Proj1, file = "Klf6_Proj1.rds")
Klf6_Proj1<-readRDS("E:/Neha/Klf6_Proj1/Save-ArchR-Project.rds")

#filter doublets from project
Klf6_Proj2 <- filterDoublets(Klf6_Proj1)
Klf6_Proj2

#Iterative Latent Semantic Indexing
Klf6_Proj2 <- addIterativeLSI(
  ArchRProj = Klf6_Proj2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

#Batch Effect Correction with Harmony
Klf6_Proj2 <- addHarmony(
  ArchRProj = Klf6_Proj2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#Clustering using Seurat’s FindClusters() function
Klf6_Proj2 <- addClusters(
  input = Klf6_Proj2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = .5,
  force = TRUE
)

head(Klf6_Proj2$Clusters)

#tabulate the number of cells present in each cluster
table(Klf6_Proj2$Clusters)

#create a cluster confusion matrix across each sample
cM <- confusionMatrix(paste0(Klf6_Proj2$Clusters), paste0(Klf6_Proj2$Sample))
cM

#plot this confusion matrix as a heatmap
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black",
  cluster_cols = TRUE
)
p

#run UMAP in ArchR we use the addUMAP() function
Klf6_Proj2 <- addUMAP(
  ArchRProj = Klf6_Proj2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

#color by “Clusters” which were identified in a previous chapter
p7 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", baseSize = 13) + theme(legend.text = element_text(size = 13)) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
))
p7

#clusters based on sample
p8 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP") + theme(legend.text = element_text(size = 13)) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
))
p8

ggAlignPlots(p7, p8, type = "h")

#run t-SNE in ArchR we use the addTSNE() function
Klf6_Proj2 <- addTSNE(
  ArchRProj = Klf6_Proj2, 
  reducedDims = "IterativeLSI", 
  name = "TSNE", 
  perplexity = 30
)

#plot the t-SNE embedding using plotEmbedding()
p9 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p9

p10 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
p10

ggAlignPlots(p9, p10, type = "h")

#Post-Harmony Dim Reductions
Klf6_Proj2 <- addUMAP(
  ArchRProj = Klf6_Proj2, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p11 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p11
p12 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p12
ggAlignPlots(p11, p12, type = "h")

Klf6_Proj2 <- addTSNE(
  ArchRProj = Klf6_Proj2, 
  reducedDims = "Harmony", 
  name = "TSNEHarmony", 
  perplexity = 30
)


p13 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p13
p14 <- plotEmbedding(ArchRProj = Klf6_Proj2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
p14
ggAlignPlots(p13, p14, type = "h")

#gene scores
markersGS <- getMarkerFeatures(
  ArchRProj = Klf6_Proj2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(Klf6_Proj2, file = "Klf6_Proj2.rds")


saveArchRProject(ArchRProj = Klf6_Proj1, outputDirectory = "Klf6_Proj1_12.16.22", load = FALSE)


#get a list of DataFrame objects, one for each of our clusters, containing the relevant marker features using the getMarkers() function
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.4 & Log2FC >= 1.25")
markerList$C1

library("WriteXLS")
write.csv(markerList, "markerlistatac.csv")


#visualize all of the marker features simultaneously
markerGenes  <- c("Notch4","Adgre5","Pdgfrb","Kdr","Meis2", "Magi2","Ptpn14","Tmem117","Aqp2", "Slc12a1", "Klf6", "Klhl3", "Cxcl1", 
                  "Nfil3", "Atp11a", "Slc34a1","Lrp2", "Gldc", "Slc5a12", "Arhgap22", "Slc22a7", "Ngef" )

markerGenes  <- c("Clu")


heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.4 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
    transpose = TRUE,
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#Marker Genes Imputation with MAGIC
Klf6_Proj2 <- addImputeWeights(Klf6_Proj2)

p15 <- plotEmbedding(
  ArchRProj = Klf6_Proj2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Klf6_Proj2)
)
p15

#plot all the marker genes at once using cowplot, Rearrange for grid plotting
####################################################################################################################
p16 <- lapply(p11, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#######################################################################################################################


#Track Plotting with ArchRBrowser
p17 <- plotBrowserTrack(
  ArchRProj = Klf6_Proj2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  baseSize = 10,
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p17$Klf6)

saveArchRProject(ArchRProj = Klf6_Proj2, outputDirectory = "Klf6_Proj2", load = FALSE)

#Launching the ArchRBrowser
ArchRBrowser(Klf6_Proj2)
#####################################################################################################################

#Cross-platform linkage of scATAC-seq cells with scRNA-seq cells

seRNA<-readRDS("Data.Combined.6.15.21.rds")
Klf6_Proj2<-readRDS("Klf6_Proj2.rds")

p1 <- DimPlot(seRNA, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(seRNA, reduction = "umap", label = TRUE)
p2
p1

head(seRNA@meta.data)
levels(seRNA)

#new.cluster.ids <- c("PT(S3)","PT(S3)-1","PT(S3)-2","PT(S3)-3","LH(AL)","DCT","PT(S1-S2)","PEC","CD-PC","EC1","IC-B","PT(S3)-4","DCT/CNT","Pod","EC2","LH(DL)","Tcell","MC","Uro","Per/Fib")  
#names(new.cluster.ids) <- levels(seRNA)
#seRNA <- RenameIdents(seRNA, new.cluster.ids)


seRNA$groupRNA <- seRNA@active.ident


#unconstrained integration
Klf6_Proj2 <- addGeneIntegrationMatrix(
  ArchRProj = Klf6_Proj2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  seRNA = seRNA,
  addToArrow = TRUE,
  force = TRUE,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  nGenes = 5000,
  query.assay = 'RNA',
  groupRNA = "groupRNA",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
table(Klf6_Proj2$predictedGroup_Un)

#Confusion matrix btw clusters and predictedGroup_un
cM<-as.matrix(confusionMatrix(Klf6_Proj2$Clusters,Klf6_Proj2$predicatedGroup_Un))
preClust <-colnames(cM)[apply(cM,1,which.max)]
cbind(preClust,rownames(cM))

unique(unique(Klf6_Proj2$predictedGroup_Un))

saveRDS(Klf6_Proj2, file = "Klf6_Proj2integrated.rds")
saveArchRProject(ArchRProj = Klf6_Proj2integrated.rds, outputDirectory = "Klf6_Proj2integrated", load = FALSE)

#Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
Klf6_Proj3 <- addGeneIntegrationMatrix(
  ArchRProj = Klf6_Proj2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "groupRNA",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

getAvailableMatrices(Klf6_Proj3)
Klf6_Proj2 <- addImputeWeights(Klf6_Proj3)

p7 <- plotEmbedding(ArchRProj = Klf6_Proj3, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", baseSize = 13) + theme(legend.text = element_text(size = 13)) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
))
p7

markerGenes  <- c("Camk1d", "Cxcl1")

plotEmbedding(
  ArchRProj = Klf6_Proj3, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Klf6_Proj3)
)

plotEmbedding(
  ArchRProj = Klf6_Proj3, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Klf6_Proj3)
)

#Labeling scATAC-seq clusters with scRNA-seq information
cM <- confusionMatrix(Klf6_Proj3$Clusters, Klf6_Proj3$predictedGroup)
labelOld <- rownames(cM)
labelOld


#Labeling scATAC-seq using unconstrained analysis
#cM <- confusionMatrix(Klf6_Proj3$Clusters, Klf6_Proj3$predictedGroup_Un)
#labelOld <- rownames(cM)
#labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew


unique(unique(Klf6_Proj3$predictedGroup))

#UMAP showing integrated RNA label
p17 <- plotEmbedding(ArchRProj = Klf6_Proj3, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", baseSize = 13) + theme(legend.text = element_text(size = 13)) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
))
p17

#clusters based on sample
p18 <- plotEmbedding(ArchRProj = Klf6_Proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAP") + theme(legend.text = element_text(size = 13)) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
))
p18

saveRDS(Klf6_Proj3, file = "Klf6_Proj3.rds")
saveArchRProject(ArchRProj = Klf6_Proj3, outputDirectory = "Klf6_Proj3", load = FALSE)

p16 <- plotEmbedding(Klf6_Proj3, colorBy = "cellColData", name = "predictedGroup", labelMeans = FALSE, baseSize = 18) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
                                                                                                                                                             ))
p16

p16 <- plotEmbedding(Klf6_Proj3, colorBy = "cellColData", name = "Sample", labelMeans = FALSE, baseSize = 18) + theme(text=element_text(family="Arial", size=18), legend.text = element_text(size = 15
))
p16

library(BSgenome.Mmusculus.UCSC.mm10)


#Making Pseudo-bulk Replicates with integrated cluster labels ( need to install the genome)
library(BSgenome.Mmusculus.UCSC.mm10)
Klf6_Proj4 <- addGroupCoverages(ArchRProj = Klf6_Proj3, groupBy = "predictedGroup")

saveRDS(Klf6_Proj4, file = "Klf6_Proj4.rds")

#tabulate the number of cells present in each cluster
table(Klf6_Proj4$predictedGroup)
table(Klf6_Proj4$Clusters2)

cM <- confusionMatrix(paste0(Klf6_Proj4$Clusters2), paste0(Klf6_Proj4$Sample))
cM


cM <- confusionMatrix(paste0(Klf6_Proj4$predictedGroup), paste0(Klf6_Proj4$Sample))
cM

pathToMacs2 <- findMACS()

#Calling Peaks w/ Macs2
Klf6_Proj4 <- addReproduciblePeakSet(
  ArchRProj = Klf6_Proj4, 
  groupBy = "predictedGroup", 
   pathToMacs2 = "/gpfs/home/ngujarati/Conda/bin/macs2"
)

getPeakSet(Klf6_Proj4)
#for predicted group
write.csv(markersPeaks@colData, file = "Klf6_markerPeaks_predictedgroup.csv")
#for Cluster2
write.csv(markersPeaks@colData, file = "Klf6_markerPeaks.csv")

saveArchRProject(ArchRProj = Klf6_Proj4, outputDirectory = "Klf6_Proj4", load = FALSE)
loadArchRProject(path = "/gpfs/projects/MallipattuGroup/Neha/01.Raw_data/Archr/Klf6_Proj5", force = FALSE, showLogo = TRUE)

Klf6_Proj5 <- addPeakMatrix(Klf6_Proj4)

getAvailableMatrices(Klf6_Proj5)

getPeakSet(Klf6_Proj4)
write.csv(markersPeaks@colData, file = "Klf6_markerPeaks.csv")


#Identifying Marker Peaks with ArchR
#First, lets remind ourselves of the cell types that we are working with in projHeme5 and their relative proportions.
#Our scRNA labels
table(Klf6_Proj4$predictedGroup)


markersPeaks <- getMarkerFeatures(
  ArchRProj = Klf6_Proj5, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks
save.image("/gpfs/projects/MallipattuGroup/Neha/01.Raw_data/Archr/Predictedgroup/markerPeaks.RData")
write.csv(markersPeaks@colData, file = "Klf6_markerPeaks_predictedgroup.csv")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.csv(markerList@colData, file = "Klf6_markerList_predictedgroup.csv")
save.image("/gpfs/projects/MallipattuGroup/Neha/01.Raw_data/Archr/Predictedgroup/markerList.RData")

#Instead of a list of DataFrame objects, we can use getMarkers() to return a GRangesList object by setting returnGR = TRUE.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
markerList$`Novel_PT`
write.csv(markerList, file = "Klf6_markerlist_predictedgroup.csv")

# Plotting Marker Peaks in ArchR
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.001 & Log2FC >= 0.5",
  transpose = TRUE
)
save.image("/gpfs/projects/MallipattuGroup/Neha/01.Raw_data/Archr/heatmapPeaks.RData")
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#Marker Peak MA and Volcano Plots
pma <- markerPlot(seMarker = markersPeaks, name = "Pod", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma

pv <- markerPlot(seMarker = markersPeaks, name = "Pod", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv

# Marker Peaks in Browser Tracks
p <- plotBrowserTrack(
  ArchRProj = Klf6_Proj5, 
  groupBy = "predictedGroup", 
  geneSymbol = c ("Clu"),
  tileSize = 500,
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["IC-A"],
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$Clu)


# Marker Peaks in Browser Tracks each genotype
p <- plotBrowserTrack(
  ArchRProj = Klf6_Proj5, 
  groupBy = "Sample", 
  geneSymbol = c("Clu"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Pod"],
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$Clu)

saveArchRProject(ArchRProj = Klf6_Proj5, outputDirectory = "Klf6_Proj5", load = FALSE)
saveRDS(Klf6_Proj5, file = "Klf6_Proj5.rds")

#Motif Enrichment in Differential Peaks
Klf6_Proj5 <- addMotifAnnotations(ArchRProj = Klf6_Proj5, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = Klf6_Proj5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifsUp

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = Klf6_Proj5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

#Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = Klf6_Proj5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.6 & Log2FC >= 0.5"
)

enrichMotifs

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

saveArchRProject(ArchRProj = Klf6_Proj5, outputDirectory = "Klf6_Proj5", load = FALSE)

#chromVAR
#Motif Deviations
if("Motif" %ni% names(Klf6_Proj5@peakAnnotation)){
  Klf6_Proj5 <- addMotifAnnotations(ArchRProj = Klf6_Proj5, motifSet = "cisbp", name = "Motif")
}

Klf6_Proj5 <- addBgdPeaks(Klf6_Proj5)

#compute per-cell deviations accross all of our motif annotations(took more than 70mins)
Klf6_Proj5 <- addDeviationsMatrix(
  ArchRProj = Klf6_Proj5, 
  peakAnnotation = "Motif",
  force = TRUE
)

#access deviations
plotVarDev <- getVarDeviations(Klf6_Proj5, name = "MotifMatrix", plot = TRUE)

plotVarDev

#extract a subset of motifs for downstream analysis
motifs <- c("Klf6","WT1","Klf15","Tcf21")
markerMotifs <- getFeatures(Klf6_Proj5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

Klf6_Proj5 <- addImputeWeights(Klf6_Proj5)

#plot the distribution of chromVAR deviation scores for each cluster
p <- plotGroups(ArchRProj = Klf6_Proj5, 
                groupBy = "predictedGroup", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(Klf6_Proj5)
)

#We can use cowplot to plot the distributions of all of these motifs in a single plot.
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

#Instead of looking at the distributions of these z-scores, we can overlay the z-scores on our UMAP embedding as we’ve done previously for gene scores.
p <- plotEmbedding(
  ArchRProj = Klf6_Proj5, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Klf6_Proj5)
)

#We can plot all of these motif UMAPs using cowplot.
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#To see how these TF deviation z-scores compare to the inferred gene expression via gene scores of the corresponding TF genes, we can overlay the gene scores for each of these TFs on the UMAP embedding.
markerRNA <- getFeatures(Klf6_Proj5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA

p <- plotEmbedding(
  ArchRProj = Klf6_Proj5, 
  colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Klf6_Proj5)
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#because we previously linked our scATAC-seq data with corresponding scRNA-seq data, we can plot the linked gene expression for each of these TFs on the UMAP embedding.
markerRNA <- getFeatures(Klf6_Proj5, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA

p <- plotEmbedding(
  ArchRProj = Klf6_Proj5, 
  colorBy = "GeneIntegrationMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  continuousSet = "blueYellow",
  imputeWeights = getImputeWeights(Klf6_Proj5)
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#Motif Footprinting
motifPositions <- getPositions(Klf6_Proj5)

motifPositions

motifs <- c("Sfpi1", "Smarcc1", "Fos", "Bcl11a")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

#get footprints
seFoot <- getFootprints(
  ArchRProj = Klf6_Proj5, 
  useGroups = c("InjuredPT-A"),
  positions = motifPositions[markerMotifs], 
  groupBy = "Sample"
)

#Normalization of Footprints for Tn5 Bias
plotFootprints(
  seFoot = seFoot,
  ArchRProj = Klf6_Proj5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

#Dividing by the Tn5 Bias
plotFootprints(
  seFoot = seFoot,
  ArchRProj = Klf6_Proj5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

#Feature Footprinting
seTSS <- getFootprints(
  ArchRProj = Klf6_Proj5, 
  positions = GRangesList(TSS = getTSS(Klf6_Proj5)), 
  groupBy = "predictedGroup",
  flank = 2000
)

#plot footprints
plotFootprints(
  seFoot = seTSS,
  ArchRProj = Klf6_Proj5, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100
)

#Co-accessibility with ArchR
Klf6_Proj5 <- addCoAccessibility(
  ArchRProj = Klf6_Proj5,
  reducedDims = "IterativeLSI"
)

#retrieve Co-accessibility
cA <- getCoAccessibility(
  ArchRProj = Klf6_Proj5,
  corCutOff = 0.5,
  resolution = 50000,
  returnLoops = TRUE
)

cA[[1]]

#Plotting browser tracks of Co-accessibility
markerGenes  <- c(
  "Clu"
)

p <- plotBrowserTrack(
  ArchRProj = Klf6_Proj5, 
  groupBy = "predictedGroup", 
  geneSymbol = markerGenes, 
  upstream = 55000,
  downstream = 55000,
  features = getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["PEC"],
  loops = getCoAccessibility(Klf6_Proj5)
)

grid::grid.newpage()
grid::grid.draw(p$Cd44)

#Peak2GeneLinkage with ArchR
Klf6_Proj5 <- addPeak2GeneLinks(
  ArchRProj = Klf6_Proj5,
  cellsToUse = "PEC",
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = Klf6_Proj5,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g
metadata(p2g)[[1]]

#set returnLoops = TRUE, then getPeak2GeneLinks() will return a loop track GRanges object that connects the peak and gene
p2g <- getPeak2GeneLinks(
  ArchRProj = Klf6_Proj5,
  corCutOff = 0.45,
  resolution = 10000,
  returnLoops = TRUE
)

p2g[[1]]

#Plotting browser tracks with peak-to-gene links

p <- plotBrowserTrack(
  ArchRProj = Klf6_Proj5, 
  groupBy = "predictedGroup", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  features = getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Pod"],
  loops = getPeak2GeneLinks(Klf6_Proj5)
)

grid::grid.newpage()
grid::grid.draw(p$Clu)


