#Script for Analysis of 4 groups  with 2 biological replicates each
library(dplyr)
library(Seurat)
library(Matrix)
library(cowplot)


# Set up control objects
ctrl1 <- Read10X(data.dir = 'SC8_2/')
ctrl1 <- CreateSeuratObject(counts = ctrl1, min.cells = 3, min.features = 200, project = "NPHS2_rtTA_1")
ctrl1$OE <- "NPHS2_rtTA"
ctrl2 <- Read10X(data.dir = 'SC15_4/')
ctrl2 <- CreateSeuratObject(counts = ctrl2, min.cells = 3, min.features = 200, project = "NPHS2_rtTA_2")
ctrl2$OE <- "NPHS2_rtTA"



# Set up OE objects
OE1 <- Read10X(data.dir = 'SC8_3/')
OE1 <- CreateSeuratObject(counts = OE1, min.cells = 3, min.features = 200, project = "hKLF6PODTA_1")
OE1$OE <- "hKLF6PODTA"
OE2 <- Read10X(data.dir = 'SC15_2/')
OE2 <- CreateSeuratObject(counts = OE2, min.cells = 3, min.features = 200, project = "hKLF6PODTA_2")
OE2$OE <- "hKLF6PODTA"

# Set up WTSHAM objects
WTSHAM1 <- Read10X(data.dir = 'SC8_4/')
WTSHAM1 <- CreateSeuratObject(counts = WTSHAM1, min.cells = 3, min.features = 200, project = "SHAM_NPHS2_rtTA_1")
WTSHAM1$OE <- "SHAM_NPHS2_rtTA"
WTSHAM2 <- Read10X(data.dir = 'SC15_3/')
WTSHAM2 <- CreateSeuratObject(counts = WTSHAM2, min.cells = 3, min.features = 200, project = "SHAM_NPHS2_rtTA_2")
WTSHAM2$OE <- "SHAM_NPHS2_rtTA"


# Set up OESHAM objects
OESHAM1 <- Read10X(data.dir = 'SC8_1/')
OESHAM1 <- CreateSeuratObject(counts = OESHAM1, min.cells = 3, min.features = 200, project = "SHAM_hKLF6PODTA_1")
OESHAM1$OE <- "SHAM_hKLF6PODTA"
OESHAM2 <- Read10X(data.dir = 'SC15_1/')
OESHAM2 <- CreateSeuratObject(counts = OESHAM2, min.cells = 3, min.features = 200, project = "SHAM_hKLF6PODTA_2")
OESHAM2$OE <- "SHAM_hKLF6PODTA"



#Merge Controls into a seurat object
WTCombined <- merge(ctrl1, y = c(ctrl2), add.cell.ids = c("NPHS2_rtTA_1", "NPHS2_rtTA_2"), project = "Combined")
WT <- WTCombined

#Merge Klf6 OEs into a seurat object
OECombined <- merge(OE1, y = c(OE2), add.cell.ids = c("hKLF6PODTA_1", "hKLF6PODTA_2"), project = "Combined")
OE <- OECombined

#Merge sham Controls into a seurat object
WTSHAMCombined <- merge(WTSHAM1, y = c(WTSHAM2), add.cell.ids = c("SHAM_NPHS2_rtTA_1", "SHAM_NPHS2_rtTA_2"), project = "Combined")
WTSHAM <- WTSHAMCombined

#Merge sham  Klf6 OEs into a seurat object
OESHAMCombined <- merge(OESHAM1, y = c(OESHAM2), add.cell.ids = c("SHAM_hKLF6PODTA_1", "SHAM_hKLF6PODTA_2"), project = "Combined")
OESHAM <- OESHAMCombined

#Merge Controls, OEs, and SHAMS into a seurat object
DataCombined <- merge(WT, y = c(OE, WTSHAM, OESHAM), add.cell.ids = c("NPHS2_rtTA", "hKLF6PODTA", "SHAM_NPHS2_rtTA", "SHAM_hKLF6PODTA"), project = "Combined")

#Assess mt content
DataCombined[["percent.mt"]] <- PercentageFeatureSet(object = DataCombined, pattern = "^mt-")
VlnPlot(DataCombined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Filter cells, only keep cells with more than 200 and less than 5000 features, and less than 10% MT-genes
DataCombined <- subset(x= DataCombined, subset= nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <5)
#DataCombined <- SCTransform(DataCombined, vars.to.regress = "percent.mt", verbose = FALSE)
DataCombined <- NormalizeData(object = DataCombined, verbose = FALSE)
DataCombined <- FindVariableFeatures(object = DataCombined, selection.method = "vst", nfeatures = 2000)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(DataCombined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DataCombined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Identify Variable Features 
obj.list <- SplitObject(DataCombined, split.by = "orig.ident")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

memory.limit(size = 10000)


# Perform Integration and Identify Anchors
Data.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
Data.Combined <- IntegrateData(anchorset = Data.anchors, dims = 1:30)


# Run the standard workflow for visualization and clustering
Data.Combined <- ScaleData(object = Data.Combined, verbose = FALSE)
Data.Combined <- RunPCA(object = Data.Combined, npcs = 30, verbose = FALSE) 

# t-SNE and Clustering
Data.Combined <- RunUMAP(object = Data.Combined, reduction = "pca", dims = 1:30) 
Data.Combined <- FindNeighbors(object = Data.Combined, reduction = "pca", dims = 1:30)
Data.Combined <- FindClusters(Data.Combined, resolution = 0.5)

Idents(object = Data.Combined) <- "integrated_snn_res.1"

# Visualization of UMAP
DimPlot(Data.Combined, reduction = "umap", group.by = "orig.ident")
DimPlot(Data.Combined, reduction = "umap", group.by = "OE")
DimPlot(Data.Combined, pt.size = 0, reduction = "umap", label = TRUE) + NoLegend()

#show WT and OE side by side
DimPlot(Data.Combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = FALSE)
DimPlot(Data.Combined, reduction = "umap", split.by = "OE", label = TRUE)

###############
## Save File ##
###############
saveRDS(Data.Combined, file = "Data.Combined.6.15.21.rds")
Data.Combined <- readRDS("Data.CombinedAC.rds")


saveRDS(Data.anchors, file = "Dataanchors1.rds")
Data.Combined.c3 <- readRDS("c0.ASC1.rds")

#####################
## Load Saved File ##
#####################
Data.Combined <- readRDS("AllSamples_subclustered.rds")
Data.Combined@meta.data

#To output average expression split by group
cluster.averages <- AverageExpression(Data.Combined, assays = "RNA", add.ident = "orig.ident")
write.csv(cluster.averages[["RNA"]], file = "Allsamples_ClusterAverages.orig.ident.csv")

# Get number of cells per cluster and per sample/group of origin
table(Data.Combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$orig.ident)
table(Data.Combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$OE)
table(Data.Combined@meta.data[["integrated_snn_res.1"]], Data.Combined@meta.data$OE)
table(Data.Combined@meta.data[["celltype"]], Data.Combined@meta.data$OE)


head(Data.Combined@meta.data)

#Feature and Violin Plots
DefaultAssay(object = Data.Combined) <- "RNA"
VlnPlot(Data.Combined, features = c("Nphs1", "Nphs2", "Synpo", "Wt1"), split.by = "OE", ncol = 1, pt.size = 0)
VlnPlot(Data.Combined, features = c("Akap12"), split.by = "OE", pt.size = 1)
FeaturePlot(Data.Combined, features = c("Nphs1"), ncol = 2, split.by = "orig.ident")
FeaturePlot(Data.Combined, features = c("Camk1d"), ncol = 2, split.by = "OE", pt.size = 0.01)
FeaturePlot(Data.Combined, features = c("Camk1d"), pt.size = 0.1)

DotPlot(Data.Combined2, features = c("Myh9"), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) + RotatedAxis()

#For Dot plot (combined)
markers.to.plot <- c("Synpo", "Ptprm", "Csmd1", "Cd44", "Pcdh7", "Slc7a8", "Slco1a6", "Cyp7b1", "Slc12a1", "Slc12a3", "Kl", "Egfem1", "Frmpd4", "Adgrf5", "Lsamp", "Fyb", "Pigr", "Klk1")
pdf("DotPlot.pdf", width=8, height=6)
DotPlot(Data.Combined, features = rev(markers.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) + RotatedAxis()
dev.off()

markers.to.plot <- c("Tmem27",  "Slc5a12", "Cntnap5a", "Keg1",  "Cyp7b1", "Slc12a1",  "Slc12a3", "Trpm6", "Egfem1",  "Frmpd4", "Clnk", "Insrr", "Ptprc","Top2a", "Cd44", "Cfh", "Flt1","Nphs1")
DotPlot(Data.Combined, features = (markers.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) + RotatedAxis()

#Assigning new cluster names
new.cluster.ids <-           c("Novel_PT", "M?-2", "M?-1", "IC-B", "IC-A", "CD-PC", "CNT", "DCT/CNT", "DCT", "LH(AL)", "PT(S3)/LH(DL)", "PT(S3)", "PT(S1-S2)", "PT(S1-S3)", "PEC/Prolif.PT", "Mes", "Endo", "Pod")
names(new.cluster.ids) <- levels(Data.Combined)
Data.Combined <- RenameIdents(Data.Combined, new.cluster.ids)
DimPlot(Data.Combined, reduction = "umap", label = TRUE, pt.size = 0.5)

#new.cluster.ids <- c("unknown1", "PT(S1-S3)", "PT(S3)", "PT(S1-S2)", "PT(S3)", "PT(S3)/LH(DL)", "LH(AL)", "DCT", "Endo", "CNT", "PEC/Prolif.PT", "CD-PC", "Mes", "IC-A", "M??", "IC-B", "DCT/CNT", "Pod", "M??")
new.cluster.ids <-       c("unknown1", "M?-2", "M?", "IC-B", "IC-A", "CD-PC", "CNT", "DCT/CNT", "DCT", "LH(AL)", "PT(S3)/LH(DL)", "PT(S3)-2", "PT(S3)", "PT(S1-S2)", "PT(S1-S3)", "PEC/Prolif.PT", "Mes", "Endo", "Pod")
new.cluster.ids <-       c("Novel_PT", "M?-2", "M?-1", "IC-B", "IC-A", "CD-PC", "CNT", "DCT/CNT", "DCT", "LH(AL)", "PT(S3)/LH(DL)", "PT(S3)", "PT(S1-S2)", "PT(S1-S3)", "PEC/Prolif.PT", "Mes", "Endo", "Pod")
new.cluster.ids <-       c("Novel_PT", "PT(S1-S2)", "PT(S1-S3)", "PT(S3)", "PT(S3)/LH(DL)", "LH(AL)", "DCT", "DCT/CNT", "CNT", "CD-PC", "IC-A", "IC-B", "M?-1", "M?-2", "PEC/Prolif.PT", "Mes", "Endo", "Pod")
new.cluster.ids <-       c("M?-1", "M?-2")



#Reorder new cluster names
levels(x = Data.Combined) <- c("Novel_PT", "PT(S1-S2)", "PT(S1-S3)", "PT(S3)", "PT(S3)/LH(DL)", "LH(AL)", "DCT", "DCT/CNT", "CNT", "CD-PC", "IC-A", "IC-B", "M?-1", "M?-2", "PEC/Prolif.PT", "Mes", "Endo", "Pod")

levels(x = Data.Combined) <- c("Novel_PT_SHAM_NPHS2_rtTA", "Novel_PT_SHAM_hKLF6PODTA", "Novel_PT_NPHS2_rtTA", "Novel_PT_hKLF6PODTA", "PT(S3)_SHAM_NPHS2_rtTA", "PT(S3)_SHAM_hKLF6PODTA", "PT(S3)_NPHS2_rtTA", "PT(S3)_hKLF6PODTA", "PT(S1-S2)_SHAM_NPHS2_rtTA", "PT(S1-S2)_SHAM_hKLF6PODTA", "PT(S1-S2)_NPHS2_rtTA", "PT(S1-S2)_hKLF6PODTA", 
                               "PT(S1-S3)_SHAM_NPHS2_rtTA", "PT(S1-S3)_SHAM_hKLF6PODTA", "PT(S1-S3)_NPHS2_rtTA", "PT(S1-S3)_hKLF6PODTA", "PT(S3)/LH(DL)_SHAM_NPHS2_rtTA", "PT(S3)/LH(DL)_SHAM_hKLF6PODTA", "PT(S3)/LH(DL)_NPHS2_rtTA", "PT(S3)/LH(DL)_hKLF6PODTA", "LH(AL)_SHAM_NPHS2_rtTA", "LH(AL)_SHAM_hKLF6PODTA", "LH(AL)_NPHS2_rtTA",
                               "LH(AL)_hKLF6PODTA", "DCT_SHAM_NPHS2_rtTA", "DCT_SHAM_hKLF6PODTA", "DCT_NPHS2_rtTA", "DCT_hKLF6PODTA", "DCT/CNT_SHAM_NPHS2_rtTA", "DCT/CNT_SHAM_hKLF6PODTA", "DCT/CNT_NPHS2_rtTA", "DCT/CNT_hKLF6PODTA", "CNT_SHAM_NPHS2_rtTA", "CNT_SHAM_hKLF6PODTA", "CNT_NPHS2_rtTA", "CNT_hKLF6PODTA", 
                               "CD-PC_SHAM_NPHS2_rtTA", "CD-PC_SHAM_hKLF6PODTA", "CD-PC_NPHS2_rtTA", "CD-PC_hKLF6PODTA", "IC-A_SHAM_NPHS2_rtTA", "IC-A_SHAM_hKLF6PODTA", "IC-A_NPHS2_rtTA", "IC-A_hKLF6PODTA", "IC-B_SHAM_NPHS2_rtTA", "IC-B_SHAM_hKLF6PODTA", "IC-B_NPHS2_rtTA", "IC-B_hKLF6PODTA", 
                               "M?-1_SHAM_NPHS2_rtTA", "M?-1_SHAM_hKLF6PODTA", "M?-1_NPHS2_rtTA", "M?-1_hKLF6PODTA", "M?-2_SHAM_NPHS2_rtTA", "M?-2_SHAM_hKLF6PODTA", "M?-2_NPHS2_rtTA", "M?-2_hKLF6PODTA", "PEC/Prolif.PT_SHAM_NPHS2_rtTA", "PEC/Prolif.PT_SHAM_hKLF6PODTA", "PEC/Prolif.PT_NPHS2_rtTA", "PEC/Prolif.PT_hKLF6PODTA",
                               "Mes_SHAM_NPHS2_rtTA", "Mes_SHAM_hKLF6PODTA", "Mes_NPHS2_rtTA", "Mes_hKLF6PODTA", "Endo_SHAM_NPHS2_rtTA", "Endo_SHAM_hKLF6PODTA", "Endo_NPHS2_rtTA", "Endo_hKLF6PODTA", "Pod_SHAM_NPHS2_rtTA", "Pod_SHAM_hKLF6PODTA", "Pod_NPHS2_rtTA", "Pod_hKLF6PODTA")

new.cluster.ids <-           c("Novel_PT", "M?-2", "M?-1", "IC-B", "IC-A", "CD-PC", "CNT", "DCT/CNT", "DCT", "LH(AL)", "PT(S3)/LH(DL)", "PT(S3)", "PT(S1-S2)", "PT(S1-S3)", "PEC/Prolif.PT", "Mes", "Endo", "Pod")
#use this order for everything that isn't Dot Plot
levels(x = Data.Combined) <- rev(c("unknown1", "Macrophages", "Macrophages2", "Macrophages3", "IC-B", "IC-A", "CD-PC", "CnT2", "CnT", "DCT","LoH-DL", "LoH-AL", "LoH-AL2", "PT (S3)2", "PT (S3)", "PT (S1-S2)2", "PT (S1-S2)", "PT (S1-S2/S3)", "MC", "EC", "Podocytes"))

#Receptors on Dot Plot
DefaultAssay(Data.Combined) <- "RNA"
receptors.to.plot <- c("Cd44", "Cd74", "Traf2", "Bsg", "Itga5", "Itgam", "Itgb2", "Ddr1", "Itga1", "Itga10", "Itga2", "Dip2a", "Bmp4", "Pappa", "Fzd8", "Igf1R", "Acvr1", "Eng", "Tgfbr1", "Tgfbr2", "Tgfbr3", "Cd47", "Cd93", "Itgav", "Itgb8", "Itgb1", "Fzd1", "Acvrl1", "Cav1", "Cd109", "Cxcr4", "Itgb6", "Sdc2", "Itgb3", "Acvr1b", "Flt4", "Il17cr", "Itga6", "Itga9", "Itgb7", "Plaur", "Ddr2", "Lrp4", "Sdc1")
DotPlot(Data.Combined, features = rev(receptors.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8, split.by = "OE") + RotatedAxis()
p = DotPlot(Data.Combined, features = rev(receptors.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.y = element_text(angle = 150, hjust = 0, vjust=0.5))

levels(x = immune.combined) <- c( "0.CTRL", "0.KO", "1.CTRL", "1.KO", "2.CTRL", "2.KO", "3.CTRL", "3.KO", "4.CTRL", "4.KO", "5.CTRL", "5.KO", "6.CTRL", "6.KO", "7.CTRL", "7.KO", "8.CTRL", "8.KO")
DefaultAssay(immune.combined) <- "RNA"
receptors.to.plot <- c("Lrp4", "Cd44", "Cd93", "Ddr1", "Ddr2", "Flt4", "Itga1", "Itga2", "Itga5", "Itgav", "Itga10", "Cd47", "Itgb8", "Itgam", "Itgb2", "Itgb3", "Itgb6", "Il17rc", "Itga6", "Itga9", "Itgb7", "Plaur", "Traf2", "Fzd8", "Acvr1b", "Tmem222", "SDC1", "Cd74")
p = DotPlot(immune.combined, features = rev(receptors.to.plot), cols = c("lightgray", "limegreen"), col.min = 0, dot.scale = 8) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.y = element_text(angle = 150, hjust = 0, vjust=0.5))


DefaultAssay(Data.Combined) <- "RNA"
receptors.to.plot <- c("Lrp4", "Cd44", "Cd93", "Ddr1", "Ddr2", "Flt4", "Itga1", "Itga2", "Itga5", "Itgav", "Itga10", "Cd47", "Itgb8", "Itgam", "Itgb2", "Itgb3", "Itgb6", "Il17rc", "Itga6", "Itga9", "Itgb7", "Plaur", "Traf2", "Fzd8", "Acvr1b", "Tmem222", "SDC1", "Cd74")
p = DotPlot(Data.Combined, features = rev(receptors.to.plot), cols = c("lightgray", "limegreen"), col.min = 0, dot.scale = 8) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + theme(axis.text.y = element_text(angle = 150, hjust = 0, vjust=0.5))

#Ligands on Dot Plot
DefaultAssay(Data.Combined) <- "RNA"
ligands.to.plot <- c("Agrn", "Anxa5", "Bmp1", "Cnn2", "Col1a1", "Col2a1", "Col4a1", "Col4a2", "Col5a2", "Cpe", "Ctgf", "Fbn1", "Fn1", "Fstl1", "Gfra1", "Gstp1", "Hapln1", "Igfbp4", "Igfbp5", "Igfbp7", "Inhba", "Lgals1", "Loxl2", "Ltbp2", "Mif", "Nucb1", "P4hb", "Ppia", "Prss23", "Sparc", "Stc2", "Tgfb2", "Tgfbi", "Thbs1", "Tkt")
DotPlot(Data.Combined, features = rev(ligands.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8, split.by = "KO") + RotatedAxis()

FeatureScatter(Data.Combined, feature1 = "Cd44", feature2 = "Itgav")

#Another way to name clusters
new.cluster.ids <- c("unknown1", "PT(S1-S3)", "PT(S3)", "PT(S1-S2)", "PT(S1-S3)/LoH(DL)", "PT(S3)/LoH(DL)", "PT(S1-S2)", "LoH(AL)", "DCT", "EC", "CnT", "PEC/Prolif.PT", "EC", "CD-PC", "MC", "IC-A", "Macrophages", "LoH(AL)", "IC-B", "DCT/CnT", "PT(S3)", "Podocytes", "Inflammatory cells")
names(new.cluster.ids) <- levels(Data.Combined)
Data.Combined <- RenameIdents(Data.Combined, new.cluster.ids)
DimPlot(Data.Combined, reduction = "umap", label = TRUE, pt.size = 0) + NoLegend()

c22m <- FindConservedMarkers(Data.Combined, ident.1 = "c22", grouping.var = "OE", verbose = FALSE)
write.csv(c22m, file = "c22markers.csv")
new.cluster.ids <- c("unknown1", "PT(S1-S3)", "PT(S3)", "PT(S1-S2)", "PT/LH(DL)", "PT(S3)/LH(DL)", "PT(S1-S2)-2", "LH(AL)", "DCT", "Endo", "CNT", "PEC/Prolif.PT", "Endo-2", "CD-PC", "Mes", "IC-A", "M??", "LH(AL)-2", "IC-B", "DCT/CNT", "PT(S3)-2", "Pod", "M??-2")
new.cluster.ids <- c("unknown1", "PT(S1-S3)", "PT(S3)", "PT(S1-S2)", "PT/LH(DL)", "PT(S3)/LH(DL)", "PT(S1-S2)", "LH(AL)", "DCT", "Endo", "CNT", "PEC/Prolif.PT", "Endo", "CD-PC", "Mes", "IC-A", "M??", "LH(AL)", "IC-B", "DCT/CNT", "PT(S3)", "Pod", "M??")
new.cluster.ids <- c("unknown1", "PT(S1-S3)", "PT(S3)", "PT(S1-S2)", "PT/LH(DL)", "PT(S3)/LH(DL)", "LH(AL)", "DCT", "Endo", "CNT", "PEC/Prolif.PT", "CD-PC", "Mes", "IC-A", "M??", "IC-B", "DCT/CNT", "Pod", "M??")


options(future.globals.maxSize= 159128960000)
new.cluster.idss <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")
new.cluster.idss <- c("12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")
new.cluster.idss <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")

for(i in new.cluster.ids) {
  clus.markers <- FindMarkers(Data.Combined, ident.1 = i, ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
  head(clus.markers)
  str1 = i
  str2 = "deg.csv"
  print(i)
  result = paste(str1,str2, sep = ".")
  write.csv(clus.markers, file = result)
}
C0 <- read.csv("0.deg.csv")
C1 <- read.csv("1.deg.csv")
C2 <- read.csv("2.deg.csv")
C3 <- read.csv("3.deg.csv")
C4 <- read.csv("4.deg.csv")
C5 <- read.csv("5.deg.csv")
C6 <- read.csv("6.deg.csv")
C7 <- read.csv("7.deg.csv")
C8 <- read.csv("8.deg.csv")
C9 <- read.csv("9.deg.csv")
C10 <- read.csv("10.deg.csv")
C11 <- read.csv("11.deg.csv")
C12 <- read.csv("12.deg.csv")
C13 <- read.csv("13.deg.csv")
C14 <- read.csv("14.deg.csv")
C15 <- read.csv("15.deg.csv")
C16 <- read.csv("16.deg.csv")
C17 <- read.csv("17.deg.csv")
C18 <- read.csv("18.deg.csv")
C19 <- read.csv("19.deg.csv")
C20 <- read.csv("20.deg.csv")
C21 <- read.csv("21.deg.csv")
C22 <- read.csv("22.deg.csv")
#C23 <- read.csv("23.deg.csv")
#C24 <- read.csv("24.deg.csv")


C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20, "C21" = C21, "C22" = C22)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18)
write.xlsx(list_of_datasets, file = "Summary.DEGforclustermarkers.xlsx")

for(each in new.cluster.idss ) {
  t.cells <- subset(Data.Combined, idents = each)
  Idents(t.cells) <- "OE"
  avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
  avg.t.cells$gene <- rownames(avg.t.cells)
  str1 = each
  str2 = "average.exp.csv"
  print(each)
  result = paste(str1,str2, sep = "")
  write.csv(avg.t.cells, file = result)
}

c0a <- read.csv("0average.exp.csv")
c1a <- read.csv("1average.exp.csv")
c2a <- read.csv("2average.exp.csv")
c3a <- read.csv("3average.exp.csv")
c4a <- read.csv("4average.exp.csv")
c5a <- read.csv("5average.exp.csv")
c6a <- read.csv("6average.exp.csv")
c7a <- read.csv("7average.exp.csv")
c8a <- read.csv("8average.exp.csv")
c9a <- read.csv("9average.exp.csv")
c10a <- read.csv("10average.exp.csv")
c11a <- read.csv("11average.exp.csv")
c12a <- read.csv("12average.exp.csv")
c13a <- read.csv("13average.exp.csv")
c14a <- read.csv("14average.exp.csv")
c15a <- read.csv("15average.exp.csv")
c16a <- read.csv("16average.exp.csv")
c17a <- read.csv("17average.exp.csv")
c18a <- read.csv("18average.exp.csv")
c19a <- read.csv("19average.exp.csv")
c20a <- read.csv("20average.exp.csv")
#c21a <- read.csv("c21average.exp.csv")
#c22a <- read.csv("c22average.exp.csv")
#c23a <- read.csv("c23average.exp.csv")
#c24a <- read.csv("c24average.exp.csv")

c0a <- data.frame(Gene = row.names(c0a), c0a)
c1a <- data.frame(Gene = row.names(c1a), c1a)
c2a <- data.frame(Gene = row.names(c2a), c2a)
c3a <- data.frame(Gene = row.names(c3a), c3a)
c4a <- data.frame(Gene = row.names(c4a), c4a)
c5a <- data.frame(Gene = row.names(c5a), c5a)
c6a <- data.frame(Gene = row.names(c6a), c6a)
c7a <- data.frame(Gene = row.names(c7a), c7a)
c8a <- data.frame(Gene = row.names(c8a), c8a)
c9a <- data.frame(Gene = row.names(c9a), c9a)
c10a <- data.frame(Gene = row.names(c10a), c10a)
c11a <- data.frame(Gene = row.names(c11a), c11a)
c12a <- data.frame(Gene = row.names(c12a), c12a)
c13a <- data.frame(Gene = row.names(c13a), c13a)
c14a <- data.frame(Gene = row.names(c14a), c14a)
c15a <- data.frame(Gene = row.names(c15a), c15a)
c16a <- data.frame(Gene = row.names(c16a), c16a)
c17a <- data.frame(Gene = row.names(c17a), c17a)
c18a <- data.frame(Gene = row.names(c18a), c18a)
c19a <- data.frame(Gene = row.names(c19a), c19a)
c20a <- data.frame(Gene = row.names(c20a), c20a)
#c21a <- data.frame(Gene = row.names(c21a), c21a)
#c22a <- data.frame(Gene = row.names(c22a), c22a)
#c23a <- data.frame(Gene = row.names(c23a), c23a)
#c24a <- data.frame(Gene = row.names(c24a), c24a)

list_of_datasets <- list("c0a" = c0a, "c1a" = c1a, "c2a" = c2a, "c3a" = c3a, "c4a" = c4a, "c5a" = c5a, "c6a" = c6a, "c7a" = c7a, "c8a" = c8a, "c9a" = c9a, "c10a" = c10a, "c11a" = c11a, "c12a" = c12a, "c13a" = c13a, "c14a" = c14a, "c15a" = c15a, "c16a" = c16a, "c17a" = c17a, "c18a" = c18a, "c19a" = c19a, "c20a" = c20a)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.Average.Exp.Genes.xlsx")

immune.combined <- Data.Combined
immune.combined$celltype.OE <- paste(Idents(immune.combined), immune.combined$OE, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.OE"
head(immune.combined@meta.data)

MC.response <- FindMarkers(immune.combined, ident.1 = "Endo_NPHS2_rtTA", ident.2 = "Endo_SHAM_NPHS2_rtTA", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(MC.response, n = 10)
write.csv(MC.response, file = "DEG_Endo_combined.WTvWTSHAM.csv")

for(i in new.cluster.idss) {
  a = paste(i, "hKLF6PODTA", sep = "_")
  print(a)
  b = paste(i, "NPHS2_rtTA", sep = "_")
  c = paste(i, "SHAM_hKLF6PODTA", sep = "_")
  d = paste(i, "SHAM_NPHS2_rtTA", sep = "_")
  MC.response <- FindMarkers(immune.combined, ident.1 = b, ident.2 = d, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}

C0 <- read.csv("0comparison.csv")
C1 <- read.csv("1comparison.csv")
C2 <- read.csv("2comparison.csv")
C3 <- read.csv("3comparison.csv")
C4 <- read.csv("4comparison.csv")
C5 <- read.csv("5comparison.csv")
C6 <- read.csv("6comparison.csv")
C7 <- read.csv("7comparison.csv")
C8 <- read.csv("8comparison.csv")
C9 <- read.csv("9comparison.csv")
C10 <- read.csv("10comparison.csv")
C11 <- read.csv("11comparison.csv")
C12 <- read.csv("12comparison.csv")
C13 <- read.csv("13comparison.csv")
C14 <- read.csv("14comparison.csv")
C15 <- read.csv("15comparison.csv")
C16 <- read.csv("16comparison.csv")
C17 <- read.csv("17comparison.csv")
C18 <- read.csv("18comparison.csv")
C19 <- read.csv("19comparison.csv")
C20 <- read.csv("20comparison.csv")
C21 <- read.csv("21comparison.csv")
C22 <- read.csv("22comparison.csv")
#C23 <- read.csv("c23comparison.csv")
#C24 <- read.csv("c24comparison.csv")

C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20, "C21" = C21, "C22" = C22)
write.xlsx(list_of_datasets, file = "Summarymt5.0.DEG_WT_v_WTSHAM.xlsx")



for(i in new.cluster.idss) {
  a = paste(i, "OE", sep = "_")
  print(a)
  b = paste(i, "WT", sep = "_")
  c = paste(i, "OE.SHAM", sep = "_")
  d = paste(i, "WT.SHAM", sep = "_")
  MC.response <- FindMarkers(immune.combined, ident.1 = c, ident.2 = d, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}

C0 <- read.csv("c0comparison.csv")
C1 <- read.csv("c1comparison.csv")
C2 <- read.csv("c2comparison.csv")
C3 <- read.csv("c3comparison.csv")
C4 <- read.csv("c4comparison.csv")
C5 <- read.csv("c5comparison.csv")
C6 <- read.csv("c6comparison.csv")
C7 <- read.csv("c7comparison.csv")
C8 <- read.csv("c8comparison.csv")
C9 <- read.csv("c9comparison.csv")
C10 <- read.csv("c10comparison.csv")
C11 <- read.csv("c11comparison.csv")
C12 <- read.csv("c12comparison.csv")
C13 <- read.csv("c13comparison.csv")
C14 <- read.csv("c14comparison.csv")
C15 <- read.csv("c15comparison.csv")
C16 <- read.csv("c16comparison.csv")
C17 <- read.csv("c17comparison.csv")
C18 <- read.csv("c18comparison.csv")
C19 <- read.csv("c19comparison.csv")
C20 <- read.csv("c20comparison.csv")
#C21 <- read.csv("c21comparison.csv")
#C22 <- read.csv("c22comparison.csv")
#C23 <- read.csv("c23comparison.csv")
#C24 <- read.csv("c24comparison.csv")

C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
#C21 <- data.frame(Gene = row.names(C21), C21)
#C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.DEG_OESHAM_vs_WTSHAM.xlsx")
for(i in new.cluster.idss) {
  a = paste(i, "OE", sep = "_")
  print(a)
  b = paste(i, "WT", sep = "_")
  c = paste(i, "OE.SHAM", sep = "_")
  d = paste(i, "WT.SHAM", sep = "_")
  MC.response <- FindMarkers(immune.combined, ident.1 = b, ident.2 = d, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}

C0 <- read.csv("c0comparison.csv")
C1 <- read.csv("c1comparison.csv")
C2 <- read.csv("c2comparison.csv")
C3 <- read.csv("c3comparison.csv")
C4 <- read.csv("c4comparison.csv")
C5 <- read.csv("c5comparison.csv")
C6 <- read.csv("c6comparison.csv")
C7 <- read.csv("c7comparison.csv")
C8 <- read.csv("c8comparison.csv")
C9 <- read.csv("c9comparison.csv")
C10 <- read.csv("c10comparison.csv")
C11 <- read.csv("c11comparison.csv")
C12 <- read.csv("c12comparison.csv")
C13 <- read.csv("c13comparison.csv")
C14 <- read.csv("c14comparison.csv")
C15 <- read.csv("c15comparison.csv")
C16 <- read.csv("c16comparison.csv")
C17 <- read.csv("c17comparison.csv")
C18 <- read.csv("c18comparison.csv")
C19 <- read.csv("c19comparison.csv")
C20 <- read.csv("c20comparison.csv")
#C21 <- read.csv("c21comparison.csv")
#C22 <- read.csv("c22comparison.csv")
#C23 <- read.csv("c23comparison.csv")
#C24 <- read.csv("c24comparison.csv")

C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
#C21 <- data.frame(Gene = row.names(C21), C21)
#C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.DEG WTvsWTSHAM.Genes.xlsx")
for(i in new.cluster.idss) {
  a = paste(i, "OE", sep = "_")
  print(a)
  b = paste(i, "WT", sep = "_")
  c = paste(i, "OE.SHAM", sep = "_")
  d = paste(i, "WT.SHAM", sep = "_")
  MC.response <- FindMarkers(immune.combined, ident.1 = c, ident.2 = d, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}

C0 <- read.csv("c0comparison.csv")
C1 <- read.csv("c1comparison.csv")
C2 <- read.csv("c2comparison.csv")
C3 <- read.csv("c3comparison.csv")
C4 <- read.csv("c4comparison.csv")
C5 <- read.csv("c5comparison.csv")
C6 <- read.csv("c6comparison.csv")
C7 <- read.csv("c7comparison.csv")
C8 <- read.csv("c8comparison.csv")
C9 <- read.csv("c9comparison.csv")
C10 <- read.csv("c10comparison.csv")
C11 <- read.csv("c11comparison.csv")
C12 <- read.csv("c12comparison.csv")
C13 <- read.csv("c13comparison.csv")
C14 <- read.csv("c14comparison.csv")
C15 <- read.csv("c15comparison.csv")
C16 <- read.csv("c16comparison.csv")
C17 <- read.csv("c17comparison.csv")
C18 <- read.csv("c18comparison.csv")
C19 <- read.csv("c19comparison.csv")
C20 <- read.csv("c20comparison.csv")
#C21 <- read.csv("c21comparison.csv")
#C22 <- read.csv("c22comparison.csv")
#C23 <- read.csv("c23comparison.csv")
#C24 <- read.csv("c24comparison.csv")

C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
#C21 <- data.frame(Gene = row.names(C21), C21)
#C22 <- data.frame(Gene = row.names(C22), C22)
#C23 <- data.frame(Gene = row.names(C23), C23)
#C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.DEG OESHAMvsWTSHAM.Genes.xlsx")





new.cluster.idss <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22", "c23", "c24")
#Finding conserved markers using a loop. Note it does not save in environment.
for(i in new.cluster.idss) {
  clus.markers <- FindConservedMarkers(Data.Combined, ident.1 = i, grouping.var = "OE", verbose = FALSE)
  head(clus.markers)
  str1 = i
  str2 = "markers.csv"
  print(i)
  result = paste(str1,str2, sep = "")
  write.csv(clus.markers, file = result)
}

c0m <- read.csv("c0markers.csv")
c1m <- read.csv("c1markers.csv")
c2m <- read.csv("c2markers.csv")
c3m <- read.csv("c3markers.csv")
c4m <- read.csv("c4markers.csv")
c5m <- read.csv("c5markers.csv")
c6m <- read.csv("c6markers.csv")
c7m <- read.csv("c7markers.csv")
c8m <- read.csv("c8markers.csv")
c9m <- read.csv("c9markers.csv")
c10m <- read.csv("c10markers.csv")
c11m <- read.csv("c11markers.csv")
c12m <- read.csv("c12markers.csv")
c13m <- read.csv("c13markers.csv")
c14m <- read.csv("c14markers.csv")
c15m <- read.csv("c15markers.csv")
c16m <- read.csv("c16markers.csv")
c17m <- read.csv("c17markers.csv")
c18m <- read.csv("c18markers.csv")
c19m <- read.csv("c19markers.csv")
c20m <- read.csv("c20markers.csv")
#c21m <- read.csv("c21markers.csv")
c22m <- read.csv("c22markers.csv")
c23m <- read.csv("c23markers.csv")
c24m <- read.csv("c24markers.csv")

c0m <- data.frame(Gene = row.names(c0m), c0m)
c1m <- data.frame(Gene = row.names(c1m), c1m)
c2m <- data.frame(Gene = row.names(c2m), c2m)
c3m <- data.frame(Gene = row.names(c3m), c3m)
c4m <- data.frame(Gene = row.names(c4m), c4m)
c5m <- data.frame(Gene = row.names(c5m), c5m)
c6m <- data.frame(Gene = row.names(c6m), c6m)
c7m <- data.frame(Gene = row.names(c7m), c7m)
c8m <- data.frame(Gene = row.names(c8m), c8m)
c9m <- data.frame(Gene = row.names(c9m), c9m)
c10m <- data.frame(Gene = row.names(c10m), c10m)
c11m <- data.frame(Gene = row.names(c11m), c11m)
c12m <- data.frame(Gene = row.names(c12m), c12m)
c13m <- data.frame(Gene = row.names(c13m), c13m)
c14m <- data.frame(Gene = row.names(c14m), c14m)
c15m <- data.frame(Gene = row.names(c15m), c15m)
c16m <- data.frame(Gene = row.names(c16m), c16m)
c17m <- data.frame(Gene = row.names(c17m), c17m)
c18m <- data.frame(Gene = row.names(c18m), c18m)
c19m <- data.frame(Gene = row.names(c19m), c19m)
c20m <- data.frame(Gene = row.names(c20m), c20m)
#c21m <- data.frame(Gene = row.names(c21m), c21m)
c22m <- data.frame(Gene = row.names(c22m), c22m)
c23m <- data.frame(Gene = row.names(c23m), c23m)
c24m <- data.frame(Gene = row.names(c24m), c24m)


list_of_datasets <- list("c0m" = c0m, "c1m" = c1m, "c2m" = c2m, "c3m" = c3m, "c4m" = c4m, "c5m" = c5m, "c6m" = c6m, "c7m" = c7m, "c8m" = c8m, "c9m" = c9m, "c10m" = c10m, "c11m" = c11m, "c12m" = c12m, "c13m" = c13m, "c14m" = c14m, "c15m" = c15m, "c16m" = c16m, "c17m" = c17m, "c18m" = c18m, "c19m" = c19m, "c20m" = c20m, "c22m" = c22m, "c23m" = c23m, "c24m" = c24m)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.MarkerGenes.xlsx")


new.cluster.idss <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22", "c23", "c24")
for(each in new.cluster.idss ) {
  t.cells <- subset(Data.Combined, idents = each)
  Idents(t.cells) <- "OE"
  avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
  avg.t.cells$gene <- rownames(avg.t.cells)
  str1 = each
  str2 = "average.exp.csv"
  print(each)
  result = paste(str1,str2, sep = "")
  write.csv(avg.t.cells, file = result)
}

c0a <- read.csv("c0average.exp.csv")
c1a <- read.csv("c1average.exp.csv")
c2a <- read.csv("c2average.exp.csv")
c3a <- read.csv("c3average.exp.csv")
c4a <- read.csv("c4average.exp.csv")
c5a <- read.csv("c5average.exp.csv")
c6a <- read.csv("c6average.exp.csv")
c7a <- read.csv("c7average.exp.csv")
c8a <- read.csv("c8average.exp.csv")
c9a <- read.csv("c9average.exp.csv")
c10a <- read.csv("c10average.exp.csv")
c11a <- read.csv("c11average.exp.csv")
c12a <- read.csv("c12average.exp.csv")
c13a <- read.csv("c13average.exp.csv")
c14a <- read.csv("c14average.exp.csv")
c15a <- read.csv("c15average.exp.csv")
c16a <- read.csv("c16average.exp.csv")
c17a <- read.csv("c17average.exp.csv")
c18a <- read.csv("c18average.exp.csv")
c19a <- read.csv("c19average.exp.csv")
c20a <- read.csv("c20average.exp.csv")
c21a <- read.csv("c21average.exp.csv")
c22a <- read.csv("c22average.exp.csv")
c23a <- read.csv("c23average.exp.csv")
c24a <- read.csv("c24average.exp.csv")

c0a <- data.frame(Gene = row.names(c0a), c0a)
c1a <- data.frame(Gene = row.names(c1a), c1a)
c2a <- data.frame(Gene = row.names(c2a), c2a)
c3a <- data.frame(Gene = row.names(c3a), c3a)
c4a <- data.frame(Gene = row.names(c4a), c4a)
c5a <- data.frame(Gene = row.names(c5a), c5a)
c6a <- data.frame(Gene = row.names(c6a), c6a)
c7a <- data.frame(Gene = row.names(c7a), c7a)
c8a <- data.frame(Gene = row.names(c8a), c8a)
c9a <- data.frame(Gene = row.names(c9a), c9a)
c10a <- data.frame(Gene = row.names(c10a), c10a)
c11a <- data.frame(Gene = row.names(c11a), c11a)
c12a <- data.frame(Gene = row.names(c12a), c12a)
c13a <- data.frame(Gene = row.names(c13a), c13a)
c14a <- data.frame(Gene = row.names(c14a), c14a)
c15a <- data.frame(Gene = row.names(c15a), c15a)
c16a <- data.frame(Gene = row.names(c16a), c16a)
c17a <- data.frame(Gene = row.names(c17a), c17a)
c18a <- data.frame(Gene = row.names(c18a), c18a)
c19a <- data.frame(Gene = row.names(c19a), c19a)
c20a <- data.frame(Gene = row.names(c20a), c20a)
c21a <- data.frame(Gene = row.names(c21a), c21a)
c22a <- data.frame(Gene = row.names(c22a), c22a)
c23a <- data.frame(Gene = row.names(c23a), c23a)
c24a <- data.frame(Gene = row.names(c24a), c24a)

list_of_datasets <- list("c0a" = c0a, "c1a" = c1a, "c2a" = c2a, "c3a" = c3a, "c4a" = c4a, "c5a" = c5a, "c6a" = c6a, "c7a" = c7a, "c8a" = c8a, "c9a" = c9a, "c10a" = c10a, "c11a" = c11a, "c12a" = c12a, "c13a" = c13a, "c14a" = c14a, "c15a" = c15a, "c16a" = c16a, "c17a" = c17a, "c18a" = c18a, "c19a" = c19a, "c20a" = c20a, "c21a" = c21a, "c22a" = c22a, "c23a" = c23a, "c24a" = c24a)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.Average.Exp.Genes.xlsx")


#For Differential Expression between Control and OE of the same cluster (ie Cluster6_CTRL vs Cluster6_KO)
immune.combined <- Data.Combined
immune.combined$celltype.OE <- paste(Idents(immune.combined), immune.combined$OE, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.OE"
head(immune.combined@meta.data)

new.cluster.idss <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c22", "c23", "c24")
#for finding DEG between the groups of the same cluster. 
for(i in new.cluster.idss) {
  a = paste(i, "OE", sep = "_")
  print(a)
  b = paste(i, "WT", sep = "_")
  c = paste(i, "OE.SHAM", sep = "_")
  d = paste(i, "WT.SHAM", sep = "_")
  MC.response <- FindMarkers(immune.combined, ident.1 = a, ident.2 = c, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
  head(MC.response, n = 10)
  str1 = i
  str2 = "comparison.csv"
  result = paste(str1,str2, sep = "")
  write.csv(MC.response, file = result)
}

C0 <- read.csv("c0comparison.csv")
C1 <- read.csv("c1comparison.csv")
C2 <- read.csv("c2comparison.csv")
C3 <- read.csv("c3comparison.csv")
C4 <- read.csv("c4comparison.csv")
C5 <- read.csv("c5comparison.csv")
C6 <- read.csv("c6comparison.csv")
C7 <- read.csv("c7comparison.csv")
C8 <- read.csv("c8comparison.csv")
C9 <- read.csv("c9comparison.csv")
C10 <- read.csv("c10comparison.csv")
C11 <- read.csv("c11comparison.csv")
C12 <- read.csv("c12comparison.csv")
C13 <- read.csv("c13comparison.csv")
C14 <- read.csv("c14comparison.csv")
C15 <- read.csv("c15comparison.csv")
C16 <- read.csv("c16comparison.csv")
C17 <- read.csv("c17comparison.csv")
C18 <- read.csv("c18comparison.csv")
C19 <- read.csv("c19comparison.csv")
C20 <- read.csv("c20comparison.csv")
#C21 <- read.csv("c21comparison.csv")
C22 <- read.csv("c22comparison.csv")
C23 <- read.csv("c23comparison.csv")
C24 <- read.csv("c24comparison.csv")

C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
#C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)
C23 <- data.frame(Gene = row.names(C23), C23)
C24 <- data.frame(Gene = row.names(C24), C24)

list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20,  "C22" = C22, "C23" = C23, "C24" = C24)
write.xlsx(list_of_datasets, file = "SummaryAllgroups.DEG OEvOESham.Genes.xlsx")


#For Differential Expression between Control and KO of the same cluster (ie Cluster6_CTRL vs Cluster6_KO)
immune.combined <- Data.Combined
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

#Cluster 0
C0 <- FindMarkers(immune.combined, ident.1 = "0_KO", ident.2 = "0_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C0, n = 15)
#Cluster 1
C1 <- FindMarkers(immune.combined, ident.1 = "1_KO", ident.2 = "1_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C1, n = 15)
#Cluster 2
C2 <- FindMarkers(immune.combined, ident.1 = "2_KO", ident.2 = "2_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C2, n = 15)
#Cluster 3
C3 <- FindMarkers(immune.combined, ident.1 = "3_KO", ident.2 = "3_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C3, n = 15)
#Cluster 4
C4 <- FindMarkers(immune.combined, ident.1 = "4_KO", ident.2 = "4_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C4, n = 15)
#Cluster 5
C5 <- FindMarkers(immune.combined, ident.1 = "5_KO", ident.2 = "5_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C5, n = 15)
#Cluster 6
C6 <- FindMarkers(immune.combined, ident.1 = "6_KO", ident.2 = "6_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C6, n = 15)
#Cluster 7
C7 <- FindMarkers(immune.combined, ident.1 = "7_KO", ident.2 = "7_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C7, n = 15)
#Cluster 8
C8 <- FindMarkers(immune.combined, ident.1 = "8_KO", ident.2 = "8_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C8, n = 15)
#Cluster 9
C9 <- FindMarkers(immune.combined, ident.1 = "9_KO", ident.2 = "9_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C9, n = 15)
#Cluster 10
C10 <- FindMarkers(immune.combined, ident.1 = "10_KO", ident.2 = "10_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C10, n = 15)
#Cluster 11
C11 <- FindMarkers(immune.combined, ident.1 = "11_KO", ident.2 = "11_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C11, n = 15)
#Cluster 12
C12 <- FindMarkers(immune.combined, ident.1 = "12_KO", ident.2 = "12_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C12, n = 15)
#Cluster 13
C13 <- FindMarkers(immune.combined, ident.1 = "13_KO", ident.2 = "13_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C13, n = 15)
#Cluster 14
C14 <- FindMarkers(immune.combined, ident.1 = "14_KO", ident.2 = "14_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C14, n = 15)
#Cluster 15
C15 <- FindMarkers(immune.combined, ident.1 = "15_KO", ident.2 = "15_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C15, n = 15)
#Cluster 16
C16 <- FindMarkers(immune.combined, ident.1 = "16_KO", ident.2 = "16_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C16, n = 15)
#Cluster 17
C17 <- FindMarkers(immune.combined, ident.1 = "17_KO", ident.2 = "17_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C17, n = 15)
#Cluster 18
C18 <- FindMarkers(immune.combined, ident.1 = "18_KO", ident.2 = "18_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C18, n = 15)
#Cluster 19
C19 <- FindMarkers(immune.combined, ident.1 = "19_KO", ident.2 = "19_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C19, n = 15)
#Cluster 20
C20 <- FindMarkers(immune.combined, ident.1 = "20_KO", ident.2 = "20_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C20, n = 15)
#Cluster 21
C21 <- FindMarkers(immune.combined, ident.1 = "21_KO", ident.2 = "21_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C21, n = 15)
#Cluster 22
C22 <- FindMarkers(immune.combined, ident.1 = "22_KO", ident.2 = "22_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C22, n = 15)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)

#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" =C20, "C21" =C21, "C22" = C22)
write.xlsx(list_of_datasets, file = "6Samp_DiffExprGenes_CTRLvsKO_Updated2.xlsx")

#For Differential Expression between a single cluster and all other cells (Cluster 1 vs all other clusters)
immune.combined <- Data.Combined
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype"

C0 <- FindMarkers(immune.combined, ident.1 = "0", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(immune.combined, ident.1 = "6", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C7 <- FindMarkers(immune.combined, ident.1 = "7", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C8 <- FindMarkers(immune.combined, ident.1 = "8", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C9 <- FindMarkers(immune.combined, ident.1 = "9", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C10 <- FindMarkers(immune.combined, ident.1 = "10", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C11 <- FindMarkers(immune.combined, ident.1 = "11", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C12 <- FindMarkers(immune.combined, ident.1 = "12", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C13 <- FindMarkers(immune.combined, ident.1 = "13", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C14 <- FindMarkers(immune.combined, ident.1 = "14", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C15 <- FindMarkers(immune.combined, ident.1 = "15", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C16 <- FindMarkers(immune.combined, ident.1 = "16", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C17 <- FindMarkers(immune.combined, ident.1 = "17", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C18 <- FindMarkers(immune.combined, ident.1 = "18", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C19 <- FindMarkers(immune.combined, ident.1 = "19", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C20 <- FindMarkers(immune.combined, ident.1 = "20", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C21 <- FindMarkers(immune.combined, ident.1 = "21", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C22 <- FindMarkers(immune.combined, ident.1 = "22", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)


#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)

#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20, "C21" = C21, "C22" = C22)
write.xlsx(list_of_datasets, file = "6Sample_ClusterMarkers_Updated2.xlsx")



#load saved file
Data.Combined <- readRDS("6Samples_subclustered.rds")

# Cell-Cycle scoring
# read in a list of cell cycle markers, from Tirosh et al, 2015
cellcyclegenes <- "regev_lab_cell_cycle_genes.lowercase.txt"
cc.genes <- readLines(con=cellcyclegenes)
length(cc.genes)
# segregate this list into markers of G2/M phase and markers of S
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]
DefaultAssay(object = Data.Combined) <- "RNA"
# assign cell-cycle scores
Data.Combined <- CellCycleScoring(Data.Combined, s.features = s.genes, g2m.features = g2m.genes, set.ident= FALSE)
head(Data.Combined@meta.data)
#Visualize results over umap
DimPlot(Data.Combined, reduction = "umap", group.by = "Phase", split.by = "OE")
#Get number of cells in each phase separated by cluster and group
table(Data.Combined@meta.data[["Phase"]], Data.Combined@meta.data$sub_cluster, Data.Combined@meta.data$OE)



Sub2 <- subset(immune.combined, idents = c("2.CTRL", "2.KO"))
saveRDS(Sub2, file = "PECSubcluster2_Only.rds")

Sub3 <- subset(immune.combined, idents = c("3.CTRL", "3.KO"))
saveRDS(Sub3, file = "PECSubcluster3_Only.rds")


###########################################################
###Subclustering of PEC cluster (Cluster 3 in 6 Sample)####
Data.Combined.c3 <- subset(Data.Combined, idents = "3")


DefaultAssay(Data.Combined) <- "integrated"
Data.Combined <- FindNeighbors(Data.Combined, dims = 1:10)
Data.Combined <- FindClusters(Data.Combined, resolution = 0.5)

# Generate a new column called sub_cluster in the metadata
Data.Combined$sub_cluster <- as.character(Idents(Data.Combined))

head(Data.Combined@meta.data)
# Change the information of cells containing sub-cluster information
Data.Combined$sub_cluster[Cells(Data.Combined)] <- paste("Novel/PTS3",Idents(Data.Combined))
DimPlot(Data.Combined, group.by = "sub_cluster")

# Visualization of UMAP
DimPlot(Data.Combined, reduction = "umap", group.by = "sub_cluster")
DimPlot(Data.Combined, reduction = "umap", split.by = "KO")
DimPlot(Data.Combined, reduction = "umap", label = FALSE)
DimPlot(Data.Combined, reduction = "umap", label = TRUE)
DimPlot(Data.combined, reduction = "umap", label = FALSE)

DefaultAssay(Data.Combined) <- "RNA"
Idents(Data.Combined) <- "celltype"
FeaturePlot(Data.Combined, features = "Egfr", split.by = "OE")
VlnPlot(Data.Combined, features = c("Itgav"), ncol = 1, split.by = "KO", pt.size = 0)+facet_wrap(~ "Receptors")
VlnPlot(Data.Combined, features = c("Camk1d", "Acy3"), ncol = 2, split.by = "OE", pt.size = 0)+facet_wrap(~ "Receptors")


VlnPlot(Data.Combined.c3, features = c("Cd44", "Slc5a12", "Havcr1", "Nphs1", "Sytl2"), ncol = 1, split.by = "KO", pt.size = 0)
RidgePlot(Data.Combined.c3, features = "Cd44", sort = "decreasing")
FeatureScatter(Data.Combined, feature1= "Cd36", feature2 = "Slc16a2")

VlnPlot(Data.Combined.c3, features = c("Cd44", "Slc5a12", "Slc7a7", "Havcr1"), ncol = 1, pt.size = 0)

#For Differential Expression between Control and KO of the same cluster
immune.combined <- Data.Combined.c3
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = ".")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

C0 <- FindMarkers(Immune.Combined, ident.1 = "0_KO", ident.2 = "0_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(Immune.Combined, ident.1 = "1_KO", ident.2 = "1_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(Immune.Combined, ident.1 = "2_KO", ident.2 = "2_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(Immune.Combined, ident.1 = "3_KO", ident.2 = "3_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(Immune.Combined, ident.1 = "4_KO", ident.2 = "4_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(Immune.Combined, ident.1 = "5_KO", ident.2 = "5_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(Immune.Combined, ident.1 = "6_KO", ident.2 = "6_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C7 <- FindMarkers(Immune.Combined, ident.1 = "7_KO", ident.2 = "7_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C8 <- FindMarkers(immune.combined, ident.1 = "8_KO", ident.2 = "8_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)

#For Differential Expression between a single cluster and all other cells (Cluster 1 vs all other clusters)
immune.combined <- Data.Combined.c3
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype"

C0 <- FindMarkers(Immune.Combined, ident.1 = "0", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(Immune.Combined, ident.1 = "1", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(Immune.Combined, ident.1 = "2", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(Immune.Combined, ident.1 = "3", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(Immune.Combined, ident.1 = "4", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(Immune.Combined, ident.1 = "5", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(Immune.Combined, ident.1 = "6", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C7 <- FindMarkers(Immune.Combined, ident.1 = "7", ident.2 = NULL, min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C8 <- FindMarkers(immune.combined, ident.1 = "8", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C6), C6)


#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7)
write.xlsx(list_of_datasets, file = "PTS3andNovelsub_markers.xlsx")

#To output average expression split by group
cluster.averages <- AverageExpression(immune.combined, assays = "RNA", add.ident = "OE")
write.csv(cluster.averages[["RNA"]], file = "C0_subclusters_CAunsplit.csv")

# Get number of cells per cluster and per sample/group of origin
table(immune.combined@meta.data[["seurat_clusters"]], immune.combined@meta.data$orig.ident)
table(Data.Combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$OE)



#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7)
write.xlsx(list_of_datasets, file = "PEC_subclusters_CTRLvsKO.xlsx")



###########################################################
###Subclustering of Podocytes cluster (Cluster 15 in 6 Sample)####
Data.Combined.c15 <- subset(Data.Combined, idents = "15")
DefaultAssay(object = Data.Combined) <- "integrated"
Data.Combined <- FindNeighbors(Data.Combined, dims = 1:10)
Data.Combined <- FindClusters(Data.Combined, resolution = 0.5)

DimPlot(Data.Combined, reduction = "umap")

# Generate a new column called sub_cluster in the metadata
Data.Combined$sub_cluster <- as.character(Idents(Data.Combined))

# Change the information of cells containing sub-cluster information
Data.Combined$sub_cluster[Cells(Data.Combined.c15)] <- paste("c15",Idents(Data.Combined.c15))
DimPlot(Data.Combined.c15, group.by = "sub_cluster")

DefaultAssay(Data.Combined.c15) <- "RNA"
FeaturePlot(Data.Combined.c15, features = "Nphs1")
VlnPlot(Data.Combined.c15, features = c("Nphs1", "Nphs2", "Fn1", "Stat3"), ncol = 1, split.by = "KO", pt.size = 0)
RidgePlot(Data.Combined.c15, features = "Klf4", sort = "decreasing")
FeatureScatter(Data.Combined.c15, feature1= "Nphs1", feature2 = "Nphs2")

#For Differential Expression between Control and KO of the same cluster
immune.combined <- Data.Combined.c15
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

C0 <- FindMarkers(immune.combined, ident.1 = "0_KO", ident.2 = "0_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1_KO", ident.2 = "1_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2_KO", ident.2 = "2_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3_KO", ident.2 = "3_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4_KO", ident.2 = "4_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5_KO", ident.2 = "5_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)

#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5)
write.xlsx(list_of_datasets, file = "Podocyte_subclusters_CTRLvsKO.xlsx")

#For Differential Expression between a single cluster and all other cells (Cluster 1 vs all other clusters)
immune.combined <- Data.Combined.c15
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype"

C0 <- FindMarkers(immune.combined, ident.1 = "0", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(immune.combined, ident.1 = "6", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)


#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5)
write.xlsx(list_of_datasets, file = "Podocyte_subcluster_markers.xlsx")

#To output average expression split by group
cluster.averages <- AverageExpression(immune.combined, assays = "RNA", add.ident = "KO")
write.csv(cluster.averages[["RNA"]], file = "Podocyte_subclusters_ClusterAverages.csv")

# Get number of cells per cluster and per sample/group of origin
table(immune.combined@meta.data[["seurat_clusters"]], immune.combined@meta.data$orig.ident)
table(immune.combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$KO)


