#Load libraries
library(dplyr)
library(ggplot2)
library(Seurat)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(celldex)
library(R.matlab)
library(scater)
library(scran)
library(gt)
library(topGO)
library(org.Mm.eg.db)
library(RColorBrewer)

#Import filtered and raw data frames
#readRDS("~/Documents/SingleCell/BrennanPaper/brenNeu.rds")
#Get file names in directory for analysis
setwd("~/Desktop/Thesis/SingleCell/BrennanPaper")
tenXDat <- list.dirs(recursive = F)

#Read in all files and turn into seurat objects
matricies_B <- sapply(tenXDat, FUN = Read10X)
seuObj_B <- lapply(matricies_B, CreateSeuratObject)


#Assign mitochondrial gene percentages to cells
seuObj_a <- lapply(X = seuObj_B, FUN = function(x) {
  x[["percent_mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x
})


#Subset objects to expected expression levels
seuObj_a <- lapply(X = seuObj_a, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 15)
})

#Normalize cells and find variable genes
seuObj_a <- lapply(X = seuObj_a, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


seuObj_a[[1]]@meta.data$time = rep(7, 478)
seuObj_a[[2]]@meta.data$time = rep(7, 6327)
seuObj_a[[3]]@meta.data$time = rep(0, 308)
seuObj_a[[4]]@meta.data$time = rep(0, 1899)
seuObj_a[[5]]@meta.data$time = rep(28, 1162)
seuObj_a[[6]]@meta.data$time = rep(28, 11309)
seuObj_a[[1]]@meta.data$treatment = rep("P", 478)
seuObj_a[[2]]@meta.data$treatment = rep("V", 6327)
seuObj_a[[3]]@meta.data$treatment = rep("P", 308)
seuObj_a[[4]]@meta.data$treatment = rep("V", 1899)
seuObj_a[[5]]@meta.data$treatment = rep("P", 1162)
seuObj_a[[6]]@meta.data$treatment = rep("V", 11309)

anchors <- FindIntegrationAnchors(object.list = seuObj_a)
seuObj.integrated_B <- IntegrateData(anchorset = anchors)
DefaultAssay(seuObj.integrated_B) <- "integrated"

FeatureScatter(seuObj.integrated_B, feature1 = "nCount_RNA", feature2 = "percent_mt")
FeatureScatter(seuObj.integrated_B, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seuObj.integrated_B <- ScaleData(seuObj.integrated_B, verbose = FALSE)
seuObj.integrated_B <- RunPCA(seuObj.integrated_B, verbose = FALSE)

ElbowPlot(seuObj.integrated_B, ndims = 40)

#seuObj.integrated_B <- RunUMAP(seuObj.integrated_B, dims = 1:35)
DimHeatmap(seuObj.integrated_B, dims = 45:50, cells = 500, balanced = TRUE)

seuObj.integrated_B <- RunUMAP(seuObj.integrated_B, dims = 1:40)
DimPlot(seuObj.integrated_B, group.by = "time")

seuObj.integrated_B <- FindNeighbors(seuObj.integrated_B, dims = 1:30)
seuObj.integrated_B <- FindClusters(seuObj.integrated_B)

DimPlot(seuObj.integrated_B, reduction = "umap", label = T)+ NoLegend()

#Cluster 12 appears to be neutrophils
VlnPlot(seuObj.integrated_B, features = c("S100a9"))
FeaturePlot(seuObj.integrated_B, features = c("Ly6g"))
FeaturePlot(seuObj.integrated_B, features = c("S100a9"))
FeaturePlot(seuObj.integrated_B, features = c("S100a8"))
FeaturePlot(seuObj.integrated_B, features = c("Mmp9"))
FeaturePlot(seuObj.integrated_B, features = c("Il1r2"))


#Celldex and singleR to label
ref.se <- ImmGenData()
sce <- as.SingleCellExperiment(seuObj.integrated_B)
pred <- SingleR(test = sce, ref = ref.se, labels = ref.se$label.main)

seuObj.integrated_B$pruned_labels = pred$pruned.labels

as.data.frame(table(seuObj.integrated_B$time)) %>% gt() 

DimPlot(seuObj.integrated_B, group.by = "pruned_labels", label = T, 
        label.size = 3) + NoLegend() 

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

tiff("BrenUmap.tiff", units = "in", width = 5, height = 5, res = 300)
DimPlot(seuObj.integrated_B, group.by = "pruned_labels", label = F,
        cols =col_vector) +
  guides(color=guide_legend(ncol =1, override.aes = list(size=2))) +
  ggtitle("Brennan: Cell Type UMAP")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text=element_text(size=10))

dev.off()

saveRDS(seuObj.integrated_B, file = "~/Desktop/Thesis/SingleCell/BrennanPaper/Bren_Dat.rds")

#subset neutrophils
Neutrophils_Bren <- subset(seuObj.integrated_B, pruned_labels == "Neutrophils" & seurat_clusters == 12)
DefaultAssay(Neutrophils_Bren) <- "integrated"


table(Neutrophils_Bren$time)
Neutrophils_Bren$time= as.factor(Neutrophils_Bren$time)
Idents(Neutrophils_Bren) <- 'time'
DimPlot(Neutrophils_Bren, split.by = "time") 
DimPlot(Neutrophils_Bren, group.by = "time") 

Neutrophils_Bren <- FindVariableFeatures(Neutrophils_Bren, selection.method = "vst", assay = "RNA")
all.genes <- rownames(Neutrophils_Bren)
Neutrophils_Bren <- ScaleData(Neutrophils_Bren, features = all.genes)
Neutrophils_Bren <- RunPCA(Neutrophils_Bren, features = VariableFeatures(object = Neutrophils_Bren))
ElbowPlot(Neutrophils_Bren, ndims = 30)

Neutrophils_Bren <- FindNeighbors(Neutrophils_Bren, dims = 1:30)
Neutrophils_Bren <- FindClusters(Neutrophils_Bren)

setwd("~/Desktop/GitHub/Neutrophils/Data")
tiff("BrenNeutrophils.tif", units = "in", width = 6, height = 6, res = 300)
DimPlot(Neutrophils_Bren, group.by = "time", label = F, cols = c("red", "gold", "aquamarine4")) +
  guides(color=guide_legend(ncol =1, override.aes = list(size=5))) +
  ggtitle("Brennan: Neutrophils UMAP")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))
dev.off()

saveRDS(Neutrophils_Bren, file = "~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/Neutrophils_Bren.rds")

FeaturePlot(Neutrophils_Bren, features = c("Ly6g"))
FeaturePlot(Neutrophils_Bren, features = c("S100a9"))
FeaturePlot(Neutrophils_Bren, features = c("S100a8"))
FeaturePlot(Neutrophils_Bren, features = c("Mmp9"))


Idents(Neutrophils_Bren) <- 'time'
Neu.Timemarkers <- FindAllMarkers(Neutrophils_Bren, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)

Neu.Timemarkers$pct.diff <- Neu.Timemarkers$pct.1 - Neu.Timemarkers$pct.2



#### GO ANALYSIS ###
Uninjured <- Neu.Timemarkers[Neu.Timemarkers$cluster == 0,]
UninjuredGenes = unlist(split(Uninjured$p_val_adj, rownames(Uninjured)))
UninjuredGenes <- UninjuredGenes[UninjuredGenes <0.05]
write.csv(c(names(UninjuredGenes), UninjuredGenes), "Uninjured", row.names = F, quote = F)

Seven <- Neu.Timemarkers[Neu.Timemarkers$cluster == 7,]
Seven = unlist(split(Seven$p_val_adj, rownames(Seven)))
Seven <- Seven[Seven <0.05]
write.csv(names(Seven), "7dpi", row.names = F, quote = F)

TwentyEight <- Neu.Timemarkers[Neu.Timemarkers$cluster == 28,]
TwentyEight = unlist(split(TwentyEight$p_val_adj, rownames(TwentyEight)))
TwentyEight <- TwentyEight[TwentyEight <0.05]
write.csv(names(TwentyEight), "28dpi", row.names = F, quote = F)



bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", 
                  "Retnlg", "Lcn2", "Olfm4", "Chil3", 
                  "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", 
                  "Lyz1"))

neuMaturationGenes <- read_csv("~/Documents/SingleCell/neuMaturationGenes.csv")
maturationGenes <- list(neuMaturationGenes[[1]])

Neutrophils_Bren <- AddModuleScore(Neutrophils_Bren, features = bmGenes, search = T, assay = "RNA",
                                   name = "BoneMarrow")

Neutrophils_Bren <- AddModuleScore(Neutrophils_Bren, features = maturationGenes, search = T, assay = "RNA",
                                   name = "Maturation")

Neutrophils_Bren[['time']] <- c(rep(7, 245), rep(0, 80), rep(28, 300))
Neutrophils_Bren$time <- factor(Neutrophils_Bren$time, levels = c(0,7,28))
Idents(Neutrophils_Bren) <- "time"
VlnPlot(Neutrophils_Bren, features = "BoneMarrow1") + ggtitle("Bone Marrow Brennan")

VlnPlot(Neutrophils_Bren, features = "Maturation1") + ggtitle("Maturation Brennan")

bmDat <- data.frame("time" = neutrophils$time, "bm" = neutrophils$BoneMarrow1)
bmDat %>% group_by(time) %>% summarise(bm = mean(bm))

matDat <- data.frame("time" = neutrophils$time, "Maturation" = neutrophils$Maturation1)
matDat %>% group_by(time) %>% summarise(mature = mean(Maturation))

saveRDS(Neutrophils_Bren, file = "~/Documents/SingleCell/BrennanPaper/BrennanDat.rds")


#Neutrotime early and late based on fig 2
earlyGenes <- list(c("Lyz2", "Retnlg", "Lgals3", "Mmp8", "Ifitm6", "Ly6g", "Pfn1", "Cybb", "Cd177", "Mgst1", "Pglyrp1",
                     "Ltf", "Prdx5", "Lcn2", "Ngp", "Adpgk", "Anxa1", "Arhgdib", "Ly6c2", "Camp", "Dstn", "Chil3",
                     "Wfdc21", "Serpinb1a"))

lateGenes <- list(c("Il1b", "Ccl6", "Fth1", "H2-D1", "Dusp1", "Junb", "Ifitm1", "Fxyd5", "Ifitm2", "Malat1", "Btg1",
                    "Tyrobp", "Jund", "Ftl1", "Srgn", "Csf3r", "Wfdc17", "Rps27", "Msrb1", "Fau", "Rps9"))

Neutrophils_Bren <- AddModuleScore(Neutrophils_Bren, features = earlyGenes, search = T, assay = "RNA",
                                   name = "neutrotimeEarly") 

Neutrophils_Bren <- AddModuleScore(Neutrophils_Bren, features = lateGenes, search = T, assay = "RNA",
                                   name = "neutrotimeLate")

VlnPlot(Neutrophils_Bren, features = "neutrotimeEarly1", group.by = "time")+ ggtitle("Early neutrotime")
VlnPlot(Neutrophils_Bren, features = "neutrotimeLate1", group.by = "time")+ ggtitle("Late neutrotime")
# Idents(nr) = "timepoint"
# nr.markers <- FindAllMarkers(nr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# nr.markers[["pct_diff"]]= nr.markers$pct.2 - nr.markers$pct.1
# nr.markers <- nr.markers[nr.markers$p_val_adj < 0.01,]
# 
# nr0 <- FindMarkers(nr, ident.1 = 0, only.pos = TRUE)
# nr0["pct.diff"] = nr0$pct.1 - nr0$pct.2
# #nr0 <- nr0[nr0$p_val_adj<0.05,]
# 
# nr7 <- FindMarkers(nr, ident.1 = 7, only.pos = TRUE)
# nr7["pct.diff"] = nr7$pct.1 - nr7$pct.2
# #nr7 <- nr7[nr7$p_val_adj<0.05,]
# 
# nr28 <- FindMarkers(nr, ident.1 = 28, only.pos = TRUE)
# nr28["pct.diff"] = nr28$pct.1 - nr28$pct.2
# #nr28 <- nr28[nr28$p_val_adj<0.05,]
# 
# geneListB = unlist(split(nr0$p_val_adj, rownames(nr0)))
# #geneListB <- ifelse(geneListB < 0.05, 1, 0)
# 
# nr0Go <- new("topGOdata",
#              description = "nr0",
#              ontology = "BP",
#              allGenes = geneListB,
#              geneSel = function(x)(x < 0.01),
#              mapping = "org.Mm.eg.db",
#              annot = annFUN.org,
#              ID = "symbol")
# 
# #Different gene selection method
# # nr0Go <- new("topGOdata",
# #              description = "nr0",
# #              ontology = "BP",
# #              allGenes = geneListB,
# #              geneSel = topDiffGenes,
# #              mapping = "org.Mm.eg.db",
# #              annot = annFUN.org,
# #              ID = "symbol")
# 
# resultFisher0 <- runTest(nr0Go, algorithm = "classic", statistic = "fisher")
# 
# resultKS <- runTest(nr0Go, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(nr0Go, algorithm = "elim", statistic = "ks")
# 
# nr0Tab <- GenTable(nr0Go, classicFisher = resultFisher0,
#          elimKS = resultKS.elim,
#          classicKS = resultKS,
#          topNodes = 20)
# 
# as.data.frame((nr0Tab[c(2, 7, 8,9)])) %>% gt() 
# 
# write.csv(names(geneListB[geneListB<0.01]), "uninjured", row.names = F, quote = F)
# 
# 
# geneList7 = unlist(split(nr7$p_val_adj, rownames(nr7)))
# #geneList7 <- ifelse(geneList7 < 0.01, 1, 0)
# 
# nr7Go <- new("topGOdata",
#              description = "nr7",
#              ontology = "BP",
#              allGenes = geneList7,
#              geneSel = function(x)(x < 0.01),
#              mapping = "org.Mm.eg.db",
#              annot = annFUN.org,
#              ID = "symbol")
# 
# resultFisher7 <- runTest(nr7Go, algorithm = "classic", statistic = "fisher")
# GenTable(nr0Go, classicFisher = resultFisher7)
# 
# write.csv(names(geneList7[geneList7 > 0]), "7dpi", row.names = F, quote = F)
# 
# 
# 
# geneList28 = unlist(split(nr28$p_val_adj, rownames(nr28)))
# #geneList28 <- ifelse(geneList28 < 0.05, 1, 0)
# 
# nr28Go <- new("topGOdata",
#              description = "nr28",
#              ontology = "BP",
#              allGenes = geneList28,
#              geneSel = function(x)(x < 0.01),
#              mapping = "org.Mm.eg.db",
#              annot = annFUN.org,
#              ID = "symbol")
# 
# resultFisher28 <- runTest(nr28Go, algorithm = "classic", statistic = "fisher")
# resultKS <- runTest(nr28Go, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(nr28Go, algorithm = "elim", statistic = "ks")
# 
# nr28Tab <- GenTable(nr28Go, classicFisher = resultFisher28,
#                    elimKS = resultKS.elim,
#                    classicKS = resultKS,
#                    topNodes = 20)
# 
# as.data.frame((nr28Tab[c(2, 7, 8,9)])) %>% gt() 
# 
# write.csv(names(geneList28[geneList28 < 0.01]), "28dpi", row.names = F, quote = F)
# 

#Bone marrow scores