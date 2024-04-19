library(Seurat)
library(SingleR)
library(celldex)
library(org.Mm.eg.db)
library(speckle)
library(SingleCellExperiment)
library(dplyr)
library(EnhancedVolcano)

#Skip the steps ahead of loading data by just loading saved rds object
seuObj_Wang <- readRDS("~/Documents/SingleCell/WangPaper/wangSeuObj.rds")
#Best place to start to skip prep steps
seuObj.integrated <- readRDS("~/Documents/SingleCell/WangPaper/seuObj.integratedWang.rds")
#Load in single cell data
setwd("~/Desktop/Thesis/SingleCell/WangPaper/data")
tenXDat <- list.dirs(recursive = F)

matricies_W <- sapply(tenXDat, FUN = Read10X)
#options(Seurat.object.assay.version = "v3")
seuObj_W <- lapply(matricies_W, CreateSeuratObject)



#Assign mitochondrial gene percentages to cells
seuObj_Wang <- lapply(X = seuObj_W, FUN = function(x) {
  x[["percent_mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x
})


#Subset objects to expected expression levels
seuObj_Wang <- lapply(X = seuObj_Wang, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 15)
})

#Normalize cells and find variable genes
seuObj_Wang <- lapply(X = seuObj_Wang, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seuObj_Wang)

seuObj_Wang <- lapply(X = seuObj_Wang, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


#Assign each sample the correct day of collection (0 for uninjured mice)
seuObj_Wang[[1]]@meta.data$time = rep(14, 10015)
seuObj_Wang[[2]]@meta.data$time = rep(14, 10029)
seuObj_Wang[[3]]@meta.data$time = rep(14, 9993)
seuObj_Wang[[4]]@meta.data$time = rep(14, 10061)
seuObj_Wang[[5]]@meta.data$time = rep(0, 3151)
seuObj_Wang[[6]]@meta.data$time = rep(0, 737)
seuObj_Wang[[7]]@meta.data$time = rep(0, 741)
seuObj_Wang[[8]]@meta.data$time = rep(0, 743)
seuObj_Wang[[9]]@meta.data$time = rep(0, 744)
seuObj_Wang[[10]]@meta.data$time = rep(3, 9278)
seuObj_Wang[[11]]@meta.data$time = rep(3, 9252)
seuObj_Wang[[12]]@meta.data$time = rep(3, 9296)
seuObj_Wang[[13]]@meta.data$time = rep(3, 9294)


#Using one uninjured mouse, one 3dpi, and one 14 dpi mouse for references for integration
anchors <- FindIntegrationAnchors(object.list = seuObj_Wang, reference = c(1,5,10), 
                                  reduction = "rpca")


seuObj.integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(seuObj.integrated) <- "integrated"


seuObj.integrated <- ScaleData(seuObj.integrated, verbose = FALSE)
seuObj.integrated <- RunPCA(seuObj.integrated, verbose = FALSE)

ElbowPlot(seuObj.integrated, ndims = 40)

seuObj.integrated <- RunUMAP(seuObj.integrated, dims = 1:30)

seuObj.integrated <- FindNeighbors(seuObj.integrated, dims = 1:30)
seuObj.integrated <- FindClusters(seuObj.integrated)

DimPlot(seuObj.integrated, group.by = "time")

DimPlot(seuObj.integrated, reduction = "umap", label = T)+ NoLegend()

#Identifying Neutrophil cluster
VlnPlot(seuObj.integrated, features = c("S100a9"))
FeaturePlot(seuObj.integrated, features = c("Ly6g"))
FeaturePlot(seuObj.integrated, features = c("S100a9"))
FeaturePlot(seuObj.integrated, features = c("S100a8"))
FeaturePlot(seuObj.integrated, features = c("Mmp9"))
FeaturePlot(seuObj.integrated, features = c("Il1r2"))


# Generate cell type predictions from database
ref.se <- ImmGenData()
sce <- as.SingleCellExperiment(seuObj.integrated)
pred <- SingleR(test = sce, ref = ref.se, labels = ref.se$label.main)

# Label clusters with the predicted cell types
seuObj.integrated$pruned_labels = pred$pruned.labels

DimPlot(seuObj.integrated, group.by = "pruned_labels", label = T, 
        label.size = 3) + NoLegend() 

setwd("~/Desktop/GitHub/Neutrophils/Data")
tiff("WangUMAP.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(seuObj.integrated, group.by = "pruned_labels", label = F) +
  guides(color=guide_legend(ncol =1, override.aes = list(size=5))) +
  ggtitle("Wang: Cell Type UMAP")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))
dev.off()

saveRDS(seuObj.integrated, file = "~/Desktop/Thesis/SingleCell/WangPaper/WangDat.rds")

#Subset neutrophils
Neutrophils_Wang <- subset(seuObj.integrated, pruned_labels == "Neutrophils")

#Neutrophils_Wang <- subset(Neutrophils_Wang, seurat_clusters %in% c(1,11,13,14,26,28,34))
DefaultAssay(Neutrophils_Wang) <- "integrated"
table(Neutrophils_Wang$seurat_clusters)

#Possible neutrophil progenitor marker/ very early marker
FeaturePlot(Neutrophils_Wang, features = c("Fcnb"))

DimPlot(Neutrophils_Wang, group.by = "time")+
  guides(color=guide_legend(ncol =1, override.aes = list(size=5))) +
  ggtitle("Wang: Neutrophils by Time")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))

set.seed(42)
Neutrophils_Wang <- FindVariableFeatures(Neutrophils_Wang, selection.method = "vst", assay = "RNA")
all.genes <- rownames(Neutrophils_Wang)
Neutrophils_Wang <- ScaleData(Neutrophils_Wang, features = all.genes)
Neutrophils_Wang <- RunPCA(Neutrophils_Wang, features = VariableFeatures(object = Neutrophils_Wang))

ElbowPlot(Neutrophils_Wang)


Neutrophils_Wang <- FindNeighbors(Neutrophils_Wang, dims = 1:20)
Neutrophils_Wang <- FindClusters(Neutrophils_Wang)

Idents(Neutrophils_Wang) <- 'seurat_clusters'
Neutrophils_Wang <- RunUMAP(Neutrophils_Wang, dims = 1:20)
DimPlot(Neutrophils_Wang, reduction = "umap")

DimPlot(Neutrophils_Wang, group.by = "time") +
  ggtitle("Wang Neutrophils by Time")

FeaturePlot(Neutrophils_Wang, features = c("Ly6g"))
FeaturePlot(Neutrophils_Wang, features = c("S100a9"))
FeaturePlot(Neutrophils_Wang, features = c("S100a8"))
FeaturePlot(Neutrophils_Wang, features = c("Mmp9"))

Neu.markers <- FindAllMarkers(Neutrophils_Wang, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Neu.markers$pct.diff <- Neu.markers$pct.1 - Neu.markers$pct.2

Idents(Neutrophils_Wang) <- 'time'
Neu.Timemarkers <- FindAllMarkers(Neutrophils_Wang, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.15)

Neu.Timemarkers$pct.diff <- Neu.Timemarkers$pct.1 - Neu.Timemarkers$pct.2



#### GO ANALYSIS ###
Uninjured <- Neu.Timemarkers[Neu.Timemarkers$cluster == 0,]
UninjuredGenes = unlist(split(Uninjured$p_val_adj, rownames(Uninjured)))
UninjuredGenes <- UninjuredGenes[UninjuredGenes <0.05]
write.csv(c(names(UninjuredGenes), UninjuredGenes), "Uninjured", row.names = F, quote = F)

Three <- Neu.Timemarkers[Neu.Timemarkers$cluster == 3,]
Three = unlist(split(Three$p_val_adj, rownames(Three)))
Three <- Three[Three <0.05]
write.csv(names(Three), "3dpi", row.names = F, quote = F)

Fourteen <- Neu.Timemarkers[Neu.Timemarkers$cluster == 14,]
Fourteen = unlist(split(Fourteen$p_val_adj, rownames(Fourteen)))
Fourteen <- Fourteen[Fourteen <0.05]
write.csv(names(Fourteen), "14dpi", row.names = F, quote = F)

#Why do so few express this y gene?

DimPlot(Neutrophils_Wang, group.by = "time") 


bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", 
                  "Retnlg", "Lcn2", "Olfm4", "Chil3", 
                  "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", 
                  "Lyz1"))

Neutrophils_Wang <- AddModuleScore(Neutrophils_Wang, features = bmGenes, search = T, assay = "RNA",
                                   name = "BoneMarrow")

Neutrophils_Wang$time <- factor(Neutrophils_Wang$time, levels = c(0,3,14))
VlnPlot(Neutrophils_Wang, features = "BoneMarrow1", group.by = "time",
        pt.size = 0.08) + ggtitle("Bone Marrow Proximity Score")+
  xlab("Time (days)")+
  ylab("Wang: BM Proximity Scores")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))

bmDat <- data.frame("time" = Neutrophils_Wang$time, "bm" = Neutrophils_Wang$BoneMarrow1)
bmDat %>% group_by(time) %>% summarise(bm = mean(bm))

neuMaturationGenes <- read_csv("~/Documents/SingleCell/neuMaturationGenes.csv")
maturationGenes <- list(neuMaturationGenes[[1]])

Neutrophils_Wang <- AddModuleScore(Neutrophils_Wang, features = maturationGenes, search = T, assay = "RNA",
                                   name = "Maturation")

VlnPlot(Neutrophils_Wang, features = "Maturation1", group.by = "time") + ggtitle("Maturation Wang")

saveRDS(Neutrophils_Wang, file = "~/Documents/SingleCell/WangPaper/WangDat.rds")


#Neutrotime
earlyGenes <- list(c("Lyz2", "Retnlg", "Lgals3", "Mmp8", "Ifitm6", "Ly6g", "Pfn1", "Cybb", "Cd177", "Mgst1", "Pglyrp1",
                     "Ltf", "Prdx5", "Lcn2", "Ngp", "Adpgk", "Anxa1", "Arhgdib", "Ly6c2", "Camp", "Dstn", "Chil3",
                     "Wfdc21", "Serpinb1a"))

lateGenes <- list(c("Il1b", "Ccl6", "Fth1", "H2-D1", "Dusp1", "Junb", "Ifitm1", "Fxyd5", "Ifitm2", "Malat1", "Btg1",
                    "Tyrobp", "Jund", "Ftl1", "Srgn", "Csf3r", "Wfdc17", "Rps27", "Msrb1", "Fau", "Rps9"))

Neutrophils_Wang <- AddModuleScore(Neutrophils_Wang, features = earlyGenes, search = T, assay = "RNA",
                                   name = "neutrotimeEarly") 

Neutrophils_Wang <- AddModuleScore(Neutrophils_Wang, features = lateGenes, search = T, assay = "RNA",
                                   name = "neutrotimeLate")

VlnPlot(Neutrophils_Wang, features = "neutrotimeEarly1", group.by = "time",
        pt.size = 0.1)+ ggtitle("Early Neutrotime Scores")+
  xlab("Time (days)")

VlnPlot(Neutrophils_Wang, features = "neutrotimeLate1", group.by = "time",
        pt.size = 0.08)+ ggtitle("Late neutrotime")+
  xlab("Time (days)")

data.frame(late = Neutrophils_Wang$neutrotimeLate1, time = Neutrophils_Wang$time) %>% group_by(time) %>% 
  summarise(meanLate = mean(late))

data.frame(early = Neutrophils_Wang$neutrotimeEarly1, time = Neutrophils_Wang$time) %>% group_by(time) %>% 
  summarise(meanLate = mean(early))


#stats tests for the above plots

#Early neutrotime
matComp <- data.frame(day = Neutrophils_Wang@meta.data$time, mat = Neutrophils_Wang@meta.data$neutrotimeEarly1) 

w <- kruskal.test(mat ~ day, data = matComp)
pairwise.wilcox.test(matComp$mat, matComp$day, p.adjust.method = "bonferroni")


#Late neutrotime
matComp <- data.frame(day = Neutrophils_Wang@meta.data$time, mat = Neutrophils_Wang@meta.data$neutrotimeLate1) 

w <- kruskal.test(mat ~ day, data = matComp)
pairwise.wilcox.test(matComp$mat, matComp$day, p.adjust.method = "bonferroni")

#Bone marrow
matComp <- data.frame(day = Neutrophils_Wang@meta.data$time, mat = Neutrophils_Wang@meta.data$BoneMarrow1) 

w <- kruskal.test(mat ~ day, data = matComp)
pairwise.wilcox.test(matComp$mat, matComp$day, p.adjust.method = "bonferroni")



#Classify by sex
#Believe none of the "Y genes" from cellxy are in this dataset

counts <- logcounts(sce)
ann <- AnnotationDbi::select(org.Mm.eg.db, keys=rownames(sce),
                             columns=c("ENSEMBL","SYMBOL"), keytype="SYMBOL")
m <- match(rownames(counts), ann$SYMBOL)
rownames(counts) <- ann$SYMBOL[m]
sex <- classifySex(counts, genome="Mm")

table(sex$prediction)



#Add volcano plots 
Idents(Neutrophils_Wang) <- 'time'
Neu.Volmarkers <- FindAllMarkers(Neutrophils_Wang, only.pos = FALSE, 
                                 min.pct = 0.25, logfc.threshold = 0.15)

Neu.Volmarkers$pct.diff <- Neu.Volmarkers$pct.1 - Neu.Volmarkers$pct.2

Uninjured <- Neu.Volmarkers[Neu.Volmarkers$cluster == 0,]

Uninjured <- Uninjured %>% 
  mutate(significant = ifelse(abs(avg_log2FC) > 1 & p_val_adj < 0.01, "Significant", "Not significant"))

top_genes <- Uninjured %>% 
  filter(p_val_adj < 0.01, abs(avg_log2FC) > 1) %>% 
  arrange(p_val_adj) %>% 
  slice_head(n = 30)
top_genes <- top_genes$gene


pdf(file = "~/Documents/SingleCellPaper/uninjuredVol.pdf", height = 6, width = 8)

EnhancedVolcano(Uninjured,
                lab = Uninjured$gene,
                selectLab = top_genes,
                x = "avg_log2FC",
                y = "p_val_adj",
                title = 'Wang: DEGs for Uninjured',
                legendPosition = 'None',
                subtitle = "")

#Add dot plots

noNA <- subset(seuObj.integrated, pruned_labels != "NA")
Idents(noNA) <- "pruned_labels"
DotPlot(noNA, features = c("S100a9", "S100a8", "Ly6g", "Ltf", "Mmp9"))+
  RotatedAxis()