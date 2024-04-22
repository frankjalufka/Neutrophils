library(RColorBrewer)
library(dplyr)
library(gt)
library(viridis)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dunn.test)
library(readr)

setwd("~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/")

#import datasets
Neutrophils_Bren <- readRDS("~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/Neutrophils_Bren.rds")
Neutrophils_Wang <- readRDS("~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/Neutrophils_Wang.rds")
Neutrophils_Lee <- readRDS("~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/Neutrophils_Lee.RDS")
#label datasets by paper
Neutrophils_Bren@meta.data$trial = "Brennan"
Neutrophils_Lee@meta.data$trial = "Lee"
Neutrophils_Wang@meta.data$trial = "Wang"

Neu_Integrated <- readRDS("~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/Neu_integrated.rds")

#Combine datasets into list for integration
Neu_Combined = list(c(Neutrophils_Bren,Neutrophils_Lee, Neutrophils_Wang ))[[1]]

#Integrate data from list
features <- SelectIntegrationFeatures(object.list = Neu_Combined)
anchors <- FindIntegrationAnchors(object.list = Neu_Combined)
Neu_Integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(Neu_Integrated) <- "integrated"

#Begin data processing
Neu_Integrated <- ScaleData(Neu_Integrated, verbose = FALSE)
Neu_Integrated <- RunPCA(Neu_Integrated, verbose = FALSE)

ElbowPlot(Neu_Integrated, ndims = 40)

Neu_Integrated <- FindNeighbors(Neu_Integrated, dims = 1:30)
Neu_Integrated <- FindClusters(Neu_Integrated)

#Plot clusters of integrated data
Idents(Neu_Integrated) <- 'seurat_clusters'
Neu_Integrated <- RunUMAP(Neu_Integrated, dims = 1:30)
DimPlot(Neu_Integrated, reduction = "umap", label = F, group.by = "seurat_clusters",
        split.by = "time") + NoLegend()

#Make labels consistent, should do this earlier for consistency
Neu_Integrated$time <- factor(Neu_Integrated$time, levels = c("0", "1", "3", "7", "14", "28"))


#Custom colors for upcoming plots
setwd("~/Desktop/GitHub/Neutrophils/Data")
tiff("CombinedNeutrophils.tif", units = "in", width = 6, height = 6, res = 300)
DimPlot(Neu_Integrated, group.by = "time", label = F, cols = c("red", "green", "blue", "gold", "purple", "aquamarine4")) +
  guides(color=guide_legend(ncol =1, override.aes = list(size=5))) +
  ggtitle("Combined Neutrophil UMAP")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))
dev.off()

saveRDS(Neu_Integrated, file = "~/Desktop/Thesis/SingleCell/NeutrophilFiles/Version_3.0/Neu_integrated.rds")

DimPlot(Neu_Integrated, split.by = "time", group.by = "time")
DimPlot(Neu_Integrated, split.by = "time", group.by = "trial")


tiff("~/Desktop/GitHub/Neutrophils/Data/CombinedUmapByTrial.tif", units = "in", width = 6, height = 6, res = 300)
DimPlot(Neu_Integrated, split.by = "trial", group.by = "trial") + NoLegend()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text=element_text(size=20),
        strip.text = element_text(size=20),
        plot.title = element_text(size=22)
  )+
  ggtitle("Neutrophils by Trial")
dev.off()

tiff("~/Desktop/GitHub/Neutrophils/Data/CombinedUmapByTime.tiff", units = "in", width = 8, height = 6, res = 300)
DimPlot(Neu_Integrated, split.by = "time", group.by = "time") + NoLegend()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text=element_text(size=20),
        strip.text = element_text(size=20),
        plot.title = element_text(size=28)
  )+
  ggtitle("Neutrophils by Time (dpi)")
dev.off()


#Find markers for each cluster
Idents(Neu_Integrated) <- 'time'
neuIntMarkers <- FindAllMarkers(Neu_Integrated, min.pct = 0.25, logfc.threshold = 0.25)
neuIntMarkers$pct.diff <- neuIntMarkers$pct.1 - neuIntMarkers$pct.2

subsetMarkers <- function(myCluster, markerList){
  markerList <- subset(markerList, cluster == myCluster)
  markerList <- subset(markerList, p_val_adj < 0.01)
  markerList %>% arrange(-avg_log2FC) %>% 
    head(n = 30)
}

topDEGByCluster <- lapply(unique(neuIntMarkers$cluster), subsetMarkers, markerList = neuIntMarkers)

tiff("~/Documents/SingleCell/IntegratedNeu/HeatMapByTime.tiff", units = "in", width = 10, height = 13, res = 300)
DoHeatmap(Neu_Integrated, features = unlist(lapply(topDEGByCluster, rownames))) 
dev.off()

DoHeatmap(Neu_Integrated, features = rownames(fourteenThreeMarks))
#+scale_fill_gradientn(colors = c("white", "#FC842B", "#A61100"))

quantile(Neu_Integrated$neutrotimeEarly1)

Neu_Integrated@meta.data <- Neu_Integrated@meta.data %>% mutate(
  earlyNeuCat = case_when(
    neutrotimeEarly1 < -0.2 ~ "Mature",
    between(neutrotimeEarly1, -0.2, 0.07) ~ "Medium Mature",
    between(neutrotimeEarly1, 0.07, 1.44) ~ "Medium Immature",
    neutrotimeEarly1> 1.44 ~ "Immature"
  )
)

Neu_Integrated$earlyNeuCat <- factor(Neu_Integrated$earlyNeuCat,
                                     levels = c("Immature", "Medium Immature", "Medium Mature", "Mature"))

tiff("~/Documents/SingleCell/IntegratedNeu/HeatMapByEarlyNeu.tiff", units = "in", width = 10, height = 13, res = 300)
DoHeatmap(Neu_Integrated, features = unlist(lapply(topDEGByCluster, rownames)), group.by = "earlyNeuCat") +
  guides(colour= "none")+
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 10))
dev.off()

tiff("~/Desktop/GitHub/Neutrophils/Data/EarlyNeutrotimeHeatmap.tif", units = "in", width = 10, height = 13, res = 300)
DoHeatmap(Neu_Integrated, features = earlyGenes[[1]], group.by = "time")+
  guides(colour= "none")+
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 10))
dev.off()

#Subset out only uninjured neutrophils
Uninjured_Neu <- subset(Neu_Integrated, time == "0")
DoHeatmap(Uninjured_Neu, features = VariableFeatures(Uninjured_Neu)[1:50], size = 4) + NoLegend()
x <- FindVariableFeatures(Uninjured_Neu)
VariableFeaturePlot(x)

Neu.markers <- FindAllMarkers(Neu_Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Neu.markers$pct.diff <- Neu.markers$pct.1 - Neu.markers$pct.2

Idents(Neu_Integrated) <- 'time'
Neu.Timemarkers <- FindAllMarkers(Neu_Integrated, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.15)

Neu.Timemarkers$pct.diff <- Neu.Timemarkers$pct.1 - Neu.Timemarkers$pct.2

#### GO ANALYSIS ###
Uninjured <- Neu.Timemarkers[Neu.Timemarkers$cluster == 0,]
UninjuredGenes = unlist(split(Uninjured$p_val_adj, rownames(Uninjured)))
UninjuredGenes <- UninjuredGenes[UninjuredGenes <0.05]
write.csv(c(names(UninjuredGenes), UninjuredGenes), "Uninjured.csv", row.names = F, quote = F)

#+ scale_fill_gradientn(colors = c("white", "#FC842B", "#A61100"))


topMarkers <- (neuIntMarkers %>%
                 arrange(p_val_adj, -avg_log2FC)  %>% 
                 group_by(cluster) %>% slice(1:5) %>% 
                 head(n = 30)) %>% select(gene)

tiff("dotplot.tiff", units = "in", width = 8, height = 5, res = 300)
DotPlot(Neu_Integrated, features = topMarkers$gene) + RotatedAxis()+
  RotatedAxis()+
  theme(legend.text = element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.size = unit(0.5, "cm"))
dev.off()

#Feature plots
topMarkers <- (neuIntMarkers %>%
                 arrange(p_val_adj, -pct.diff) %>% 
                 head(n = 30))

DefaultAssay(Neu_Integrated) <- "RNA"
FeaturePlot(Neu_Integrated, features = c("Ace"))


DefaultAssay(Neu_Integrated) <- "integrated"
#Heatmap


DoHeatmap(Neu_Integrated, features = VariableFeatures(Neu_Integrated)[1:50],  size = 4,
          angle = 90) + NoLegend()

DoHeatmap(Neu_Integrated, features = topMarkers$gene, cells = 1:10, size = 4,
          angle = 90) + NoLegend()


#BM Scores
bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", 
                  "Retnlg", "Lcn2", "Olfm4", "Chil3", 
                  "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", 
                  "Lyz1"))

Neu_Integrated <- AddModuleScore(Neu_Integrated, features = bmGenes, search = T, assay = "RNA",
                                 name = "BoneMarrow")

tiff("~/Desktop/GitHub/Neutrophils/Data/BMIntegrated.tif", units = "in", width = 5, height = 5, res = 300)

VlnPlot(Neu_Integrated, features = "BoneMarrow1", group.by = "time",
        pt.size = 0, cols = c("red", "green", "blue", "gold", "purple", "aquamarine4")) + 
  ggtitle("Bone Marrow Proximity Score")+
  xlab("Time (days)")+
  ylab("BM Proximity Feature Scores")+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))+
  ylim(c(-1,3.5))+
  geom_jitter(size = 0.01, alpha = 0.2) +
  stat_compare_means(label.y = 3.2, label.x = 1.5)

dev.off()

#BM Kruskal wallis
bmKruskal <- kruskal.test(BoneMarrow1 ~ time, data = Neu_Integrated@meta.data)
bmDunn <- dunn.test(x = Neu_Integrated@meta.data$BoneMarrow1, g = Neu_Integrated@meta.data$time)

as.data.frame(table(Neu_Integrated$time)) %>% gt()

#Early and late neutrotime
earlyGenes <- list(c("Lyz2", "Retnlg", "Lgals3", "Mmp8", "Ifitm6", "Ly6g", "Pfn1", "Cybb", "Cd177", "Mgst1", "Pglyrp1",
                     "Ltf", "Prdx5", "Lcn2", "Ngp", "Adpgk", "Anxa1", "Arhgdib", "Ly6c2", "Camp", "Dstn", "Chil3",
                     "Wfdc21", "Serpinb1a"))

lateGenes <- list(c("Il1b", "Ccl6", "Fth1", "H2-D1", "Dusp1", "Junb", "Ifitm1", "Fxyd5", "Ifitm2", "Malat1", "Btg1",
                    "Tyrobp", "Jund", "Ftl1", "Srgn", "Csf3r", "Wfdc17", "Rps27", "Msrb1", "Fau", "Rps9"))

Neu_Integrated <- AddModuleScore(Neu_Integrated, features = earlyGenes, search = T, assay = "RNA",
                                 name = "neutrotimeEarly") 

Neu_Integrated <- AddModuleScore(Neu_Integrated, features = lateGenes, search = T, assay = "RNA",
                                 name = "neutrotimeLate")

tiff("~/Desktop/GitHub/Neutrophils/Data/EarlyNeutro.tiff", units = "in", width = 5, height = 5, res = 300)

VlnPlot(Neu_Integrated, features = "neutrotimeEarly1", group.by = "time",
        pt.size = 0, cols = c("red", "green", "blue", "gold", "purple", "aquamarine4"))+ 
  ggtitle("Early Neutrotime Scores")+
  xlab("Time (days)")+
  ylim(-1,3.5)+
  geom_jitter(size = 0.01, alpha = 0.2) +
  stat_compare_means(label.y = 3.2, label.x = 1.5)

dev.off()

DoHeatmap(Neu_Integrated, features = earlyGenes[[1]], group.by = "neutrotimeEarly1")

tiff("~/Desktop/GitHub/Neutrophils/Data/LateNeutro.tiff", units = "in", width = 5, height = 5, res = 300)
VlnPlot(Neu_Integrated, features = "neutrotimeLate1", group.by = "time",
        pt.size = 0, cols = c("red", "green", "blue", "gold", "purple", "aquamarine4"))+ 
  ggtitle("Late neutrotime")+
  xlab("Time (days)")+
  ylim(-0.5,3.5)+
  geom_jitter(size = 0.01, alpha = 0.2) +
  stat_compare_means(label.y = 3.2, label.x = 1.5)
dev.off()

#Dotplots

tiff("~/Desktop/GitHub/Neutrophils/Data/DotplotImmature.tiff", units = "in", width = 5, height = 6, res = 300)
DotPlot(Neu_Integrated, features = c("Lyz2", "Ngp","Wfdc21", "Camp", "Ifitm6", "Lcn2",
                                     "Ltf", "Arhgdib", "Mmp8", "Anxa1"), assay = "RNA",
        group.by = "time") +
  coord_flip()+ 
  scale_colour_gradient2(low = "#B0FAFF", mid = "#FFD167", high = "#CB0000")+
  ggtitle("Early Neutrotime Genes")

dev.off()

tiff("~/Desktop/GitHub/Neutrophils/Data/DotplotMature.tiff", units = "in", width = 5, height = 6, res = 300)
DotPlot(Neu_Integrated, features = c("Hbb-bs", "Fxyd5", "Rps27", "Malat1",
                                     "Ifitm2", "Btg1", "Msrb1", "Srgn", "Tyrobp","Il1b"),
        assay = "RNA",
        group.by = "time" )+
  coord_flip()+ 
  scale_colour_gradient2(low = "#B0FAFF", mid = "#FFD167", high = "#CB0000")+
  ggtitle("Late Neutrotime Genes")

dev.off()
#Heatmap average expression
avgexp = AverageExpression(Neu_Integrated, assays = "RNA", return.seurat = T)
avgexp$orig.ident = c(0,1,3,7,14,28)
DoHeatmap(avgexp, features = earlyGenes[[1]], 
          assay = "RNA", slot = "data",
          group.by = "orig.ident",
          draw.lines = F)+ 
  scale_fill_viridis()+
  guides(colour=FALSE)

DoHeatmap(avgexp, features = lateGenes[[1]], 
          assay = "RNA", slot = "data",
          group.by = "orig.ident",
          draw.lines = F)+
  scale_fill_viridis()+
  guides(colour=FALSE)


#Write genes for gene ontology
geneMarkerList <- lapply(clusterNums <- seq(0:13), 
                        function(x) neuIntMarkers[neuIntMarkers$cluster == x,])

uninjuredMarkers <- neuIntMarkers[neuIntMarkers$cluster == 0,]


#writeMarkers <- function(geneList, name){
#  geneList = unlist(split(geneList$p_val_adj, rownames(geneList)))
#  geneList = geneList[geneList < 0.01]
#  write.csv(c(names(geneList)), name, row.names = F, quote = F)
  
#}

#writeMarkers(geneMarkerList[[0]], "uninjured")
write.csv(c(names(uninjuredMarkers), uninjuredMarkers), "UninjuredMarkers.csv", row.names =F, quote = F)

#Making nice table of neutrophils by trial and time
time_table <- Neu_Integrated@meta.data %>% group_by(time, trial) %>% summarise(total = n()) %>% 
  pivot_wider(names_from = time, values_from = total) %>% 
  replace_na(list(`0` = 0, `1` = 0, `3` = 0, `7` = 0, `14` = 0, `28` = 0))

time_table_gt <- gt(time_table, rowname_col  = "trial") %>% cols_label(`0` = "Uninjured") %>% 
  cols_width(
    everything() ~ px(80))

gtsave(time_table_gt, "gtTimes.pdf")

#Correlation heatmap

av.exp <- AverageExpression(Neu_Integrated, assays = "integrated", group.by = "time")$integrated
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- factor(rownames(cor.exp), levels = c('0','1','3','7','14','28'))
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c('0', '1', '3', '7', '14', '28'))
cor.df$y = factor(cor.df$y, levels = c('0','1','3','7','14','28'))
tiff("~/Desktop/GitHub/Neutrophils/Data/Correlation.tiff", units = "in", width = 5, height = 5, res = 300)
ggplot(cor.df, aes(x, reorder(y, x), fill = correlation)) +
  geom_tile()
dev.off()
