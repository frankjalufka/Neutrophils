library(CellChat)
library(NMF)
library(patchwork)
library(Seurat)
library(readr)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(ggalluvial)

#Read in Seurat objects
wang <- readRDS("WangDat.rds")
bren <- readRDS("BrenDat.rds")
lee <- readRDS("LeeDat.rds")

#subset uninjured data from whole datasets
wang.uninjured <- subset(wang, time == "0")
bren.uninjured <- subset(bren, time == "0")
lee.uninjured <- subset(lee, time == 'Uninjured')

lee.uninjured@meta.data$label = case_when(lee.uninjured$celltype == "Neutrophil" ~ "Neutrophils",
                                          lee.uninjured$celltype == "Macrophage" ~ "Macrophages",
                                          lee.uninjured$celltype == "Monocyte" ~ "Monocytes",
                                          lee.uninjured$celltype == "Fibroblast" ~ "Fibroblasts",
                                          lee.uninjured$celltype == "Endothelial" ~ "Endothelial cells",
                                          .default = lee.uninjured$celltype)
wang.uninjured@meta.data$label = case_when(wang.uninjured$pruned_labels == "Neutrophil" ~ "Neutrophils",
                                           wang.uninjured$pruned_labels == "Macrophage" ~ "Macrophages",
                                           wang.uninjured$pruned_labels == "Monocyte" ~ "Monocytes",
                                           wang.uninjured$pruned_labels == "Fibroblast" ~ "Fibroblasts",
                                           wang.uninjured$pruned_labels == "Endothelial" ~ "Endothelial cells",
                                           .default = wang.uninjured$pruned_labels)
bren.uninjured@meta.data$label = case_when(bren.uninjured$pruned_labels == "Neutrophil" ~ "Neutrophils",
                                           bren.uninjured$pruned_labels == "Macrophage" ~ "Macrophages",
                                           bren.uninjured$pruned_labels == "Monocyte" ~ "Monocytes",
                                           bren.uninjured$pruned_labels == "Fibroblast" ~ "Fibroblasts",
                                           bren.uninjured$pruned_labels == "Endothelial" ~ "Endothelial cells",
                                           .default = bren.uninjured$pruned_labels)

unique(wang.uninjured$label)
#merge seurat objects
uninjured.combined <- merge(wang.uninjured, y = c(bren.uninjured, lee.uninjured))

#Remove rows from dataset with NA in the pruned_labels field
unique(uninjured.combined$label)
uninjured.combined <- subset(uninjured.combined, subset = label != "NA")
unique(uninjured.combined$label)

saveRDS(uninjured.combined, file = "D:/Frank/uninjured_combined.rds")

#create cellchat object
cellchat <- createCellChat(object = uninjured.combined, group.by = "label", assay = "RNA")

#Check for NA values that will mess up downstream analysis
unique(cellchat@idents)

#create object that contains the ligand-receptor database for mice
CellChatDB <- CellChatDB.mouse

#choose all of CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB #simply use the default CellChatDB

#set the used database in the object
cellchat@DB <- CellChatDB.use

#subset the expression data of signaling genes to save computation cost
cellchat <- subsetData(cellchat) #this is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#compute the communication probability and infer cellular communication network
#longest step, requires a lot of RAM
cellchat <- computeCommunProb(cellchat, type = "triMean")

#Filter out cell-cell communication if fewer than 10 cells in a group
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "D:/Frank/cellchat_combined.rds")

cellchat <- readRDS("~/Desktop/Thesis/scRNA_Chapter/NeutrophilFiles/Version_3.0/cellchat_combined.rds")

groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = "Number of Interactions")

tiff("~/Desktop/GitHub/Neutrophils/Data/Neu_Target_circle.tif", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 targets.use = "Neutrophils",
                 top = 0.5,
                 margin = 0)

dev.off()

tiff("~/Desktop/GitHub/Neutrophils/Data/Neu_source_circle.tif", units = "in", width = 9, height = 9, res = 300)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 sources.use = "Neutrophils",
                 top = 0.5,
                 margin = 0)

dev.off()

tiff("~/Desktop/GitHub/Neutrophils/Data/Neu_source_bubble.tif", units = "in", width = 8, height = 8, res = 300)
netVisual_bubble(cellchat,
                 sources.use = "Neutrophils",
                 targets.use = c("Astrocyte", "B cells", "Dendridic", "Endothelial cells", "Ependymal", "Macrophages", "Microglia", "Neuron", "Neutrophils", "Oligodendrocyte", "Pericyte"),
                 remove.isolate = TRUE)
dev.off()

tiff("~/Desktop/GitHub/Neutrophils/Data/Neu_target_bubble.tif", units = "in", width = 8, height = 8, res = 300)
netVisual_bubble(cellchat,
                 sources.use = c("Dendridic", "Endothelial cells", "Ependymal", "Microglia", "Neuron", "Neutrophils", "NK cells", "Oligodendrocyte", "Pericyte"),
                 targets.use = "Neutrophils",
                 remove.isolate = TRUE)

dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.5, font.size = 10)

gg1 <- netAnalysis_signalingRole_scatter(cellchat)

pathways.show <- c("SELL")

netVisual_heatmap(cellchat, signaling = pathways.show)

dev.off()

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("TNF", "NOTCH"))
selectK(cellchat, pattern = "outgoing")

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
