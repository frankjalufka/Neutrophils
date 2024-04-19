library(Seurat)
library(ggplot2)
library(dplyr)
library(destiny)
library(cowplot)
library(monocle3)
#library(SeuratWrappers)
library(gt)


#options(buildtools.check = function(action) TRUE )

sci <- readRDS(file = '~/Desktop/Thesis/SingleCell/LeePaper/Data/sci.rds')
#Neutrophils_Lee <- readRDS("/Users/matthewbradley/Documents/SingleCell/LeePaper/LeeNeu.rds")
########Find DE genes in Neutrophils_Lee. Wilcox is default, nonparametric t-test

Idents(sci) <- 'celltype'

seurat_markers <- FindAllMarkers(object = sci,
                                 assay = 'RNA',
                                 logfc.threshold = 1,
                                 only.pos = TRUE)



seurat_markers = seurat_markers[seurat_markers$cluster == "Neutrophil",]
seurat_markers$pct.diff = seurat_markers$pct.1 - seurat_markers$pct.2
topMarks <- head(seurat_markers[order(seurat_markers$pct.diff, decreasing = TRUE),], n = 240)

neuMarkers <- FindMarkers(sci, ident.1 = "Neutrophil", test.use = "wilcox",
                          only.pos = TRUE)
neuMarkers$pct.diff = neuMarkers$pct.1 - neuMarkers$pct.2
topneu <- head(neuMarkers[order(neuMarkers$pct.diff, decreasing = TRUE),], n = 240)

head(seurat_markers)
head(neuMarkers)


write.table(rownames(topneu), file="~/Documents/topneu", sep="\n", row.names = FALSE, quote=FALSE)
write.table(rownames(topMarks), file="~/Documents/topMarks", sep="\n", row.names = FALSE, quote=FALSE)

# What we see if findallmarkers is much stricter, not sure what parameter is causing this though
both = intersect(rownames(neuMarkers), rownames(seurat_markers))
diff = setdiff(rownames(neuMarkers), rownames(seurat_markers))

unique(rownames(topneu)[! rownames(topneu) %in% rownames(seurat_markers)])
unique(rownames(seurat_markers)[! rownames(seurat_markers) %in% rownames(topneu)])


######UMAP just Neutrophil cells

#subset data to just Neutrophils_Lee
Neutrophils_Lee = subset(sci, celltype == "Neutrophil")


Neutrophils_Lee$time <- factor(x = Neutrophils_Lee$time, levels = c('Uninjured','1dpi','3dpi','7dpi'), labels = c("0", "1", "3", "7"))
Neutrophils_Lee$sample_id <- factor(Neutrophils_Lee$sample_id, levels = c('uninj_sample1', 'uninj_sample2', 'uninj_sample3',
                                                                          '1dpi_sample1', '1dpi_sample2', '1dpi_sample3',
                                                                          '3dpi_sample1', '3dpi_sample2', '7dpi_sample1',
                                                                          '7dpi_sample2'))


DefaultAssay(Neutrophils_Lee) <- 'integrated'

Neutrophils_Lee <- RunPCA(Neutrophils_Lee, npcs = 50)
ElbowPlot(Neutrophils_Lee, ndims = 50) +
  geom_vline(mapping = aes(xintercept = 15), linewidth = 1, color = 'red', linetype = 'dashed') +
  labs(title = paste0('PCA of full SCI dataset (', ncol(Neutrophils_Lee), ' cells)'),
       subtitle = 'PCs up to dashed line used for downstream analysis')



npcs <- 1:15
Neutrophils_Lee <- FindNeighbors(Neutrophils_Lee, dims = npcs)
Neutrophils_Lee <- FindClusters(Neutrophils_Lee, resolution = 0.8) # default resolution
Neutrophils_Lee <- RunUMAP(Neutrophils_Lee, dims = npcs, min.dist = 0.3, n.neighbors = 30L, umap.method = 'uwot')
bad_cells <- colnames(Neutrophils_Lee)[Neutrophils_Lee$seurat_clusters == 28] # No

Neutrophils_Lee$default_cluster <- Neutrophils_Lee$seurat_clusters


umap_theme <- theme(panel.background = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = 'black'),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_text(size = 12, color = 'black'),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.text = element_text(size = 12, color = 'black'))


tmp <- FetchData(object = Neutrophils_Lee, vars = c('UMAP_1','UMAP_2','default_cluster')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = default_cluster), size = 1, alpha = 0.5) +
  umap_theme +
  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))+
  ggtitle("Neutrophil UMAP")
umap_defaultcluster <- LabelClusters(plot = tmp, id = 'default_cluster', size = 6, repel = FALSE)




#Choose color for each time
time_cols <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(time_cols) <- c('Uninjured','1dpi','3dpi','7dpi')

### Plot I was trying to recreate
umap_time <- FetchData(object = Neutrophils_Lee, vars = c('UMAP_1','UMAP_2','time')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = time), size = 1, alpha = 0.9) +
  scale_color_manual(values = time_cols) +
  umap_theme +
  theme(legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Time after\ninjury', override.aes = list(size = 7, alpha = 1)))+
  ggtitle("Neutrophil UMAP, labeled by time")

DimPlot(Neutrophils_Lee, split.by = "time", group.by = "time")
DimPlot(Neutrophils_Lee, group.by = "time")

#Values for original sample UMAP
sample_cols <- RColorBrewer::brewer.pal(n = 10, name = 'Spectral')
sample_cols <- c('firebrick','dodgerblue','goldenrod','firebrick','dodgerblue','goldenrod','firebrick','dodgerblue','firebrick','dodgerblue')
names(sample_cols) <- c('uninj_sample1', 'uninj_sample2', 'uninj_sample3','1dpi_sample1', '1dpi_sample2', '1dpi_sample3', '3dpi_sample1', '3dpi_sample2', '7dpi_sample1', '7dpi_sample2')

# Original sample UMAP
umap_sample <- FetchData(object = Neutrophils_Lee, vars = c('UMAP_1','UMAP_2','sample_id', 'time')) %>%
  .[sample(1:nrow(.), size = nrow(.)),] %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(color = sample_id), size = 0.2, alpha = 0.5) +
  facet_wrap(. ~ time, nrow = 2) +
  scale_color_manual(values = sample_cols) +
  umap_theme +
  theme(strip.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 12, color = 'black')) +
  guides(color = guide_legend(title = 'Sample ID\n(colored by\nreplicate)', override.aes = list(size = 7, alpha = 1)))

# 1 dpi most, uninjured least
sample_nums <- table(Neutrophils_Lee@meta.data[['time']])



#Find markers differentiating neutrophil groups
neu_markers <- FindAllMarkers(object = Neutrophils_Lee,
                              assay = 'RNA',
                              only.pos = TRUE)

neu_markers
neu_markers$pct.diff = neu_markers$pct.1 - neu_markers$pct.2

Uninjured <- neu_markers[neu_markers$cluster == 0,]
UninjuredGenes = unlist(split(Uninjured$p_val_adj, rownames(Uninjured)))
UninjuredGenes <- UninjuredGenes[UninjuredGenes <0.05]
write.csv(c(names(UninjuredGenes), UninjuredGenes), "Uninjured", row.names = F, quote = F)


OneGenes <- neu_markers[neu_markers$cluster == 1,]
OneGenes = unlist(split(OneGenes$p_val_adj, rownames(OneGenes)))
OneGenes <- OneGenes[OneGenes <0.05]
write.csv(c(names(OneGenes), OneGenes), "OneGenes", row.names = F, quote = F)

ThreeGenes <- neu_markers[neu_markers$cluster == 3,]
ThreeGenes = unlist(split(ThreeGenes$p_val_adj, rownames(ThreeGenes)))
ThreeGenes <- ThreeGenes[ThreeGenes <0.05]
write.csv(c(names(ThreeGenes), ThreeGenes), "ThreeGenes", row.names = F, quote = F)

SevenGenes <- neu_markers[neu_markers$cluster == 7,]
SevenGenes = unlist(split(SevenGenes$p_val_adj, rownames(SevenGenes)))
SevenGenes <- SevenGenes[SevenGenes <0.05]
write.csv(c(names(SevenGenes), SevenGenes), "SevenGenes", row.names = F, quote = F)

neu_markers_top <- neu_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = pct.diff)

View(neu_markers_top[,c(7,1,2,3,4,5,6,8)])

#Nice table of markers
neu_markers_top$cluster = as.numeric(as.character(neu_markers_top$cluster))
neu_markers_top %>% arrange(desc(pct.diff)) %>% 
  select(c("p_val_adj", "cluster", "gene", "pct.1", "pct.2", "pct.diff")) %>% 
  head(n = 10) %>% 
  gt(groupname_col = NULL) %>% 
  tab_header(title = "Differentially expressed genes", subtitle = "Neutrophils_Lee")


VlnPlot(Neutrophils_Lee, features = c("Camp"), assay = "RNA")

#Find markers separating Neutrophils_Lee by day

Idents(Neutrophils_Lee) <- 'time'

neuTime <- FindAllMarkers(object = Neutrophils_Lee,
                          assay = 'RNA',
                          only.pos = TRUE)

neuTime$pct.diff = neuTime$pct.1 - neuTime$pct.2
View(neuTime[order(neuTime$pct.diff, decreasing = TRUE),])

#Split by day
uNeu = Neutrophils_Lee[,Neutrophils_Lee$time == "Uninjured"]
oNeu = Neutrophils_Lee[,Neutrophils_Lee$time == "1dpi"]
tNeu = Neutrophils_Lee[,Neutrophils_Lee$time == "3dpi"]
sNeu = Neutrophils_Lee[,Neutrophils_Lee$time == "7dpi"]

#uGenes <- uNeu@assays$RNA@counts[rowSums(uNeu@assays$RNA@counts) > 40,]

uGenes <- FindMarkers(Neutrophils_Lee, ident.1 = "Uninjured", test.use = "wilcox",
                      only.pos = TRUE)
uGenes$pct.diff = (uGenes$pct.1 - uGenes$pct.2)
uGenes = subset(uGenes, p_val_adj < 0.05)
uGenes <- rownames(uGenes)

#oGenes <- oNeu@assays$RNA@counts[rowSums(oNeu@assays$RNA@counts) > 400,]

oGenes <- FindMarkers(Neutrophils_Lee, ident.1 = "1dpi", test.use = "wilcox",
                      only.pos = TRUE)
oGenes$pct.diff = (oGenes$pct.1 - oGenes$pct.2)
oGenes = subset(oGenes, p_val_adj < 0.05)
oGenes <- rownames(oGenes)

#tGenes <- tNeu@assays$RNA@counts[rowSums(tNeu@assays$RNA@counts) > 40,]

tGenes <- FindMarkers(Neutrophils_Lee, ident.1 = "3dpi", test.use = "wilcox",
                      only.pos = TRUE)
tGenes$pct.diff = (tGenes$pct.1 - tGenes$pct.2)
tGenes = subset(tGenes, p_val_adj < 0.05)
tGenes <- rownames(tGenes)



#sGenes <- sNeu@assays$RNA@counts[rowSums(sNeu@assays$RNA@counts) > 40,]

sGenes <- FindMarkers(Neutrophils_Lee, ident.1 = "7dpi", test.use = "wilcox",
                      only.pos = TRUE)
sGenes$pct.diff = (sGenes$pct.1 - sGenes$pct.2)
sGenes = subset(sGenes, p_val_adj < 0.05)
sGenes <- rownames(sGenes)

(c(length(uGenes), length(oGenes), length(tGenes), length(sGenes)))

write.table((oGenes), file="~/Documents/oGenes", sep="\n", row.names = FALSE, quote=FALSE)
write.table((uGenes), file="~/Documents/uGenes", sep="\n", row.names = FALSE, quote=FALSE)
write.table((tGenes), file="~/Documents/tGenes", sep="\n", row.names = FALSE, quote=FALSE)
write.table((sGenes), file="~/Documents/sGenes", sep="\n", row.names = FALSE, quote=FALSE)


#Get gene data for Neutrophils_Lee
GetAssayData(Neutrophils_Lee, slot = "counts")

#Plot log of gene occurences
hist(log10(rowSums(Neutrophils_Lee@assays$RNA@counts)+1), main = "Number of transcripts",
     xlab = "log10 counts (+1)")

hist(log10(rowSums(Neutrophils_Lee@assays$RNA@counts)[rowSums(Neutrophils_Lee@assays$RNA@counts)>0]),
     main = "counts of each gene", xlab = "log10 counts")


#Get gene names for gene ontology
presentGenes <- Neutrophils_Lee@assays$RNA@counts[rowSums(Neutrophils_Lee@assays$RNA@counts) > 0,]
presentGenes <- rownames(presentGenes)
length(presentGenes)
write.table(presentGenes, file="~/Documents/presentGenes", sep="\n", row.names = FALSE, quote=FALSE)

write.table(rownames(neuMarkers), file="~/Documents/neuMarkers", sep="\n", row.names = FALSE, quote=FALSE)
nrow(neuMarkers)

#Dot plot code for just Neutrophils_Lee
topneu <- head(neu_markers[order(neu_markers$pct.diff, decreasing = TRUE),], n = 49)
de_markers <- list("Neutrophil" = rownames(topneu))
de_markers <- append(de_markers$Neutrophil, "Ly6g")

DefaultAssay(sci) <- 'SCT'
avg_exp <- ScaleData(sci[['SCT']]@data, features = unlist(de_markers, use.names = FALSE))
avg_exp <- cbind(t(avg_exp), sci@meta.data[,c('celltype','time')]) %>%
  reshape2::melt(id.vars = c('celltype','time')) %>%
  group_by(celltype, time, variable) %>%
  summarise(avg.exp = mean(value))
pct_exp <- sci[['RNA']]@counts[which(rownames(sci[['RNA']]@counts) %in% de_markers) ,]
pct_exp = as.matrix(pct_exp)
pct_exp <- cbind(t(pct_exp), sci@meta.data[,c('celltype','time')]) %>%
  reshape2::melt(id.vars = c('celltype','time')) %>%
  group_by(celltype, time, variable) %>%
  summarise(pct.exp = mean(value > 0) * 100)
expr_colors <- colorRampPalette(colors = c('grey85', 'red3'))(100)

de_markers_dotplot <- merge(avg_exp, pct_exp) %>%
  filter(celltype == "Neutrophil") %>%
  mutate(time = factor(time, levels = rev(levels(time)))) %>%
  ggplot(mapping = aes(x = variable, y = time)) +
  geom_point(mapping = aes(size = pct.exp, fill = avg.exp), color = 'black', pch = 21) +
  scale_size(range = c(0,7), limits = c(0,100)) +
  scale_fill_gradientn(colors = expr_colors,
                       limits = c(NA, 3.5),
                       breaks = seq(-5, 10, 1),
                       na.value = expr_colors[length(expr_colors)]) +
  scale_y_discrete(position = 'right') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 16, color = 'black'),
        strip.placement = 'outside',
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 16, color = 'black', hjust = 0),
        legend.title = element_text(size = 16, angle = 90, color = 'black', hjust = 0.5),
        legend.margin = margin(5,0,10,0),
        legend.position = 'right',
        legend.spacing.x = unit(x = 2, units = 'mm'),
        panel.border = element_rect(fill = NA, size = 1),
        panel.background = element_rect(fill = NA)) +
  guides(fill = guide_colorbar(title = 'z-score',
                               barwidth = 1.25,
                               frame.colour = 'black',
                               frame.linewidth = 1.25,
                               ticks.colour = 'black',
                               ticks.linewidth = 1.25,
                               title.position = 'left'),
         size = guide_legend(title = '% expression',
                             override.aes = list(fill = 'black'),
                             title.position = 'left'))

ggsave(filename = '~/Documents/celltype_markers_dotplot.tiff',
       plot = de_markers_dotplot, device = 'tiff', height = 10, width = 16.8)


#Heatmap of gene, ordered by time

gene1 = "Camp"
testValues = as.data.frame(Neutrophils_Lee@assays$RNA@counts[gene1,])
testValues$Names = rownames(testValues)
colnames(testValues)[[1]] = "Count"
testValues$time = factor(Neutrophils_Lee@meta.data$time, ordered = T, levels = c("Uninjured", "1dpi", "3dpi", "7dpi"))
testValues <- testValues[order(testValues$time),]
rownames(testValues) = seq(1:nrow(testValues))
testValues$index = seq(1, nrow(testValues))
#testValues <- testValues[50:75,]

# ggplot(data = testValues, aes(x = Names, y = gene1, fill = Count))+
#   geom_tile()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   facet_grid(.~time)+
#   scale_fill_gradient(low = "#0481B0",
#                       high = "#FF6900",
#                       guide = "colorbar")

ggplot(data = testValues, aes(x = Names, y = gene1, fill = Count))+
  geom_tile()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_grid(.~time)+
  scale_fill_gradient(low = "white",
                      high = "black",
                      guide = "colorbar")


#Violin plot multiple markers
violFun <- function(gene){
  if(gene %in% rownames(Neutrophils_Lee@assays$RNA@counts)){
    VlnPlot(Neutrophils_Lee, features = gene, assay = "RNA")+
      theme(legend.position = 'none')
  }
}

timePlusList = c("Malat1", "Tyrobp", "Ccl6", "Il1b", "Fth1", "Msrb1", "Btg1", "Fth1")
timeMinusList = c("Ngp", "Camp", "Ltf", "Chil3", "lfitm6", "Cybb", "Lcn2", "Mmp8")
newMinus = c("Mmp9", "Itgam", "Fcn1", "Camp", "Cybb", "Cst3", "Csf3r")

plusPlots <- lapply(timePlusList, violFun)
plusPlots <- plusPlots[lengths(plusPlots) != 0]
minusPlots <- lapply(newMinus, violFun)
minusPlots <- minusPlots[lengths(minusPlots) != 0]

cowPlus <- cowplot::plot_grid(plotlist = plusPlots,  align = "hv")
cowMinus <- cowplot::plot_grid(plotlist = minusPlots, align = "hv")

ggsave(filename = "~/Documents/cowPlus.pdf", plot = cowPlus, device = 'pdf', height = 15, width = 16.8)
ggsave(filename = "~/Documents/cowMinus.pdf", plot = cowMinus, device = 'pdf', height = 15, width = 16.8)

#diffusion maps

neu_sce <- as.SingleCellExperiment(Neutrophils_Lee, assay = 'RNA')
#dm <- DiffusionMap(neu_sce, verbose = TRUE)


#Janky Expression Scores
jankExp <- function(gene, ageGroup){
  gene <- as.character(gene)
  totCells <- ncol(ageGroup)
  totExp <- sum(ageGroup@assays$RNA@data[gene,])
  expLevel <- totExp / totCells
  return(expLevel)
}

#Genes from neutrotime paper
early <- c("Lyz2", "Retnlg", "Lgals3", "Mmp8", "Ifitm6", "Ly6g", "Pfn1", "Cybb", "Cd177", "Mgst1", "Pglyrp1",
           "Ltf", "Prdx5", "Lcn2", "Ngp", "Adpgk", "Anxa1", "Arhgdib", "Ly6c2", "Camp", "Dstn", "Chil3",
           "Wfdc21", "Serpinb1a")

late <- c("Il1b", "Ccl6", "Fth1", "H2-D1", "Dusp1", "Junb", "Ifitm1", "Fxyd5", "Ifitm2", "Malat1", "Btg1",
          "Tyrobp", "Jund", "Ftl1", "Srgn", "Csf3r", "Wfdc17", "Rps27", "Msrb1", "Fau", "Rps9")

late %in% rownames(Neutrophils_Lee@assays$RNA@data)

ages <- c(uNeu, oNeu, tNeu, sNeu)
earlyRes <- unlist(lapply(ages, jankExp, gene = early))
lateRes <- unlist(lapply(ages, jankExp, gene = late))

#Plot early and late genes
agesChar <- c("uNeu", "oNeu", "tNeu", 'sNeu')
useAges <- factor((rep(agesChar, times =1, each = 2)), levels = c("uNeu", "oNeu", "tNeu", "sNeu"))
timePlotDat <- data.frame(Ages = useAges, times = rep(c("early", "late"), 4), values = c(earlyRes[1], lateRes[1],
                                                                                         earlyRes[2], lateRes[2],
                                                                                         earlyRes[3], lateRes[3],
                                                                                         earlyRes[4], lateRes[4]))

ggplot(timePlotDat, aes(x = Ages, y = values, fill = times))+
  geom_bar(position = "dodge", stat = "identity")+
  ggtitle("Gene time signatures by timepoint")+
  ylab("Gene Expression")


#Plot differences
diff <- lateRes - earlyRes
diffDat <- data.frame(ages = factor(agesChar, levels = c("uNeu", "oNeu", "tNeu", "sNeu")), diffs = diff)
ggplot(data = diffDat, aes(x = ages, y = diffs, fill = diffs))+
  geom_bar(stat = "identity")+
  ylab("Late - Early")+
  ggtitle("Gene time signatures by timepoint")



#Larger early and late (50 genes each)

#not included: H2afz
earlyLong <- c("Lyz2", "Ngp", "Wfdc21", "Camp", "Ifitm6", "Lcn2", "Ltf", "Arhgdib", "Mmp8", "Anxa1", "Prdx5", "Pglyrp1",
               "Cybb", "Ly6c2", "Chil3", "Dstn", "Cd177", "Pfn1", "Lgals3", "Ly6g", "Retnlg", "Adpgk",
               "Mgst1", "mt-Co1", "Serpinb1a", "Tkt", "Aldh2", "Capg", "Lamtor4", "Hmgn2", "Cd24a",
               "AA467197", "Cebpe", "Scp2", "Stmn1", "Plaur", "mt-Co3", "Tecr", "Sri", "Mcemp1", "Mmp9",
               "Tmsb10", "mt-Atp6", "Cd63", "mt-Co2", "Ltb4r1", "Chil1", "Limd2", "Ube2c")

earlyLong %in% rownames(Neutrophils_Lee@assays$RNA@data)
earlyLongRes <- unlist(lapply(ages, jankExp, gene = earlyLong))




features = c("Lyz2", "Camp", "Cybb", "Ngp") 

# featPlot <- RidgePlot(Neutrophils_Lee, features = features, ncol = 2)
# 
# ggsave(filename = "~/Documents/featPlot.jpg", plot = featPlot, device = 'jpeg', height = 15, width = 16.8)

DotPlot(Neutrophils_Lee, features = features) + RotatedAxis()


#Monocle integration

cds <- as.cell_data_set(Neutrophils_Lee)
cds <- cluster_cells(cds, resolution=1e-3)

list.cluster <- Neutrophils_Lee@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster


p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

cds <- learn_graph(cds, use_partition = F, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=T,
           group_label_size = 5,
           label_leaves=FALSE,
           label_branch_points=T)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 0]))

pseudotimePlot <- plot_cells(cds,
                             color_cells_by = "pseudotime",
                             group_cells_by = "cluster",
                             label_cell_groups = FALSE,
                             label_groups_by_cluster=FALSE,
                             label_leaves=FALSE,
                             label_branch_points=FALSE,
                             label_roots = FALSE,
                             trajectory_graph_color = "black")

p1 = umap_time
p2 = pseudotimePlot
combPlot <- wrap_plots(p1, p2)
ggsave(filename = "~/Documents/combPlot.pdf", plot = combPlot, device = 'pdf', height = 8, width = 11)


#Bone marrow proximity
#From supplement methods of Dynamics of Cardiac Neutrophil Diversity in Murine Myocardial Infarction

bmGenes <- list(c("Mmp8", "Ifitm6", "Mmp25", 
                  "Retnlg", "Lcn2", "Olfm4", "Chil3", 
                  "Itgb2", "Fpr1", "Ltf", "Camp", "Lyz2", 
                  "Lyz1"))

Neutrophils_Lee <- AddModuleScore(Neutrophils_Lee, features = bmGenes, search = T, assay = "RNA",
                                  name = "BoneMarrow")

Idents(Neutrophils_Lee) <- "time"
VlnPlot(Neutrophils_Lee, features = "BoneMarrow1") +ggtitle("Bone Marrow Lee")

bmComp <- data.frame(day = Neutrophils_Lee@meta.data$time, bm = Neutrophils_Lee@meta.data$BoneMarrow1)
kw <- kruskal.test(bm ~ day, data = bmComp)
pairwise.wilcox.test(bmComp$bm, bmComp$day, p.adjust.method = "bonferroni")


#Maturation Score from "Single-cell transcriptome profiling 
#reveals neutrophil heterogeneity in homeostasis and infection" Supp table 4
neuMaturationGenes <- read_csv("~/Documents/SingleCell/neuMaturationGenes.csv")
maturationGenes <- list(neuMaturationGenes[[1]])

Neutrophils_Lee <- AddModuleScore(Neutrophils_Lee, features = maturationGenes, search = T, assay = "RNA",
                                  name = "Maturation")


VlnPlot(Neutrophils_Lee, features = "Maturation1") + ggtitle("Maturation Scores Lee")
matComp <- data.frame(day = Neutrophils_Lee@meta.data$time, mat = Neutrophils_Lee@meta.data$Maturation1) 
#%>% group_by(day) %>%  summarise(mean = mean(mat))

w <- kruskal.test(mat ~ day, data = matComp)
pairwise.wilcox.test(matComp$mat, matComp$day, p.adjust.method = "bonferroni")

#Neutrotime early and late based on fig 2
earlyGenes <- list(c("Lyz2", "Retnlg", "Lgals3", "Mmp8", "Ifitm6", "Ly6g", "Pfn1", "Cybb", "Cd177", "Mgst1", "Pglyrp1",
                     "Ltf", "Prdx5", "Lcn2", "Ngp", "Adpgk", "Anxa1", "Arhgdib", "Ly6c2", "Camp", "Dstn", "Chil3",
                     "Wfdc21", "Serpinb1a"))

lateGenes <- list(c("Il1b", "Ccl6", "Fth1", "H2-D1", "Dusp1", "Junb", "Ifitm1", "Fxyd5", "Ifitm2", "Malat1", "Btg1",
                    "Tyrobp", "Jund", "Ftl1", "Srgn", "Csf3r", "Wfdc17", "Rps27", "Msrb1", "Fau", "Rps9"))

Neutrophils_Lee <- AddModuleScore(Neutrophils_Lee, features = earlyGenes, search = T, assay = "RNA",
                                  name = "neutrotimeEarly") 

Neutrophils_Lee <- AddModuleScore(Neutrophils_Lee, features = lateGenes, search = T, assay = "RNA",
                                  name = "neutrotimeLate")

VlnPlot(Neutrophils_Lee, features = "neutrotimeEarly1", group.by = "time")+ ggtitle("Early neutrotime")
VlnPlot(Neutrophils_Lee, features = "neutrotimeLate1", group.by = "time")+ ggtitle("Late neutrotime")


#Gene ontology files


