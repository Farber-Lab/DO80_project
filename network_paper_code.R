############################################################
############################################################ 
############################################################

library(gridExtra)
library(Seurat)
library(speckle)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(readxl)
library(scCustomize)
library(Polychrome)
library(devtools)
library(scCustomize)
library(readxl)
library(tradeSeq)
library(slingshot)
library(igraph)
library(bnlearn)
library(graph)
library(qtl2)
library(ggsignif)
library(stringr)
library(limma)
library(GenomicRanges)
library(DESeq2)
library(preprocessCore)

############################################################
############################################################ 
############################################################

#make idents the cluster numbers
Idents(do80.combined_all_clusters) <- do80.combined_all_clusters@meta.data$seurat_clusters

#change numbers to annotations
do80.combined_all_clusters <- RenameIdents(object = do80.combined_all_clusters, 
                                           "0" = "OBP",
                                           "1" = "MF1",    
                                           "2" = "OB2",
                                           "3" = "LMP",
                                           "4" = "MALP",
                                           "5" = "LMP",
                                           "6" = "MF2",               
                                           "7" = "OB1",
                                           "8" = "MPC",
                                           "9" = "MALP",
                                           "10" = "Ocy",
                                           "11" = "LMP",
                                           "12" = "MO",
                                           "13" = "GC",
                                           "14" = "TC",
                                           "15" = "BC",
                                           "16" = "EC",
                                           "17" = "OC")

#store cell cluster annotations
do80.combined_all_clusters$celltype_annotation <- Idents(do80.combined_all_clusters)

#legend
do80.combined_all_clusters[["legend"]] <- recode(do80.combined_all_clusters$celltype_annotation, 
                                                 "MPC" = paste0("8   - Mesenchymal Progenitor                    - MPC ", "  (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MPC")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MPC")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "LMP" = paste0("4   - Late Mesenchymal Progenitor             - LMP ", "  (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "LMP")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "LMP")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "OBP" = paste0("1   - Osteoblast Progenitor                         - OBP ", " (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OBP")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OBP")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "OB1" = paste0("7   - Osteoblast Population 1                      - OB1 ", "  (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OB1")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OB1")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "OB2" = paste0("3   - Osteoblast Population 2                      - OB2 ", " (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OB2")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OB2")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "Ocy" = paste0("9   - Osteocyte-like cell                              - Ocy ", "    (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "Ocy")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "Ocy")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "MALP" = paste0("5   - Marrow Adipogenic Lineage Precursor - MALP ", "(", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MALP")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MALP")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "MF1" = paste0("2   - Macrophage Population 1                    - MF1 ", " (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MF1")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MF1")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "MF2" = paste0("6   - Macrophage Population 2                    - MF2 ", "  (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MF2")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MF2")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "MO" = paste0("10 - Monocyte                                           - MO ", "    (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MO")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "MO")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "TC" = paste0("12 - T-cell                                                  - TC ", "     (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "TC")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "TC")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "GC" = paste0("11 - Granulocyte                                       - GC ", "     (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "GC")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "GC")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "BC" = paste0("13 - B-cell                                                 - BC ", "     (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "BC")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "BC")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "EC" = paste0("14 - Endothelial cell                                   - EC ", "     (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "EC")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "EC")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"),
                                                 "OC" = paste0("15 - Osteoclast-like cell                             - OC ", "     (", paste(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OC")), ", ", paste(round(sum(do80.combined_all_clusters@meta.data$celltype_annotation == "OC")/ncol(do80.combined_all_clusters)*100, digits = 2)), "%", ")"))


############################################################          ############################################################
############################################################  Fig 1A  ############################################################
############################################################          ############################################################


DimPlot(do80.combined_all_clusters, reduction = 'umap', group.by = "legend", label = F, raster=F, pt.size = 0.25) +
  ggtitle("") + 
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  scale_colour_manual(name = "Cell Clusters",
                      values =  c('#9A6324',
                                  '#FDA4BA', 
                                  "#f58231", 
                                  "#f7c564",
                                  "#8762a6",
                                  "#82a653", 
                                  "#469990",
                                  '#4363d8',
                                  "#db4105", 
                                  '#52B2BF', 
                                  "#3e8cc6",
                                  '#7F7D9C', 
                                  '#601A35', 
                                  '#234F1E', 
                                  '#dcbeff'))

############################################################
############################################################
############################################################

#cell cluster annotations for cell cluster proportion analysis
do80.combined_all_clusters$celltype_annotation_props <- Idents(do80.combined_all_clusters)
do80.combined_all_clusters[["celltype_annotation_props"]] <- recode(do80.combined_all_clusters$celltype_annotation_props, 
                                                                    "MF1" = "Hem",
                                                                    "MF2" = "Hem",
                                                                    "MO" = "Hem",
                                                                    "GC"  = "Hem",
                                                                    "TC"  = "Hem",
                                                                    "BC"  = "Hem",
                                                                    "EC" = "Hem", 
                                                                    "OC" = "Hem")

#subset to include only mesenchmyal lineage cell types
do80.combined_mes_subset <- subset(x = do80.combined_all_clusters, subset = celltype_annotation_props != "Hem")
do80.combined_mes_subset$celltype_annotation <- Idents(do80.combined_mes_subset)

#aggregate expression from Bglap genes

bglap_sub <- mes_subset_markers[grep("^Bglap", mes_subset_markers$gene), ]
bglap_sub_OB1 <- mean((subset(bglap_sub, subset = cluster == "OB1"))$avg_log2FC)
bglap_sub_OB2 <- mean((subset(bglap_sub, subset = cluster == "OB2"))$avg_log2FC)
bglap_sub_OBP <- mean((subset(bglap_sub, subset = cluster == "OBP"))$avg_log2FC)
bglap_sub$avg_log2FC <- recode(bglap_sub$cluster, "OB1" = bglap_sub_OB1, "OB2" = bglap_sub_OB2, "OBP" = bglap_sub_OBP)
bglap_sub_final <- subset(bglap_sub, subset = gene == "Bglap")
mes_subset_markers <- mes_subset_markers[grep("^Bglap", mes_subset_markers$gene, invert = T), ]
mes_subset_markers <- rbind(mes_subset_markers, bglap_sub_final)

top5_markers <- Extract_Top_Markers(marker_dataframe = mes_subset_markers, 
                                    num_genes = 5,
                                    named_vector = FALSE,
                                    make_unique = TRUE)


############################################################          ############################################################
############################################################  Fig 1B  ############################################################
############################################################          ############################################################

Clustered_DotPlot(seurat_object = do80.combined_mes_subset, 
                  features = top5_markers, 
                  show_parent_dend_line = F,
                  colors_use_idents = c('#9A6324',"#f58231","#f7c564","#8762a6","#469990",'#4363d8',"#db4105"),
                  column_label_size = 25,
                  row_label_size = 18,
                  legend_label_size = 15,
                  legend_title_size = 12)

############################################################          ############################################################
############################################################  Fig 1C  ############################################################
############################################################          ############################################################

plotCellTypeProps(x = do80.combined_all_clusters, clusters = do80.combined_all_clusters$celltype_annotation_props, sample = do80.combined_all_clusters$mouse) +
  geom_bar(position='stack', stat='identity') +
  xlab(label = "DO Mouse") +
  ylab(label = "Proportion") +
  labs(fill = "Cell Clusters") +
  theme(axis.title=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 0.25),
        axis.ticks.length.x =  unit(0.25, "cm"),
        legend.position="bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=15),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_manual(name='Cell Clusters',
                     breaks=c('MPC', 'LMP', 'OBP', 'OB1', 'OB2', 'Ocy', 'MALP', 'Hem'),
                     values=c('MPC'='#4363d8', 
                              'LMP'='#f7c564', 
                              'OBP'='#9A6324',
                              'OB1'='#469990',
                              'OB2'='#f58231',
                              'Ocy'='#db4105',
                              'MALP'='#8762a6',
                              'Hem'='#c5c5c5')) +

  coord_flip()

############################################################
############################################################ 
############################################################

#make a subset of the object to include mouse 50 and mouse 233
sub_50 <- do80.combined_mes_subset[ , rownames(do80.combined_mes_subset@meta.data)[do80.combined_mes_subset@meta.data$mouse == 50]]
sub_233 <- do80.combined_mes_subset[ , rownames(do80.combined_mes_subset@meta.data)[do80.combined_mes_subset@meta.data$mouse == 233]]

sub_50[["legend"]] <- recode(sub_50$celltype_annotation, 
                                     "MPC" = paste0("MPC ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "MPC")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "MPC")/ncol(sub_50)*100, digits = 2)), "%", ")"),
                                     "LMP" = paste0("LMP ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "LMP")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "LMP")/ncol(sub_50)*100, digits = 2)), "%", ")"),
                                     "OBP" = paste0("OBP ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "OBP")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "OBP")/ncol(sub_50)*100, digits = 2)), "%", ")"),
                                     "OB1" = paste0("OB1 ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "OB1")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "OB1")/ncol(sub_50)*100, digits = 2)), "%", ")"),
                                     "OB2" = paste0("OB2 ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "OB2")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "OB2")/ncol(sub_50)*100, digits = 2)), "%", ")"),
                                     "Ocy" = paste0("Ocy ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "Ocy")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "Ocy")/ncol(sub_50)*100, digits = 2)), "%", ")"),
                                     "MALP" = paste0("MALP ", "(", paste(sum(sub_50@meta.data$celltype_annotation == "MALP")), ", ", paste(round(sum(sub_50@meta.data$celltype_annotation == "MALP")/ncol(sub_50)*100, digits = 2)), "%", ")"))

sub_233[["legend"]] <- recode(sub_233$celltype_annotation, 
                                      "MPC" = paste0("MPC ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "MPC")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "MPC")/ncol(sub_233)*100, digits = 2)), "%", ")"),
                                      "LMP" = paste0("LMP ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "LMP")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "LMP")/ncol(sub_233)*100, digits = 2)), "%", ")"),
                                      "OBP" = paste0("OBP ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "OBP")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "OBP")/ncol(sub_233)*100, digits = 2)), "%", ")"),
                                      "OB1" = paste0("OB1 ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "OB1")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "OB1")/ncol(sub_233)*100, digits = 2)), "%", ")"),
                                      "OB2" = paste0("OB2 ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "OB2")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "OB2")/ncol(sub_233)*100, digits = 2)), "%", ")"),
                                      "Ocy" = paste0("Ocy ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "Ocy")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "Ocy")/ncol(sub_233)*100, digits = 2)), "%", ")"),
                                      "MALP" = paste0("MALP ", "(", paste(sum(sub_233@meta.data$celltype_annotation == "MALP")), ", ", paste(round(sum(sub_233@meta.data$celltype_annotation == "MALP")/ncol(sub_233)*100, digits = 2)), "%", ")"))


cols_50 <- c('#4363d8','#9A6324','#f7c564', "#8762a6", "#db4105")
names(cols_50) <- c(paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend), "MPC")]),
                       paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend) , "OBP")]),
                       paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend) , "OB2")]),
                       paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend) , "LMP")]),
                       paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend) , "MALP")]),
                       paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend) , "OB1")]),
                       paste(levels(sub_50$legend)[startsWith(levels(sub_50$legend) , "Ocy")]))

cols_233 <- c('#4363d8','#9A6324',"#f58231",'#f7c564', "#8762a6", "#469990","#db4105")
names(cols_233) <- c(paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend), "MPC")]),
                        paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend) , "OBP")]),
                        paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend) , "OB2")]),
                        paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend) , "LMP")]),
                        paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend) , "MALP")]),
                        paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend) , "OB1")]),
                        paste(levels(sub_233$legend)[startsWith(levels(sub_233$legend) , "Ocy")]))

############################################################          ############################################################
############################################################  Fig 1D  ############################################################
############################################################          ############################################################

grid.arrange(
  DimPlot(sub_50, reduction = 'umap', group.by = "legend", label = F, raster=F, pt.size = 2, cols = cols_50) +
    ggtitle("DO Mouse 50") + 
    theme(title = element_text(size = 25),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 21),
          legend.title = element_text(size = 23, hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=5))), 
  DimPlot(sub_233, reduction = 'umap', group.by = "legend", label = F, raster=F, pt.size = 2, cols = cols_233) +
    ggtitle("DO Mouse 233") + 
    theme(title = element_text(size = 25),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 21),
          legend.title = element_text(size = 23, hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=5))),
  nrow = 1)

############################################################
############################################################
############################################################

CELLECT$`Adjusted P-value` <- stats::p.adjust(CELLECT$`P-value`, method = "bonferroni")

CELLECT_df <- data.frame()

for (i in 1:nrow(CELLECT)) {

  sub <- CELLECT[i,]
  
  if (sub$`Adjusted P-value` < 0.05) {
    sub$group <- "^"
  } else
    sub$group <- "-"
  
  CELLECT_df <- rbind(CELLECT_df, sub)
  
}

CELLECT_df$log <- -log10(CELLECT_df$`Adjusted P-value`)
x <- CELLECT_df[order(CELLECT_df$`Adjusted P-value`), ] 

colors <- numeric(2)
colors[x$group == "^"] <- "red"
colors[x$group == "-"] <- "black"


############################################################          ############################################################
############################################################  Fig 1E  ############################################################
############################################################          ############################################################

dotchart(x$log, 
         labels = x$`Cell cluster`,
         groups = as.factor(x$group),
         pt.cex = 1.5,
         pch = 19,
         color = colors,
         xlab = expression("-log"[10]*"(P"["adj"]*")"),
         ylab = "Cell Cluster") +
  abline(v = -log10(0.05), 
         col = "red",            
         lty = "dashed",         
         lwd = 1) 

############################################################
############################################################
############################################################

#subset to include only mesenchymal cell clusters without small outlier population of cells for trajectory analysis
umap_coord <- as.data.frame(do80.combined_mes_subset@reductions$umap@cell.embeddings)
umap_coord$cluster <- do80.combined_mes_subset@meta.data$seurat_clusters
nrow(umap_coord)
umap_coord <- subset(umap_coord, subset = UMAP_1 >-6 & UMAP_2 >-7)
nrow(umap_coord)

umap_coord_MPC <- subset(umap_coord, subset = cluster == 8)
nrow(umap_coord_MPC)
umap_coord_MPC <- subset(umap_coord_MPC, subset = UMAP_2 < -5)
nrow(umap_coord_MPC)

nrow(umap_coord)
umap_coord <- umap_coord[!rownames(umap_coord) %in% rownames(umap_coord_MPC),]
nrow(umap_coord)

do80.combined_mes_subset <- do80.combined_mes_subset[,colnames(do80.combined_mes_subset) %in% rownames(umap_coord)]

sce <- as.SingleCellExperiment(do80.combined_mes_subset)

sce <- slingshot(sce,
                 clusterLabels = sce$celltype_annotation,
                 reducedDim = reducedDim(sce, 'PCA')[,1:15],
                 allow.breaks = TRUE, 
                 start.clus = "MPC",
                 approx_points = 100)

umap_curve_embedding <- embedCurves(sce, newDimRed = reducedDims(sce)$UMAP)

############################################################          ############################################################
############################################################  Fig 3A  ############################################################
############################################################          ############################################################

colors <- c('#4363d8','#9A6324',"#f58231",'#f7c564', "#8762a6", "#469990","#db4105")
names(colors) <- c("MPC", "OBP", "OB2", "LMP", "MALP", "OB1", "Ocy")


plot(reducedDims(sce)$UMAP, 
     col = colors[as.character(sce$ident)],
     pch = 19, cex = 0.2, asp =0.75,
     ylab = "", yaxt='n',
     xlab="", xaxt='n')
lines(SlingshotDataSet(umap_curve_embedding), pch = 25, lwd=3.0,col = c("black"))

#legend
names <- c(paste0("MPC ", "(", paste(sum(sce$ident == "MPC")), ", ", paste(round(sum(sce$ident == "MPC")/ncol(sce)*100, digits = 2)), "%", ")"),
           paste0("OBP ", "(", paste(sum(sce$ident == "OBP")), ", ", paste(round(sum(sce$ident == "OBP")/ncol(sce)*100, digits = 2)), "%", ")"),
           paste0("OB2 ", "(", paste(sum(sce$ident == "OB2")), ", ", paste(round(sum(sce$ident == "OB2")/ncol(sce)*100, digits = 2)), "%", ")"),
           paste0("LMP ", "(", paste(sum(sce$ident == "LMP")), ", ", paste(round(sum(sce$ident == "LMP")/ncol(sce)*100, digits = 2)), "%", ")"),
           paste0("MALP ", "(", paste(sum(sce$ident == "MALP")), ", ", paste(round(sum(sce$ident == "MALP")/ncol(sce)*100, digits = 2)), "%", ")"),
           paste0("OB1 ", "(", paste(sum(sce$ident == "OB1")), ", ", paste(round(sum(sce$ident == "OB1")/ncol(sce)*100, digits = 2)), "%", ")"),
           paste0("Ocy ", "(", paste(sum(sce$ident == "Ocy")), ", ", paste(round(sum(sce$ident == "Ocy")/ncol(sce)*100, digits = 2)), "%", ")"))
clrs <- c('#4363d8','#9A6324',"#f58231",'#f7c564', "#8762a6", "#469990","#db4105")
ltype <- c(1, 1, 1, 1, 1, 1)
plot(NULL, xaxt='n', yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", title="Mesenchymal\nCell Clusters", legend = names, 
       pch=19, cex=1.25, pt.cex = 1.75, x.intersp = 1.25, y.intersp = 0.8,
       bty='n', col = clrs)

############################################################
############################################################
############################################################

#subset the data in preparation for the tradeSeq
sce_sub <- sce[,c(sample(x = colnames(sce), size = round(ncol(sce)/10, digits = 0), replace = F))]

counts <- as.matrix(assays(sce_sub)$counts)
pseudotime <- slingPseudotime(SlingshotDataSet(sce_sub), na = FALSE)
cellWeights <- slingCurveWeights(sce_sub)
cellWeightsFiltered <- cellWeights[colnames(counts), ]
pseudotimeFiltered <- pseudotime[colnames(counts), ]
countsFiltered <- subset.matrix(counts, select=rownames(cellWeightsFiltered))

evalK_10_sub <- evaluateK(counts = as.matrix(countsFiltered),
                          pseudotime = as.matrix(pseudotimeFiltered), 
                          cellWeights = as.matrix(cellWeightsFiltered),
                          k= 3:10, 
                          nGenes = 100,
                          verbose = TRUE)

gamFIT_10knots <- fitGAM(counts = as.matrix(countsFiltered), 
                         pseudotime = as.matrix(pseudotimeFiltered), 
                         cellWeights = as.matrix(cellWeightsFiltered), 
                         nknots = 10, 
                         sce = TRUE,
                         verbose = TRUE)

#establish the boundaries for the tradeseq analysis
#for lineage 1 (terminal cluster: Ocy)
MPC_to_LMP_lin1  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(0, 20), global = T, lineages = T, l2fc = 0.5)
LMP_to_OBP_lin1  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(20, 42), global = T, lineages = T, l2fc = 0.5)
OBP_to_OB1_lin1  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(42, 53), global = T, lineages = T, l2fc = 0.5)
OB1_to_Ocy_lin1  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(53, 62), global = T, lineages = T, l2fc = 0.5)
Ocy_to_end_lin1  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(62, 85), global = T, lineages = T, l2fc = 0.5)

#for lineage 2 (terminal cluster: OB2; branch point at OBP)
OBP_to_OB2_lin2  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(42, 56), global = T, lineages = T, l2fc = 0.5)
OB2_to_end_lin2  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(56, 85), global = T, lineages = T, l2fc = 0.5)

#for lineage 3 (terminal cluster: MALP; branch point between lin1 and lin3 at LMP)
LMP_to_MALP_lin3  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(20, 36), global = T, lineages = T, l2fc = 0.5)
MALP_to_end_lin3  <- startVsEndTest(gamFIT_10knots, pseudotimeValues = c(36, 134), global = T, lineages = T, l2fc = 0.5)

#for each tradeseq dfs, subset to include genes with an adjusted P-value < 0.05 
#lineage 1 as example

all_df_lin1 <- c(names(mget(ls(pattern = "_lin1"))))

for (a in all_df_lin1){
  data <- get(a)
  data$gene <- data$X
  data$pseudotime_boundary <- gsub(", ", "_", toString(strsplit(a, "_")[[1]][1:3]))
  data_polished <- subset(as.data.frame(data), select = c(pseudotime_boundary, gene, waldStat_lineage1, logFClineage1, pvalue_lineage1))
  data_polished$lin1_adj_pvalue <- p.adjust(data_polished$pvalue_lineage1, method = c("fdr"))
  data_polished <- subset(data_polished, subset = lin1_adj_pvalue <= 0.05)
  assign(paste(a, "_polished", sep = "_"), data_polished)
  rm(data, data_polished)
}

#combine all lin 1
all_df_lin1 <- c(names(mget(ls(pattern = "lin1_polished"))))

all_lin1_specific_tradeseq_polished_df <- data.frame()

for (a in all_df_lin1) {
  data <- get(a)
  all_lin1_specific_tradeseq_polished_df <- rbind(all_lin1_specific_tradeseq_polished_df, data)
}

#same to other lineages

############################################################          ############################################################
############################################################  Fig 3C  ############################################################
############################################################          ############################################################

Lars2_feature <- FeaturePlot(do80.combined_mes_subset, 
                             features = "Lars2", 
                             cols =c("lightgrey", "#4363d8"), 
                             pt.size = 0.6, 
                             label = F, raster = F, label.size = 2) + 
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")

Aspn_feature <- FeaturePlot(do80.combined_mes_subset, 
                            features = "Aspn", 
                            cols =c("lightgrey", "#FFDA03"), 
                            pt.size = 0.6, 
                            label = F, raster = F, label.size = 2) + 
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")

Ifitm5_feature <- FeaturePlot(do80.combined_mes_subset, 
                             features = "Ifitm5", 
                             cols =c("lightgrey", "#9A6324"), pt.size = 0.6, 
                             label = F, raster = F, label.size = 2) + 
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")


Bglap_feature <- FeaturePlot(do80.combined_mes_subset, 
                            features = "Bglap", 
                            cols =c("lightgrey", "#469990"),pt.size = 0.6, 
                            label = F, raster = F, label.size = 2) +  
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")

Phex_feature <- FeaturePlot(do80.combined_mes_subset, 
                            features = "Phex", 
                            cols =c("lightgrey", "#db4105"),pt.size = 0.6, 
                            label = F, raster = F, label.size = 2) +  
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")


Car12_feature <- FeaturePlot(do80.combined_mes_subset, 
                            features = "Car12", 
                            cols =c("lightgrey", "#f58231"),pt.size = 0.6, 
                            label = F, raster = F, label.size = 2) +  
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")

Adipoq_feature <- FeaturePlot(do80.combined_mes_subset, 
                             features = "Adipoq",
                             cols =c("lightgrey", "#8762a6"), pt.size = 0.6, 
                             label = F, raster = F, label.size = 2) +
  theme(plot.title = element_text(size = 40, face = "italic"), 
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(1, 'cm'),
        legend.text=element_text(size=20),
        legend.position = "right")

genes_plot_smooths <- c("Lars2", "Aspn", "Ifitm5", "Bglap", "Phex", "Car12", "Adipoq")

for (g in genes_plot_smooths) {

  models = gamFIT_10knots
  counts = countsFiltered
  gene = g
  nPoints = 100
  lwd = 2
  size = 2/3 
  xlab = "Pseudotime"
  ylab = "Log(expression + 1)" 
  border = FALSE
  alpha = 2/3 
  sample = 1 
  pointCol = NULL
  curvesCols = NULL
  plotLineages = TRUE 
  lineagesToPlot = c(1,2,3)
  
  # Predicting fits ----
  # lpmatrix given X and design
  predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL) {
    # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
    # INPUT:
    # lpmatrix is the linear predictor matrix of the GAM model
    # df is a data frame of values for which we want the lpmatrix
    # pseudotime is the n x l matrix of pseudotimes
    # conditions is the vector of conditions, if present.
    
    # if pseudotime is vector, make it a matrix.
    if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)
    
    condPresent <- !is.null(conditions)
    if(condPresent) nConditions <- nlevels(conditions)
    
    # for each curve, specify basis function IDs for lpmatrix
    allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")
    
    if(!condPresent){
      lineages <- sub(pattern = "s\\(", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
      lineages <- sub(pattern = "\\):.*", replacement = "",
                      x = lineages)
      nCurves <- length(unique(lineages))
      for (ii in seq_len(nCurves)) {
        assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
      }
    } else if(condPresent){
      lineages <- sub(pattern = "s\\(t", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
      lineages <- sub(pattern = "\\):.*", replacement = "",
                      x = lineages)
      nLineages <- length(unique(lineages))
      curves <- sub(pattern = ".*:l", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
      curves <- sub(pattern = "\\..*", replacement = "",
                    x = curves)
      nCurves <- length(unique(curves))
      for (ii in seq_len(nLineages)) {
        for(kk in seq_len(nConditions))
          assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
      }
    }
    
    
    # specify lineage assignment for each cell (i.e., row of lpmatrix)
    if(!condPresent){
      lineageID <- apply(lpmatrix, 1, function(x){
        for (ii in seq_len(nCurves)) {
          if (!all(x[get(paste0("id", ii))] == 0)) {
            return(ii)
          }
        }
      })
    } else if(condPresent){
      # first number is lineage, second number is condition.
      lineageID <- apply(lpmatrix, 1, function(x){
        for (ii in seq_len(nLineages)) {
          # loop over lineages
          for(kk in seq_len(nConditions)){
            # loop over conditions
            if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
              return(as.numeric(paste0(ii, kk)))
            }
          }
        }
      })
    }
    
    
    # fit splinefun for each basis function based on assigned cells
    if(!condPresent) {
      for (ii in seq_len(nCurves)) { # loop over curves
        for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
          assign(paste0("l",ii,".",jj),
                 stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                  y = lpmatrix[lineageID == ii, #only cells for lineage
                                               get(paste0("id", ii))[jj]],
                                  ties = mean)) #basis function
        }
      }
    } else if(condPresent) {
      for (ii in  seq_len(nLineages)) {
        # loop over curves
        for(kk in seq_len(nConditions)){
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
            #within curve, loop over basis functions
            assign(paste0("l",ii, "_", kk,".",jj),
                   stats::splinefun(
                     x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                     y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                  get(paste0("id", ii, "_", kk))[jj]],
                     ties = mean)) #basis function
          }
        }
      }
    }
    
    
    # use input to estimate X for each basis function
    Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
    if(!condPresent){
      for (ii in seq_len(nCurves)) { # loop over curves
        if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
            f <- get(paste0("l", ii, ".", jj))
            Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    } else if(condPresent){
      # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
      for (ii in seq_len(nLineages)) {
        # loop over curves
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
            for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
              # within curve, loop over basis functions
              f <- get(paste0("l", ii, "_", kk, ".", jj))
              Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
            }
          }
        }
      }
    }
    
    
    # add fixed covariates as in df
    dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
    dfOffsetID <- grep(x = colnames(df), pattern = "offset")
    Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]
    
    # return
    colnames(Xout) <- colnames(lpmatrix)
    return(Xout)
  }
  
  # get the first non-errored fit in models
  .getModelReference <- function(models){
    for (i in seq_len(length(models))) {
      m <- models[[i]]
      if (is(m)[1] != "try-error") return(m)
    }
    stop("All models errored")
  }
  
  .getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
    vars <- dm[1, ]
    if ("y" %in% colnames(vars)) {
      vars <- vars[!colnames(vars) %in% "y"]
      off <- 1
    } else {
      off <- 0
    }
    offsetId <- grep(x = colnames(vars), pattern = "offset")
    offsetName <- colnames(vars)[offsetId]
    offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
    names(vars)[offsetId] <- offsetName
    # set all times on 0
    vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
    # set all lineages on 0
    vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
    # duplicate to nPoints
    vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
    # set range of pseudotime for lineage of interest
    if (is.null(conditionId)) {
      lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
    } else {
      lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                          "_", conditionId, "$"))
    }
    if (length(lineageIds) == 1){
      lineageData <- dm[dm[, lineageIds + off] == 1,
                        paste0("t", lineageId)]
    } else {
      lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                        paste0("t", lineageId)]
    }
    # make sure lineage starts at zero
    if(min(lineageData) / max(lineageData) < .01) {
      lineageData[which.min(lineageData)] <- 0
    }
    vars[, lineageIds] <- 1 / length(lineageIds)
    # set lineage
    vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                          max(lineageData),
                                          length = nPoints)
    # set offset
    vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                         pattern = "offset")])
    return(vars)
  }
  
  
  #input is singleCellExperiment object.
  if (is.null(names(models))) {
    rownames(models) <- rownames(counts) <- seq_len(nrow(models))
    message(paste0(
      "The sce object has no rownames. Assuming that the counts and the sce ",
      "objects are ordered in the same way"))
  }
  if (length(gene) > 1) stop("Only provide a single gene's ID with the ",
                             "gene argument.")
  # check if all gene IDs provided are present in the models object.
  if (is(gene, "character")) {
    if (!all(gene %in% names(models))) {
      stop("The gene ID is not present in the models object.")
    }
    id <- which(names(models) %in% gene)
  } else id <- gene
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  y <- unname(counts[names(models),][id,])
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  # pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
  #                                     pattern = "pseudotime")]
  pseudotime <- pseudotimeFiltered
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id,]
  
  
  #construct time variable based on cell assignments.
  lcol <- timeAll <- rep(0, nrow(dm))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(dm))) {
      if (dm[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- dm[ii, paste0("t", jj)]
        lcol[ii] <- jj
      } else {
        next
      }
    }
  }
  
  if (!is.null(pointCol)) {
    if (length(pointCol) == 1) {
      col <- colData(models)[,pointCol]
    } else if (length(pointCol) == ncol(models)) {
      col <- pointCol
    } else {
      col <- lcol
      message(paste("pointCol should have length of either 1 or the number of cells,",
                    "reverting to default color scheme."))
    }
  } else {
    col <- lcol
  }
  
  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "pCol" = as.character(col),
                   "lineage" = as.character(lcol))
  rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
  df <- df[rows, ]
  if(!is.null(lineagesToPlot)){
    df <- df[df$lineage %in% lineagesToPlot,]
  }
  p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if(is.null(pointCol)){
    p <- p +
      geom_point(size = size, aes(col = lineage)) +
      #add this to the geom_point parameters if necessary: colour="#ffffff"
      scale_color_viridis_d(alpha = alpha)
  } else {
    p <- p +
      geom_point(size = size, alpha = alpha, aes(col = pCol)) +
      scale_color_discrete() +
      labs(col = "Cell labels")
  }
  
  
  
  # predict and plot smoothers across the range
  if (plotLineages) {
    if (!is.null(curvesCols)) {
      if (length(curvesCols) != nCurves) {
        curvesCols <- viridis::viridis(nCurves)
        message("Incorrect number of lineage colors. Default to viridis")
      }
    } else {
      curvesCols <- viridis::viridis(nCurves)
    }
    if(is.null(lineagesToPlot)){
      lineagesToPlot <- seq_len(nCurves)
    }
    for (jj in lineagesToPlot) {
      df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime)
      yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
      if (border) {
        p <- p +
          geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                      "gene_count" = yhat,
                                      "lineage" = as.character(jj),
                                      "pCol" = as.character(jj)),
                    lwd = lwd + 1, colour = "white")
      }
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj),
                                    "pCol" = as.character(jj)),
                  lwd = lwd, col = curvesCols[jj])
    }
  }
  
  assign(paste(gene, "smooth_plot", sep = "_"), p)
}


genes_plot_smooths_graph  <- c("Lars2_smooth_plot", "Aspn_smooth_plot", "Ifitm5_smooth_plot", 
                               "Bglap_smooth_plot", "Phex_smooth_plot", "Car12_smooth_plot", "Adipoq_smooth_plot")

for (i in genes_plot_smooths_graph) {

  final_plot <- get(i) +
  #ggtitle(paste(gene)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2", "Adipogenic"), values= c("#552586","#009999","#FFFF33")) +
  #scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2"), values= c("#552586","#009999")) +
  #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 1"), values= c("#552586")) +
  #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 2"), values= c("#009999")) +
  #scale_color_manual(name = "Cell Lineages", labels = c("Adipogenic"), values= c("#FFFF33")) + 
  scale_x_continuous(limits = c(0, 85)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = "italic"),
        axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.position = "none",
        legend.title = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",color = NA)) + 
  
  geom_vline(xintercept = 0, color = 'red')  +
  annotate("label", x=0, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MPC", size = 5) +
  
  geom_vline(xintercept = 20, color = 'red')  +
  annotate("label", x=20, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "LMP", size = 5) +
  
  geom_vline(xintercept = 36, color = 'red')  +
  annotate("label", x=32, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MALP", size = 5) +
  
  geom_vline(xintercept=42, color = 'red') + 
  annotate("label", x=42, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "OBP", size = 5) +
  
  geom_vline(xintercept=55, color = 'red') +
  annotate("label", x=54, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.35), label= "OB1/\nOB2", size = 5) +
  
  geom_vline(xintercept=62, color = 'red') +
  annotate("label", x=62, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "Ocy", size = 5)

  #xlim(NA, NA)
  assign(paste(strsplit(i, split = "_")[[1]][1], "final", sep = "_"), final_plot)

}

grid.arrange(Lars2_feature, Aspn_feature, Ifitm5_feature, Bglap_feature, Phex_feature, Car12_feature, Adipoq_feature, 
             Lars2_final, Aspn_final, Ifitm5_final, Bglap_final, Phex_final, Car12_final, Adipoq_final, nrow = 2)

############################################################          ############################################################
############################################################  Fig 4A  ############################################################
############################################################          ############################################################

genes_plot_smooths <- c("Pkm", "S100a1")

for (g in genes_plot_smooths) {
  
    models = gamFIT_10knots
    counts = countsFiltered
    gene = g
    nPoints = 100
    lwd = 2
    size = 2/3 
    xlab = "Pseudotime"
    ylab = "Log(expression + 1)" 
    border = FALSE
    alpha = 2/3 
    sample = 1 
    pointCol = NULL
    curvesCols = NULL
    plotLineages = TRUE 
    lineagesToPlot = c(1,2,3)
    
    # Predicting fits ----
    # lpmatrix given X and design
    predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL) {
      # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
      # INPUT:
      # lpmatrix is the linear predictor matrix of the GAM model
      # df is a data frame of values for which we want the lpmatrix
      # pseudotime is the n x l matrix of pseudotimes
      # conditions is the vector of conditions, if present.
      
      # if pseudotime is vector, make it a matrix.
      if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)
      
      condPresent <- !is.null(conditions)
      if(condPresent) nConditions <- nlevels(conditions)
      
      # for each curve, specify basis function IDs for lpmatrix
      allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")
      
      if(!condPresent){
        lineages <- sub(pattern = "s\\(", replacement = "",
                        x = colnames(lpmatrix[,allBs]))
        lineages <- sub(pattern = "\\):.*", replacement = "",
                        x = lineages)
        nCurves <- length(unique(lineages))
        for (ii in seq_len(nCurves)) {
          assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
        }
      } else if(condPresent){
        lineages <- sub(pattern = "s\\(t", replacement = "",
                        x = colnames(lpmatrix[,allBs]))
        lineages <- sub(pattern = "\\):.*", replacement = "",
                        x = lineages)
        nLineages <- length(unique(lineages))
        curves <- sub(pattern = ".*:l", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
        curves <- sub(pattern = "\\..*", replacement = "",
                      x = curves)
        nCurves <- length(unique(curves))
        for (ii in seq_len(nLineages)) {
          for(kk in seq_len(nConditions))
            assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
        }
      }
      
      
      # specify lineage assignment for each cell (i.e., row of lpmatrix)
      if(!condPresent){
        lineageID <- apply(lpmatrix, 1, function(x){
          for (ii in seq_len(nCurves)) {
            if (!all(x[get(paste0("id", ii))] == 0)) {
              return(ii)
            }
          }
        })
      } else if(condPresent){
        # first number is lineage, second number is condition.
        lineageID <- apply(lpmatrix, 1, function(x){
          for (ii in seq_len(nLineages)) {
            # loop over lineages
            for(kk in seq_len(nConditions)){
              # loop over conditions
              if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
                return(as.numeric(paste0(ii, kk)))
              }
            }
          }
        })
      }
      
      
      # fit splinefun for each basis function based on assigned cells
      if(!condPresent) {
        for (ii in seq_len(nCurves)) { # loop over curves
          for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
            assign(paste0("l",ii,".",jj),
                   stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                    y = lpmatrix[lineageID == ii, #only cells for lineage
                                                 get(paste0("id", ii))[jj]],
                                    ties = mean)) #basis function
          }
        }
      } else if(condPresent) {
        for (ii in  seq_len(nLineages)) {
          # loop over curves
          for(kk in seq_len(nConditions)){
            for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
              #within curve, loop over basis functions
              assign(paste0("l",ii, "_", kk,".",jj),
                     stats::splinefun(
                       x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                       y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                    get(paste0("id", ii, "_", kk))[jj]],
                       ties = mean)) #basis function
            }
          }
        }
      }
      
      
      # use input to estimate X for each basis function
      Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
      if(!condPresent){
        for (ii in seq_len(nCurves)) { # loop over curves
          if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
            for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
              f <- get(paste0("l", ii, ".", jj))
              Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
            }
          }
        }
      } else if(condPresent){
        # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
        for (ii in seq_len(nLineages)) {
          # loop over curves
          for(kk in seq_len(nConditions)){
            # loop over conditions
            if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
              for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
                # within curve, loop over basis functions
                f <- get(paste0("l", ii, "_", kk, ".", jj))
                Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
              }
            }
          }
        }
      }
      
      
      # add fixed covariates as in df
      dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
      dfOffsetID <- grep(x = colnames(df), pattern = "offset")
      Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]
      
      # return
      colnames(Xout) <- colnames(lpmatrix)
      return(Xout)
    }
    
    # get the first non-errored fit in models
    .getModelReference <- function(models){
      for (i in seq_len(length(models))) {
        m <- models[[i]]
        if (is(m)[1] != "try-error") return(m)
      }
      stop("All models errored")
    }
    
    .getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
      vars <- dm[1, ]
      if ("y" %in% colnames(vars)) {
        vars <- vars[!colnames(vars) %in% "y"]
        off <- 1
      } else {
        off <- 0
      }
      offsetId <- grep(x = colnames(vars), pattern = "offset")
      offsetName <- colnames(vars)[offsetId]
      offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
      names(vars)[offsetId] <- offsetName
      # set all times on 0
      vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
      # set all lineages on 0
      vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
      # duplicate to nPoints
      vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
      # set range of pseudotime for lineage of interest
      if (is.null(conditionId)) {
        lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
      } else {
        lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                            "_", conditionId, "$"))
      }
      if (length(lineageIds) == 1){
        lineageData <- dm[dm[, lineageIds + off] == 1,
                          paste0("t", lineageId)]
      } else {
        lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                          paste0("t", lineageId)]
      }
      # make sure lineage starts at zero
      if(min(lineageData) / max(lineageData) < .01) {
        lineageData[which.min(lineageData)] <- 0
      }
      vars[, lineageIds] <- 1 / length(lineageIds)
      # set lineage
      vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                            max(lineageData),
                                            length = nPoints)
      # set offset
      vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                           pattern = "offset")])
      return(vars)
    }
    
    
    #input is singleCellExperiment object.
    if (is.null(names(models))) {
      rownames(models) <- rownames(counts) <- seq_len(nrow(models))
      message(paste0(
        "The sce object has no rownames. Assuming that the counts and the sce ",
        "objects are ordered in the same way"))
    }
    if (length(gene) > 1) stop("Only provide a single gene's ID with the ",
                               "gene argument.")
    # check if all gene IDs provided are present in the models object.
    if (is(gene, "character")) {
      if (!all(gene %in% names(models))) {
        stop("The gene ID is not present in the models object.")
      }
      id <- which(names(models) %in% gene)
    } else id <- gene
    
    dm <- colData(models)$tradeSeq$dm # design matrix
    y <- unname(counts[names(models),][id,])
    X <- colData(models)$tradeSeq$X # linear predictor
    slingshotColData <- colData(models)$crv
    # pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
    #                                     pattern = "pseudotime")]
    pseudotime <- pseudotimeFiltered
    if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
    betaMat <- rowData(models)$tradeSeq$beta[[1]]
    beta <- betaMat[id,]
    
    
    #construct time variable based on cell assignments.
    lcol <- timeAll <- rep(0, nrow(dm))
    for (jj in seq_len(nCurves)) {
      for (ii in seq_len(nrow(dm))) {
        if (dm[ii, paste0("l", jj)] == 1) {
          timeAll[ii] <- dm[ii, paste0("t", jj)]
          lcol[ii] <- jj
        } else {
          next
        }
      }
    }
    
    if (!is.null(pointCol)) {
      if (length(pointCol) == 1) {
        col <- colData(models)[,pointCol]
      } else if (length(pointCol) == ncol(models)) {
        col <- pointCol
      } else {
        col <- lcol
        message(paste("pointCol should have length of either 1 or the number of cells,",
                      "reverting to default color scheme."))
      }
    } else {
      col <- lcol
    }
    
    # plot raw data
    df <- data.frame("time" = timeAll,
                     "gene_count" = y,
                     "pCol" = as.character(col),
                     "lineage" = as.character(lcol))
    rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
    df <- df[rows, ]
    if(!is.null(lineagesToPlot)){
      df <- df[df$lineage %in% lineagesToPlot,]
    }
    p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
      labs(x = xlab, y = ylab) +
      theme_classic()
    if(is.null(pointCol)){
      p <- p +
        geom_point(size = size, aes(col = lineage)) +
        #add this to the geom_point parameters if necessary: colour="#ffffff"
        scale_color_viridis_d(alpha = alpha)
    } else {
      p <- p +
        geom_point(size = size, alpha = alpha, aes(col = pCol)) +
        scale_color_discrete() +
        labs(col = "Cell labels")
    }
    
    
    
    # predict and plot smoothers across the range
    if (plotLineages) {
      if (!is.null(curvesCols)) {
        if (length(curvesCols) != nCurves) {
          curvesCols <- viridis::viridis(nCurves)
          message("Incorrect number of lineage colors. Default to viridis")
        }
      } else {
        curvesCols <- viridis::viridis(nCurves)
      }
      if(is.null(lineagesToPlot)){
        lineagesToPlot <- seq_len(nCurves)
      }
      for (jj in lineagesToPlot) {
        df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
        Xdf <- predictGAM(lpmatrix = X,
                          df = df,
                          pseudotime = pseudotime)
        yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
        if (border) {
          p <- p +
            geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                        "gene_count" = yhat,
                                        "lineage" = as.character(jj),
                                        "pCol" = as.character(jj)),
                      lwd = lwd + 1, colour = "white")
        }
        p <- p +
          geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                      "gene_count" = yhat,
                                      "lineage" = as.character(jj),
                                      "pCol" = as.character(jj)),
                    lwd = lwd, col = curvesCols[jj])
      }
    }
    
    assign(paste(gene, "smooth_plot", sep = "_"), p)
    
    final_plot <- p +
      ggtitle(paste(gene)) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2", "Adipogenic"), values= c("#552586","#009999","#FFFF33")) +
      #scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2"), values= c("#552586","#009999")) +
      #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 1"), values= c("#552586")) +
      #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 2"), values= c("#009999")) +
      #scale_color_manual(name = "Cell Lineages", labels = c("Adipogenic"), values= c("#FFFF33")) + 
      scale_x_continuous(limits = c(0, 85)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 40, face = "italic"),
            axis.title.x = element_text(size = 25),
            axis.text.x = element_text(size = 25),
            axis.title.y = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.ticks = element_line(linewidth = 1),
            axis.ticks.length  = unit(0.25, "cm"),
            legend.position = "none",
            legend.title = element_text(size =15),
            legend.text = element_text(size = 15),
            legend.background = element_rect(fill = "transparent"),
            panel.background = element_rect(fill = "transparent"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent",color = NA)) + 
      
      geom_vline(xintercept = 0, color = 'red')  +
      annotate("label", x=0, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MPC", size = 7) +
      
      geom_vline(xintercept = 20, color = 'red')  +
      annotate("label", x=20, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "LMP", size = 7) +
      
      geom_vline(xintercept = 36, color = 'red')  +
      annotate("label", x=33, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MALP", size = 7) +
      
      geom_vline(xintercept=42, color = 'red') + 
      annotate("label", x=42, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "OBP", size = 7) +
      
      geom_vline(xintercept=55, color = 'red') +
      annotate("label", x=54, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.35), label= "OB1/\nOB2", size = 7) +
      
      geom_vline(xintercept=62, color = 'red') +
      annotate("label", x=62, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "Ocy", size = 7)
    
    #xlim(NA, NA)
    final_plot
    assign(paste(gene, "final_plot", sep = "_"), final_plot)
}

grid.arrange(Pkm_final_plot, S100a1_final_plot, nrow = 1)

############################################################
############################################################ 
############################################################

#for each mouse, generate a df which contains the average RNA counts for all cells derived from each mouse for each cluster.

for (i in unique(do80.combined_mes_subset@meta.data$mouse)) {
  mouse <- subset(do80.combined_mes_subset, subset = mouse == i)
  av <- AverageExpression(mouse, slot = "counts")
  av <- as.data.frame(av)
  assign(paste("av", i, sep = "_"), av)
  rm(av)
}

av_list <- mget(ls(pattern="^av"))

for (i in c(names(av_list[[1]]))) {
  
  for (x in names(av_list)) {
    sub <- av_list[[x]]
    
    if (i %in% colnames(sub) == T) {
      sub <- subset(sub, select = i)
      assign(paste(i, x, "sub", sep = "_"), sub)
    } else
      sub
  }
}

av_list_cluster_OB1 <- mget(ls(pattern="^RNA.OB1"))
av_list_cluster_OB2 <- mget(ls(pattern="^RNA.OB2"))
av_list_cluster_LMP <- mget(ls(pattern="^RNA.LMP"))
av_list_cluster_Ocy <- mget(ls(pattern="^RNA.Ocy"))
av_list_cluster_MALP <- mget(ls(pattern="^RNA.MALP"))
av_list_cluster_MPC <- mget(ls(pattern="^RNA.MPC"))
av_list_cluster_OBP <- mget(ls(pattern="^RNA.OBP"))

av_list_clusters <- mget(ls(pattern="^av_list_cluster_"))

for (a in 1:length(names(av_list_clusters))) {
  bound <- do.call(cbind, av_list_clusters[[a]])
  colnames(bound) <- names(av_list_clusters[[a]])
  bound <- data.table::setnames(bound, sub('RNA.', '', names(bound), fixed = TRUE))
  bound <- data.table::setnames(bound, sub('_sub', '', names(bound), fixed = TRUE))
  assign(paste("bound", print((stringr::str_split(colnames(bound), "_", simplify=T)[,1])[1]), sep = "_"), bound)
  rm(bound)
}

bounds <- c(names(mget(ls(pattern = "bound_"))))

#clean up column names. 
for (b in bounds) {
  data <- get(b)
  colnames(data) <- sub(".*_", "", colnames(data))
  assign(paste(b), data)
}

#make a table of cell count (mouse x cluster)
cellcount_table <- table(do80.combined_mes_subset@meta.data$mouse, do80.combined_mes_subset@active.ident)
cellcount_table_df <- as.data.frame.matrix(cellcount_table)
cellcount_table_df$mouse <- rownames(cellcount_table_df)

#get df of metadata, but only include columns for mouse, pool, sex
meta <- subset(do80.combined_mes_subset@meta.data, select = c(mouse, pool, sex))
meta <- meta[!duplicated(meta), ]
#for the NA in place of the pool information for the 5DO mice, add 0
meta$pool[is.na(meta$pool)] <- 0

#merge dfs
cellcount_table_df <- merge(meta, cellcount_table_df, by="mouse")

MPC_meta <- subset(cellcount_table_df, subset = MPC >= 5)
LMP_meta <- subset(cellcount_table_df, subset = LMP >= 5)
OBP_meta <- subset(cellcount_table_df, subset = OBP >= 5)
OB1_meta <- subset(cellcount_table_df, subset = OB1 >= 5)
OB2_meta <- subset(cellcount_table_df, subset = OB2 >= 5)
Ocy_meta <- subset(cellcount_table_df, subset = Ocy >= 5)
MALP_meta <- subset(cellcount_table_df, subset = MALP >= 5)

#check
bound_MPC <- subset(bound_MPC, select = unique(MPC_meta$mouse))
bound_LMP <- subset(bound_LMP, select = unique(LMP_meta$mouse))
bound_MALP <- subset(bound_MALP, select = unique(MALP_meta$mouse))
bound_OBP <- subset(bound_OBP, select = unique(OBP_meta$mouse))
bound_OB1 <- subset(bound_OB1, select = unique(OB1_meta$mouse))
bound_OB2 <- subset(bound_OB2, select = unique(OB2_meta$mouse))
bound_Ocy <- subset(bound_Ocy, select = unique(Ocy_meta$mouse))

#round matrices
for (b in bounds){
  data <- get(b)
  data <- round(data, digits = 2)
  assign(paste(b, "rounded", sep = "_"), data)
  rm(data, b)
}

#filter matrices to only include those genes that have at least 15 non-zero values across all mice. 
rounds <- c(names(mget(ls(pattern = "_rounded"))))

for (r in rounds) {
  matrix <- get(r)
  #add a column that counts the non-zero values
  matrix$count_nonzeros <- rowSums(matrix != 0)
  #subset the matrix to only include those genes with a non-zero count greater than or equal to 15
  matrix <- subset(matrix, count_nonzeros >= 15)
  #remove the count_nonzero column
  matrix <- subset(matrix, select = -c(count_nonzeros))
  assign(paste(r, "filtered", sep = "_"), matrix)
  rm(matrix, r)
}

#generate CPMs 
filtered <- c(names(mget(ls(pattern = "_filtered"))))

for (f in filtered) {
  cpm_count <- get(f)
  
  for (c in colnames(cpm_count)) {
    cpm_count[, c] <- (cpm_count[, c] / (sum(cpm_count[, c])/1e6))
  }
  
  assign(paste(f, "CPM", sep = "_"), cpm_count)
  rm(cpm_count)
}

cpms <- c(names(mget(ls(pattern = "_CPM"))))

for (c in cpms) {
  data <- get(c)
  data <- round(data, digits = 0)
  assign(paste(c, "integer", sep = "_"), data)
  rm(data, c)
}

#perform vst on the matrices
#columns = samples, rows = genes
integers <- c(names(mget(ls(pattern = "_integer"))))
for (i in integers) {
  vst <- DESeq2::varianceStabilizingTransformation(as.matrix(get(i)))
  assign(paste(i, "vst", sep = "_"), vst)
  rm(vst, i)
}

#perform quantile normalization on matrices
#rows = genes, columns = samples
vsts <- c(names(mget(ls(pattern = "_vst"))))

for (v in vsts) {
  data <- get(v)
  data_QN <- normalize.quantiles(as.matrix(data))
  rownames(data_QN) <- rownames(data)
  colnames(data_QN) <- colnames(data)
  assign(paste(v, "QN", sep = "_"), data_QN)
  rm(data_QN)
}

############################################################
############################################################ 
############################################################

#average expression matrices for each cluster were used as input to the eQTL analysis using the R package, qtl2
#example below for MPC cell cluster, but performed for each average expression matrix for each of the mesenchmyal cell clusters

celltype <- "MPC"

#transpose the average expression count matrix generated for the cluster such that the DO mouse are the rows and the genes are the column names
matrix <- t(bound_MPC_rounded_filtered_CPM_integer_vst_QN)

#use cross file generated previously and subset to only include those mice contributing to the average expression matrix for the cluster
subset_cross <- subset(cross_basic, ind = rownames(matrix))

#use allele probabilities generated previously and subset to only include those mice contributing to the average expression matrix for the cluster
subset_apr <- apr[rownames(matrix),]

#use kinship matrix generated previously and subset to only include those contributing to the average expression matrix for the cluster
subset_k_loco <- list()
for (k in names(k_loco)) {
  subset_k_loco[[k]] <- k_loco[[k]][rownames(matrix), rownames(matrix)]
}

#add the matrix to the cross files
subset_cross$scRNA_seq_matrix <- matrix

#get the X chrom covars from the cross file
#1 = male, 0 = female
Xcovar <- get_x_covar(subset_cross)

#create a covar object from covariates in cross file
#must be numeric
covar = as.matrix(subset_cross$covar)

#convert sex to 1's and 0's
#1 = male, 0 = female
covar[,"sex"] = (covar[,"sex"] == "M")*1

#sac date to factors
covar[,1] = as.factor(covar[,1])

#generation to factors
covar[,6] = as.factor(covar[,6])
covar = apply(covar,2,as.numeric)

#make the rownames of the covariate file have the same rownames as the cross file
rownames(covar) <- rownames(subset_cross$covar)

out_eqtl_local <- scan1(genoprobs = subset_apr,
                        pheno = subset_cross$scRNA_seq_matrix,
                        kinship = subset_k_loco,
                        Xcovar = Xcovar,
                        addcovar = covar,
                        cores = 15)

local_eqtl_peaks <- find_peaks(scan1_output = out_eqtl_local,
                               map = subset_cross$pmap,
                               threshold=4,
                               drop=1.5)

assign(paste(celltype, "out_eqtl_local", sep = "_"), out_eqtl_local)
assign(paste(celltype, "local_eqtl_peaks", sep = "_"), local_eqtl_peaks)

scan1_object <- get(paste(celltype, "out_eqtl_local", sep = "_"))
local_eqtl_peaks <- get(paste(celltype, "local_eqtl_peaks", sep = "_"))

#load the annotation file which contains coordinate information about genes
chr = c(seq(1:19),"X")
annot_file = annot_file[which(annot_file$Reference %in% chr),]

#define function for calculating eqtl

getLocalEqtl = function(peaks,geneAnnot,lodThreshAuto, lodThreshX, localDistFromStart,geneCol1,geneCol2){
  #merge the annotation file and the output of the find_peaks function using the gene name to filter genes only present in both files
  out = merge(peaks,geneAnnot,by.x = geneCol1, by.y = geneCol2)
  #multiply the peak positions by a million to get in base pairs (bp)
  #then calculate distance (in bp) of the starting position of each gene to the position of the peak
  out$dist_start = abs((as.numeric(out$pos*1000000)) - as.numeric(out$Start))
  #keep only the peaks for the genes that correspond to the chromosome on which the gene is actually located
  out = out[which(out$chr == out$Reference),]
  #keep only the genes that have peaks within x (bp) from the TSS of the gene
  out = out[which(out$dist_start <= localDistFromStart),] 
  #get a vector of indexes (corresponding to the rows of the df) which keeps only the genes located on autosomal chromosomes that also meet the specified LOD threshold
  idx = which((out$chr != "X") & (out$lod >= lodThreshAuto))
  #append to that vector the genes located on the X chromosome that also meet the specified LOD threshold
  idx = c(idx, which((out$chr == "X") & (out$lod >= lodThreshX)))
  #filter from the df only those genes in the vector that list the genes meeting the specified LOD thresholds for both the autosomal and X chromosomes
  out = out[idx,]
  return(out)
}

#load("./MPC_local_eqtl_peaks.Rdata")

cluster_peaks <- c(names(mget(ls(pattern = "_eqtl_peaks"))))

for (c in cluster_peaks) {
  peaks <- get(c)
  local_eqtl <- getLocalEqtl(peaks,
                             annot_file,
                             lodThreshAuto = 9.68,
                             lodThreshX = 9.49,
                             localDistFromStart = 1000000,
                             geneCol1 = "lodcolumn",
                             geneCol2 = "Gene.Name")
  assign(paste(c, "pruned", sep = "_"), local_eqtl)
}

pruned_peaks <- c(names(mget(ls(pattern = "_pruned"))))

#add a column that says the celltype for each peak
pruned_peaks_df = data.frame()
for (p in pruned_peaks){
  data <- get(p)
  data$celltype <- as.character((strsplit(p, "_")[[1]][1]))
  pruned_peaks_df <- rbind(pruned_peaks_df, data)
}

############################################################              ############################################################
############################################################  Fig 4B Pkm  ############################################################
############################################################              ############################################################

cell_type_to_examine <- "LMP"
gene <- "Pkm"
chr <- 9
pos <- 59.74122

coef <- scan1blup(genoprobs = LMP_subset_apr[,chr], 
                  pheno =  LMP_subset_cross$scRNA_seq_matrix[,gene], 
                  kinship = LMP_subset_k_loco[[chr]],
                  addcovar = LMP_covar)

plot_coefCC(x = coef, 
            map = LMP_subset_cross$pmap, 
            lodcolumn = gene, 
            scan1_output = LMP_scan1_object,
            main = paste(cell_type_to_examine, "eQTL for", gene, sep = " "),
            legend = NULL)

locus_genotypes <- maxmarg(geno_probs, 
                           LMP_subset_cross$pmap, 
                           chr = chr, 
                           pos= pos, 
                           return_char=TRUE, 
                           tol = 0.5)

############################################################              ############################################################
############################################################  Fig 4C Pkm  ############################################################
############################################################              ############################################################

plot_pxg(locus_genotypes, LMP_subset_cross$scRNA_seq_matrix[,gene], ylab = paste(gene, "Expression"))

############################################################              ############################################################
############################################################  Fig 4D Pkm  ############################################################
############################################################              ############################################################

props <- getTransformedProps(clusters = do80.combined_all_clusters$celltype_annotation_props, 
                             sample = do80.combined_all_clusters$mouse, 
                             transform="asin")
mouse_props <- as.data.frame.matrix(t(props$TransformedProps))

DO_table <- data.frame(letter = c("A", "B", "C", "D", "E", "F", "G", "H"),
                       founder = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"))

interest <- "G"

locus_genotypes <- as.data.frame(locus_genotypes)
locus_genotypes <- na.omit(locus_genotypes)
nrow(locus_genotypes)

propeller_geno <- merge(locus_genotypes, mouse_props, by = "row.names", all.x = T)
propeller_geno$mouse <- propeller_geno$Row.names
propeller_geno <- subset(propeller_geno, select = -c(Row.names))
head(propeller_geno)

#add M/F and pool data
meta <- as.data.frame(do80.combined_all_clusters@meta.data)
deduped.data <- unique(meta[ , c("mouse", "pool", "sex")] )
#check
nrow(deduped.data)
#make DO80 pool 0
deduped.data[is.na(deduped.data)] <- "other"

#get props and geno in the same df
propeller_geno <- merge(propeller_geno, deduped.data, by = "mouse", all.x = T)
nrow(propeller_geno)
head(propeller_geno)

#add columns to df
propeller_geno$pool <- paste0("pool_", propeller_geno$pool)
propeller_geno$pool <- as.factor(propeller_geno$pool)
propeller_geno$sex <- as.factor(propeller_geno$sex)

#add new column 
propeller_geno$geno_group <- NA


for (r in rownames(propeller_geno)) {
  if (str_detect(propeller_geno[r,"locus_genotypes"], paste(strsplit(interest[1], "")[[1]][1])) == TRUE) {
    propeller_geno[r,"geno_group"] <- paste( (subset(DO_table, subset = letter == paste(strsplit(interest[1], "")[[1]][1])))[1,2],
                                             "*",
                                             sep = "/" )
  } else 
    propeller_geno[r,"geno_group"] <- "other"
  
}

head(propeller_geno)
table(propeller_geno$geno_group)
alleles <- levels(as.data.frame(table(propeller_geno$geno_group))$Var1)
tab <- as.data.frame(table(propeller_geno$geno_group))
rownames(tab) <- tab$Var1

for (a in rownames(tab)){
  
  for (r in rownames(propeller_geno)) {
    
    if (propeller_geno[r,"geno_group"] == a) {
      propeller_geno[r,"geno_group"] <- paste0(a, " ", "(N=", tab[a,"Freq"], ")")
    }
  }
}

propeller_geno$propeller_geno <- recode(propeller_geno$geno_group, "PWK/* (N=15)"="PWK",
                                        "other (N=65)" = "other")

#define group as a vector (order/name of sample stays the same as it was listed in the df)
group <- propeller_geno[,"propeller_geno"]

#define pool from which the samples were sequenced as a vector
pool <- propeller_geno$pool

#define sex of the samples as a vector
sex <- propeller_geno$sex

#Please note that the way that `propeller` has been designed is such that the 
#group information is *always* first in the design matrix specification, and 
#there is NO intercept. this must be a design matrix with rows corresponding to 
#samples and columns to coefficients to be estimated

design <- model.matrix(~ 0 + group + pool + sex)
rownames(design) <- propeller_geno$mouse

#check
head(propeller_geno[,c("mouse", "geno_group")])
head(design)

#make contrasts
mycontr <- makeContrasts(groupPWK-groupother, levels=design)


Pkm_ttest_results <- propeller.ttest(props, 
                                     design, 
                                     contrasts = mycontr,
                                     robust=TRUE, 
                                     trend=FALSE, 
                                     sort=TRUE)

cols_Pkm_PWK <- c("PWK/* (N=15)" = "#CC5801", 
                  "other (N=65)" = "grey")


Pkm_LMP_box_PWK <- ggplot(propeller_geno, 
                          aes(x = reorder(geno_group, LMP),
                              y = LMP,
                              group = geno_group,
                              fill = geno_group)) + 
  scale_fill_manual(breaks = propeller_geno$geno_group,
                    values = cols_Pkm_PWK) +
  geom_boxplot() +
  geom_jitter() +
  xlab(label = "DO Mouse") + 
  ylab(label = "LMP Proportion") + 
  annotate(geom= "text", 
           x = 1.5, 
           y = (max(propeller_geno$LMP) + 0.075), 
           label = paste("p =", round(Pkm_ttest_results["LMP", "P.Value"], digits = 3)),
           size = 5, 
           color = "black") +
  annotate("segment", 
           x = 1, 
           xend = 2, 
           y = (max(propeller_geno$LMP) + 0.045), 
           yend = (max(propeller_geno$LMP) + 0.045)) + 
  theme(#plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.y= element_text(size = 20),
    axis.text.y=element_text(size = 15),
    legend.position = "none")
Pkm_LMP_box_PWK

Pkm_OB1_box_PWK <- ggplot(propeller_geno, 
                          aes(x = reorder(geno_group, OB1),
                              y = OB1,
                              group = geno_group,
                              fill = geno_group)) + 
  scale_fill_manual(breaks = propeller_geno$geno_group,
                    values = cols_Pkm_PWK) +
  geom_boxplot() +
  geom_jitter() +
  xlab(label = "DO Mouse") + 
  ylab(label = "OB1 Proportion") + 
  annotate(geom= "text", 
           x = 1.5, 
           y = (max(propeller_geno$OB1) + 0.075), 
           label = paste("p =", round(Pkm_ttest_results["OB1", "P.Value"], digits = 3)),
           size = 5, 
           color = "red") +
  annotate("segment", 
           x = 1, 
           xend = 2, 
           y = (max(propeller_geno$OB1) + 0.045), 
           yend = (max(propeller_geno$OB1) + 0.045)) + 
  theme(#plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.y= element_text(size = 20),
    axis.text.y=element_text(size = 15),
    legend.position = "none")
Pkm_OB1_box_PWK

Pkm_Ocy_box_PWK <- ggplot(propeller_geno, 
                          aes(x = reorder(geno_group, Ocy),
                              y = Ocy,
                              group = geno_group,
                              fill = geno_group)) + 
  scale_fill_manual(breaks = propeller_geno$geno_group,
                    values = cols_Pkm_PWK) +
  geom_boxplot() +
  geom_jitter() +
  xlab(label = "DO Mouse") + 
  ylab(label = "Ocy Proportion") + 
  annotate(geom= "text", 
           x = 1.5, 
           y = (max(propeller_geno$Ocy) + 0.075), 
           label = paste("p =", round(Pkm_ttest_results["Ocy", "P.Value"], digits = 3)),
           size = 5, 
           color = "red") +
  annotate("segment", 
           x = 1, 
           xend = 2, 
           y = (max(propeller_geno$Ocy) + 0.045), 
           yend = (max(propeller_geno$Ocy) + 0.045)) + 
  theme(#plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.y= element_text(size = 20),
    axis.text.y=element_text(size = 15),
    legend.position = "none")
Pkm_Ocy_box_PWK

grid.arrange(Pkm_LMP_box_PWK, Pkm_OB1_box_PWK, Pkm_Ocy_box_PWK, nrow = 1)

############################################################                  ############################################################
############################################################  Fig 4B S100a1   ############################################################ 
############################################################                  ############################################################

cell_type_to_examine <- "OBP"
gene <- "S100a1"
chr <- 3
pos <- 90.441235

coef <- scan1blup(genoprobs = OBP_subset_apr[,chr], 
                  pheno =  OBP_subset_cross$scRNA_seq_matrix[,gene], 
                  kinship = OBP_subset_k_loco[[chr]],
                  addcovar = OBP_covar)

plot_coefCC(x = coef, 
            map = OBP_subset_cross$pmap, 
            lodcolumn = gene, 
            scan1_output = OBP_scan1_object,
            main = paste(cell_type_to_examine, "eQTL for", gene, sep = " "),
            legend = NULL)

locus_genotypes <- maxmarg(geno_probs, 
                           OBP_subset_cross$pmap, 
                           chr = chr, 
                           pos= pos, 
                           return_char=TRUE, 
                           tol = 0.5)

############################################################                  ############################################################
############################################################  Fig 4C S100a1   ############################################################  
############################################################                  ############################################################

plot_pxg(locus_genotypes, OBP_subset_cross$scRNA_seq_matrix[,gene], ylab = paste(gene, "Expression"))

############################################################                  ############################################################
############################################################  Fig 4D S100a1   ############################################################
############################################################                  ############################################################

props <- getTransformedProps(clusters = do80.combined_all_clusters$celltype_annotation_props, 
                             sample = do80.combined_all_clusters$mouse, 
                             transform="asin")

mouse_props <- as.data.frame.matrix(t(props$TransformedProps))

DO_table <- data.frame(letter = c("A", "B", "C", "D", "E", "F", "G", "H"),
                       founder = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"))

interest <- "C"

locus_genotypes <- as.data.frame(locus_genotypes)
locus_genotypes <- na.omit(locus_genotypes)
nrow(locus_genotypes)

propeller_geno <- merge(locus_genotypes, mouse_props, by = "row.names", all.x = T)
propeller_geno$mouse <- propeller_geno$Row.names
propeller_geno <- subset(propeller_geno, select = -c(Row.names))
head(propeller_geno)

#add M/F and pool data
meta <- as.data.frame(do80.combined_all_clusters@meta.data)
deduped.data <- unique(meta[ , c("mouse", "pool", "sex")] )
#check
nrow(deduped.data)
#make DO80 pool 0
deduped.data[is.na(deduped.data)] <- "other"

#get props and geno in the same df
propeller_geno <- merge(propeller_geno, deduped.data, by = "mouse", all.x = T)
nrow(propeller_geno)
head(propeller_geno)

#add columns to df
propeller_geno$pool <- paste0("pool_", propeller_geno$pool)
propeller_geno$pool <- as.factor(propeller_geno$pool)
propeller_geno$sex <- as.factor(propeller_geno$sex)

#add new column 
propeller_geno$geno_group <- NA


for (r in rownames(propeller_geno)) {
  if (str_detect(propeller_geno[r,"locus_genotypes"], paste(strsplit(interest[1], "")[[1]][1])) == TRUE) {
    propeller_geno[r,"geno_group"] <- paste( (subset(DO_table, subset = letter == paste(strsplit(interest[1], "")[[1]][1])))[1,2],
                                             "*",
                                             sep = "/" )
  } else 
    propeller_geno[r,"geno_group"] <- "other"
  
}

head(propeller_geno)
table(propeller_geno$geno_group)
alleles <- levels(as.data.frame(table(propeller_geno$geno_group))$Var1)
tab <- as.data.frame(table(propeller_geno$geno_group))
rownames(tab) <- tab$Var1

for (a in rownames(tab)){
  
  for (r in rownames(propeller_geno)) {
    
    if (propeller_geno[r,"geno_group"] == a) {
      propeller_geno[r,"geno_group"] <- paste0(a, " ", "(N=", tab[a,"Freq"], ")")
    }
  }
}

propeller_geno$propeller_geno <- recode(propeller_geno$geno_group, "129/* (N=30)"=	"129",
                                        "other (N=50)" = "other")

#define group as a vector (order/name of sample stays the same as it was listed in the df)
group <- propeller_geno[,"propeller_geno"]

#define pool from which the samples were sequenced as a vector
pool <- propeller_geno$pool

#define sex of the samples as a vector
sex <- propeller_geno$sex

#Please note that the way that `propeller` has been designed is such that the 
#group information is *always* first in the design matrix specification, and 
#there is NO intercept. this must be a design matrix with rows corresponding to 
#samples and columns to coefficients to be estimated

design <- model.matrix(~ 0 + group + pool + sex)
rownames(design) <- propeller_geno$mouse

#check
head(propeller_geno[,c("mouse", "geno_group")])
head(design)

#make contrasts
mycontr <- makeContrasts(group129-groupother, levels=design)


S100a1_ttest_results_129 <- propeller.ttest(props, 
                                           design, 
                                           contrasts = mycontr,
                                           robust=TRUE, 
                                           trend=FALSE, 
                                           sort=TRUE)

cols_S100a1_129 <- c("129/* (N=30)" = "#FCAE1E", 
                     "other (N=50)" = "grey")

S100a1_LMP_box_129 <- ggplot(propeller_geno, 
                             aes(x = reorder(geno_group, LMP),
                                 y = LMP,
                                 group = geno_group,
                                 fill = geno_group)) + 
  scale_fill_manual(breaks = propeller_geno$geno_group,
                    values = cols_S100a1_129) +
  geom_boxplot() +
  geom_jitter() +
  xlab(label = "DO Mouse") + 
  ylab(label = "LMP Proportion") + 
  annotate(geom= "text", 
           x = 1.5, 
           y = (max(propeller_geno$LMP) + 0.075), 
           label = paste("p =", round(S100a1_ttest_results_129["LMP", "P.Value"], digits = 3)),
           size = 5, 
           color = "red") +
  annotate("segment", 
           x = 1, 
           xend = 2, 
           y = (max(propeller_geno$LMP) + 0.045), 
           yend = (max(propeller_geno$LMP) + 0.045)) + 
  theme(#plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.y= element_text(size = 20),
    axis.text.y=element_text(size = 15),
    legend.position = "none")
S100a1_LMP_box_129

S100a1_OBP_box_129 <- ggplot(propeller_geno, 
                             aes(x = reorder(geno_group, OBP),
                                 y = OBP,
                                 group = geno_group,
                                 fill = geno_group)) + 
  scale_fill_manual(breaks = propeller_geno$geno_group,
                    values = cols_S100a1_129) +
  geom_boxplot() +
  geom_jitter() +
  xlab(label = "DO Mouse") + 
  ylab(label = "OBP Proportion") + 
  annotate(geom= "text", 
           x = 1.5, 
           y = (max(propeller_geno$OBP) + 0.075), 
           label = paste("p =", round(S100a1_ttest_results_129["OBP", "P.Value"], digits = 3)),
           size = 5, 
           color = "black") +
  annotate("segment", 
           x = 1, 
           xend = 2, 
           y = (max(propeller_geno$OBP) + 0.045), 
           yend = (max(propeller_geno$OBP) + 0.045)) + 
  theme(#plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.y= element_text(size = 20),
    axis.text.y=element_text(size = 15),
    legend.position = "none")
S100a1_OBP_box_129

S100a1_OB1_box_129 <- ggplot(propeller_geno, 
                             aes(x = reorder(geno_group, OB1),
                                 y = OB1,
                                 group = geno_group,
                                 fill = geno_group)) + 
  scale_fill_manual(breaks = propeller_geno$geno_group,
                    values = cols_S100a1_129) +
  geom_boxplot() +
  geom_jitter() +
  xlab(label = "DO Mouse") + 
  ylab(label = "OB1 Proportion") + 
  annotate(geom= "text", 
           x = 1.5, 
           y = (max(propeller_geno$OB1) + 0.075), 
           label = paste("p =", round(S100a1_ttest_results_129["OB1", "P.Value"], digits = 3)),
           size = 5, 
           color = "red") +
  annotate("segment", 
           x = 1, 
           xend = 2, 
           y = (max(propeller_geno$OB1) + 0.045), 
           yend = (max(propeller_geno$OB1) + 0.045)) + 
  theme(#plot.title = element_text(hjust = 0.5, size = 20, face = "italic"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.y= element_text(size = 20),
    axis.text.y=element_text(size = 15),
    legend.position = "none")
S100a1_OB1_box_129

grid.arrange(S100a1_LMP_box_129, S100a1_OBP_box_129, S100a1_OB1_box_129, nrow = 1)

############################################################
############################################################
############################################################

#average expression matrices for each cluster were converted to .tsv files and used as input to the iterativeWGCNA analysis, which is run via a Docker container

#read in all of the iterativeWGCNA data for each cluster.
#MPC_it <- read.delim(file = "./MPC_merged-0.05-membership.txt", header = T, sep = "\t")
#LMP_it  <- read.delim(file = "./LMP_merged-0.05-membership.txt", header = T, sep = "\t")
#OBP_it  <- read.delim(file = "./OBP_merged-0.05-membership.txt", header = T, sep = "\t")
#OB1_it  <- read.delim(file = "./OB1_merged-0.05-membership.txt", header = T, sep = "\t")
#OB2_it  <- read.delim(file = "./OB2_merged-0.05-membership.txt", header = T, sep = "\t")
#Ocy_it  <- read.delim(file = "./Ocy_merged-0.05-membership.txt", header = T, sep = "\t")
#MALP_it  <- read.delim(file = "./MALP_merged-0.05-membership.txt", header = T, sep = "\t")

all_its <- c(names(mget(ls(pattern = "_it"))))

for (i in all_its) {
  data <- get(i)
  #subset the data to include only those genes that were assigned to a module
  data_sub <- subset(data, subset = Module != "UNCLASSIFIED")
  df <- data.frame(
    celltype = str_split(i, "_")[[1]][1],
    gene_id = data_sub$Gene,
    colors = data_sub$Module)
  assign(paste("module_df", i, sep = "_"), df)
}

mods <- c(names(mget(ls(pattern = "^module_df_"))))

module_df_all = data.frame()

for (m in mods) {
  data <- get(m)
  module_df_all <- rbind(module_df_all, data)
}

rownames(module_df_all) <- NULL

#for every cell type and every module generated for the cell type, generate the bn graphs.
celltypes <- unique(module_df_all$celltype)

for (c in celltypes) {
  #subset to include only the celltype
  celltype_sub <- subset(module_df_all, subset = celltype == c)
  #get vector of all modules
  modules <- unique(celltype_sub$colors)
  #get the appropriate celltype matrix
  matrix <- get(paste("bound", c, "rounded_filtered_CPM_integer_vst_QN", sep = "_"))
  
  for (m in modules) {
    #subset to include only the module 
    celltype_module_sub <- subset(celltype_sub, subset = colors == m)
    #subset the celltype count matrix to only include the genes within the module
    bn_exp <- matrix[rownames(matrix) %in% celltype_module_sub$gene_id, ]
    #transpose such that the genes are the columns are the samples are the rows
    bn_exp <- t(bn_exp)
    #bn learn using the mmhc model
    bn = mmhc(as.data.frame(bn_exp))
    #rename
    assign(paste(c, "module", m, "bn", sep ="_"), bn)
    #save
    save(bn, file = paste0("/bn_outputs/",
                           paste(c),
                           "_module_", 
                           paste(m),
                           "_bn",
                           ".Rdata"))
  }
}

############################################################
############################################################
############################################################

#create a dataframe with genes and their neighborhoods, degrees, number of bone genes in neighborhood, etc. to start network analysis

#gather df output from tradeseq analysis
#make subsets of the df for use in network analysis

#lin1
MPC_to_LMP_tradeseqs <- subset(all_lin1_specific_tradeseq_polished_df, subset = pseudotime_boundary == "MPC_to_LMP")
LMP_to_OBP_tradeseqs <- subset(all_lin1_specific_tradeseq_polished_df, subset = pseudotime_boundary == "LMP_to_OBP")
OBP_to_OB1_tradeseqs <- subset(all_lin1_specific_tradeseq_polished_df, subset = pseudotime_boundary == "OBP_to_OB1")
OB1_to_Ocy_tradeseqs <- subset(all_lin1_specific_tradeseq_polished_df, subset = pseudotime_boundary == "OB1_to_Ocy")
Ocy_to_end_tradeseqs <- subset(all_lin1_specific_tradeseq_polished_df, subset = pseudotime_boundary == "Ocy_to_end")
#lin2
OBP_to_OB2_tradeseqs <- subset(all_lin2_specific_tradeseq_polished_df, subset = pseudotime_boundary == "OBP_to_OB2")
OB2_to_end_tradeseqs <- subset(all_lin2_specific_tradeseq_polished_df, subset = pseudotime_boundary == "OB2_to_end")
#lin3
LMP_to_MALP_tradeseqs <- subset(all_lin3_specific_tradeseq_polished_df, subset = pseudotime_boundary == "LMP_to_MALP")
MALP_to_end_tradeseqs <- subset(all_lin3_specific_tradeseq_polished_df, subset = pseudotime_boundary == "MALP_to_end")

trades <- c(names(mget(ls(pattern = "_tradeseqs"))))

#load all the bns
file_names <- list.files(path= "./bn_outputs", 
                         pattern="*.Rdata", full.names=TRUE)

for (f in file_names) {
  load(paste(f))
  assign(paste(str_extract(paste(str_split(f, pattern = "/")[[1]][9]),  "[^.]+")), bn)
}

bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}

bns <- c(names(mget(ls(pattern = "_bn"))))

for (t in trades) {
  #get the celltype
  celltype <- strsplit(t,"_")[[1]][1]
  #get the data
  trade <- get(t) 
  #get only bns of that cell type
  bns_interest <- str_subset(bns, celltype)
  
  out = list()
  counter = 1
  
  for(net in bns_interest) {
    print(counter)
    #get the module number that was used to generate this specific bn
    color = str_c(strsplit(net, "_")[[1]][3:5], collapse = "_")
    z = bn2igraph(get(net))
    #subset the large df that contains the gene, module, and cell type information for all iterativeWGCNA modules
    #subset based on celltype and then the module being examined here
    mod_genes = subset(module_df_all, subset = celltype == strsplit(net,"_")[[1]][1] & 
                         colors == color)
    #get vector of the genes within the module examined here
    mod_genes = mod_genes$gene_id
    #make an empty df with the number of rows being the length of the genes within the module being examined
    kda = as.data.frame(matrix(nrow=length(mod_genes),ncol = 8))
    
    i=0
    
    for (gene in mod_genes) {
      #get all genes except the gene being examined
      mod_genes_cur = mod_genes[-which(mod_genes == gene)]
      #get all genes that are neighbors (order 3 means that its the vertex immediate neighbors, and 
      #then their immediate neighbors, and then theirs too) to the gene being examined
      neighbors = names(unlist(neighborhood(z, nodes = gene, order = 3)))
      #get the degree of connection of all genes in neighborhood, except the gene being examined
      degree = length(names(unlist(neighborhood(z, nodes = gene, order = 1))))-1
      #how many neighbors?
      num_neib = length(neighbors)
      #list of neighbors
      neib_list <- toString(neighbors)
      #how many neighbors are tradeseq genes?
      num_trade_neib = length(which(tolower(neighbors) %in% tolower(trade$gene)))
      #how many neighbors are NOT tradeseq genes?
      in_neib_not_trade = num_neib - num_trade_neib
      #of all the genes, how many are tradeseq genes?
      not_neib_in_trade = length(which(tolower(mod_genes_cur) %in% tolower(trade$gene) & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
      #of all the genes, how many are NOT radeseq genes?
      not_neib_not_trade = length(which(tolower(mod_genes_cur) %in% tolower(trade$gene) == FALSE & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
      num_in = length(neighbors(z, gene, mode="in"))
      num_out = length(neighbors(z, gene, mode="out"))
      
      i=i+1
      kda[i,1] = gene
      kda[i,2] = color
      kda[i,3] = num_neib
      kda[i,4] = num_trade_neib
      kda[i,5] = num_in
      kda[i,6] = degree
      kda[i,7] = kda[i,4]/kda[i,3]
      kda[i,8] = length(mod_genes_cur)
      kda[i,9] = not_neib_in_trade
      kda[i,10] = in_neib_not_trade
      kda[i,11] = not_neib_not_trade
      
      num_trade_genes_inMod = length(which(tolower(mod_genes_cur) %in% tolower(trade$gene)))
      kda[i,12] = num_trade_genes_inMod
      kda[i,13] = num_out
      kda[i,14] = neib_list
      
    }
    
    colnames(kda) = c("gene","color","num_neib","num_trade_neib","num_in_arcs","degree","ratio_trade_neib",
                      "num_genes_inMod", "not_neib_in_trade","in_neib_not_trade","not_neib_not_trade","num_trade_genes_inMod", "num_out", "neib_list")
    
    out[[counter]] = kda
    names(out)[counter] = color
    counter = counter+1
  }
  
  kda_full = out
  #save
  assign(paste(t, "kda_full", sep = "_"), kda_full)
  
}

kdas <- c(names(mget(ls(pattern = "_kda_full$"))))

for(k in kdas) {
  
  kda_full = get(k)
  celltype <- strsplit(k,"_")[[1]][1]
  test <- get(sub("_kda_full", "", k))
  
  for(i in 1:length(kda_full)) { 
    #candidate driver: neib > mu-bar + sd(mu)
    kda_full[[i]]$driver = 0
    kda_full[[i]]$driver[which(kda_full[[i]]$num_neib > mean(kda_full[[i]]$num_neib) + sd(kda_full[[i]]$num_neib))] = 1
    #2 sd
    kda_full[[i]]$driver2 = 0
    kda_full[[i]]$driver2[which(kda_full[[i]]$num_neib > mean(kda_full[[i]]$num_neib) + 2*(sd(kda_full[[i]]$num_neib)))] = 1
    #hub genes (d-bar + 2sd(d))
    #d = out-degree
    kda_full[[i]]$hub = 0
    kda_full[[i]]$hub[which(kda_full[[i]]$num_out > (mean(kda_full[[i]]$num_out) + 2*(sd(kda_full[[i]]$num_out))))] = 1
  }
  
  zhang = bind_rows(kda_full)
  kda_full = zhang
  all = kda_full
  all = all[-which(all$num_neib<=2),] #remove unconnected genes or those connected to only 1 
  thresh = mean(all$num_neib) - sd(all$num_neib) #threshold based on all networks combined mean neighborhoods 
  all = all[-which(all$num_neib < thresh),]
  
  #hypergeometric test to define key drivers
  
  for(i in 1:nrow(all)) {
    all$hyper[i] = phyper(q = all$num_trade_neib[i]-1, # How many tradeseq genes are in the BN? Minus one here because P(observed less than the number of tradeseq genes) and lower.tail = F because of it
                          m = length(test$gene), # What is the total number of possible tradeseq genes for this cluster?
                          n = nrow(zhang) - length(test$gene),  # What is the total number of genes for this cluster that are NOT tradeseq genes? (total number of genes across all modules made for this cluster - the total number of tradeseq genes for this cluster)
                          k = all$num_neib[i], # How many genes are in the BN? 
                          lower.tail = FALSE)
  }
  
  # FDR correction
  all$hyper_bonf = NA
  all$hyper_bonf = p.adjust(all$hyper, method="fdr")
  
  #save
  assign(paste(k, "zhang", sep = "_"), all)
}

############################################################
############################################################
############################################################

#get IMPC genes
##open conn to IMPC "experiment" core
impc = SolrClient$new(host = "www.ebi.ac.uk",path="/mi/impc/solr/experiment/select", scheme = "https",errors = "complete",port=NULL)

#query core for all "experimental" mice with BMD parameter name, and fact by gene symbol.  
#faceting gives counts of occurance of gene symbol, with a mincount of 14, which is the minimum (7 males + 7 females) 
#for DXA phenotyping protocol
res = impc$facet(params = list(q="parameter_name:Bone*Mineral*Density*", 
                               facet.field='gene_symbol',
                               facet.limit=-1, 
                               facet.sort="count", 
                               facet.mincount=14,rows=0),
                 progress = httr::progress())

#these are the genes with experimental mice for BMD query term
impc_genes = (res$facet_fields$gene_symbol)
impc_genes_vec = unique(as.character(impc_genes$term))

#remove gene names with parentheses. There are two of these and were causing issues with getting the results below.
#"Gt(ROSA)26Sor" "Tg(Thy1-MAPT*P301S)2541Godt"
impc_genes_vec = impc_genes_vec[-grep(pattern = "\\(",impc_genes_vec)] 

#For each of the genes above, get stats report.
stat = SolrClient$new(host = "www.ebi.ac.uk",
                      path = "/mi/impc/solr/statistical-result/select", 
                      scheme = "https",
                      errors = "complete",
                      port=NULL)

#check
q = paste0("parameter_name:Bone*Mineral*Density* AND marker_symbol:", as.character(impc_genes_vec[2]))
res_stat = stat$search(params = list(q=q, rows=-1), progress = httr::progress())

#get stats output for each gene and put in a df
out = as.data.frame(matrix(ncol=92,nrow=0))
colnames(out) = colnames(res_stat)

for(i in 1:length(impc_genes_vec)){
  geneName =  impc_genes_vec[i]
  q = paste0("parameter_name:Bone*Mineral*Density* AND marker_symbol:",as.character(geneName))
  res_stat = as.data.frame(stat$search(params = list(q=q, rows=-1)))
  out = merge(out, res_stat, all.x = T, all.y = T)
  print(i)
}

#last accessed May. 28, 2023
#keep only "successful" analyses
out = out[which(out$status == "Successful"),]
#some are just "not processed".

#remove BMD excluding skull
table(out$parameter_name)
out = out[-grep(pattern = "including skull", x = out$parameter_name),]

#NOTE: SPTBN1 is NOT PROCESSED. WHEN LOOKING, THERE ARE NO CONTROL MICE AND NO WEIGHTS FOR THE EXPERIMENTAL MICE
#genes with significant terms at 5e-2 pval
out_5e2 = out[which(as.numeric(out$male_ko_effect_p_value)<=0.05 | as.numeric(out$female_ko_effect_p_value)<=0.05 | as.numeric(out$genotype_effect_p_value)<= 0.05),]
#write.csv(out_5e2, file = "./impc_results_raw_processed_noskull_0.05_pvalue.csv")

############################################################
############################################################
############################################################

#using genes from Arby Abood et al. (2023), Basel Al-Barghouthi et al. (2022), IMPC genes, and eQTLs identified from the 80 DO mice
#abood <- read.csv("./abood_sQTL_MGI.csv")
#basel <- read.csv("./basel_TWASeQTL_MGI.csv")
#impc <- read.csv("./impc_results_raw_processed_noskull_0.05_pvalue.csv")
#eQTL_all_clusters <- read.csv("./eQTL_all_clusters.csv")

zhangs <- c(names(mget(ls(pattern = "_zhang"))))

for (z in zhangs) {
  data <- get(z)
  cluster <- strsplit(z,"_")[[1]][1]
  ts_genes <- get(sub("_kda_full_zhang", "", z))
  eQTL_df_sub <- subset(eQTL_all_clusters, subset = celltype == cluster)
  data$celltype <- strsplit(z,"_")[[1]][1]
  data$tradeseq_gene <- data$gene %in% c(ts_genes$gene)
  data$basel_TWAS_gene <- data$gene %in% c(unique(basel$MGI.symbol))
  data$abood_sQTL_gene <- data$gene %in% c(unique(abood$MGI.symbol))
  data$IMPC_BMD_gene <- data$gene %in% c(unique(impc$marker_symbol))
  data$DO_eQTL_gene <- data$gene %in% c(eQTL_df_sub$lodcolumn)
  
  assign(paste(z, "info", sep = "_"), data)
}

infos <- c(names(mget(ls(pattern = "_info"))))

for (i in infos) {
  data <- get(i)
  cluster <- strsplit(i,"_")[[1]][1]
  eQTL_df_sub <- subset(eQTL_all_clusters, subset = celltype == cluster)
  ts_genes <- get(sub("_kda_full_zhang_info", "", i))
  
  for (r in 1:nrow(data)) {
    
    data[r,"neib_basel_TWAS_gene"] <- toString(Reduce(intersect, list(str_split(str_remove_all(data[r, "neib_list"], ","), " ")[[1]],
                                                                      c(unique(basel$MGI.symbol)))))
    data[r,"neib_abood_sQTL_gene"] <- toString(Reduce(intersect, list(str_split(str_remove_all(data[r, "neib_list"], ","), " ")[[1]],
                                                                      c(unique(abood$MGI.symbol)))))
    data[r,"neib_IMPC_BMD_gene"] <- toString(Reduce(intersect, list(str_split(str_remove_all(data[r, "neib_list"], ","), " ")[[1]],
                                                                    c(unique(impc$marker_symbol)))))
    data[r,"neib_DO_eQTL_gene"] <- toString(Reduce(intersect, list(str_split(str_remove_all(data[r, "neib_list"], ","), " ")[[1]],
                                                                   c(unique(eQTL_df_sub$lodcolumn)))))
    data[r,"neib_tradeseq_gene"] <- toString(Reduce(intersect, list(str_split(str_remove_all(data[r, "neib_list"], ","), " ")[[1]],
                                                                    c(ts_genes$gene))))
  }
  
  assign(paste(i, "polished", sep = "_"), data)
}

polishes <- c(names(mget(ls(pattern = "_info_polished"))))

for (p in polishes) {
  data <- get(p)
  data_sub <- subset(data, subset = driver == 1 &
                       hyper <= 0.05)
  assign(paste(p, "filter", sep = "_"), data_sub)
  rm(data, data_sub)
}

all_filter_lin1_dfs <- rbind(MPC_to_LMP_tradeseqs_kda_full_zhang_info_polished_filter,
                              LMP_to_OBP_tradeseqs_kda_full_zhang_info_polished_filter,
                              OBP_to_OB1_tradeseqs_kda_full_zhang_info_polished_filter,
                              OB1_to_Ocy_tradeseqs_kda_full_zhang_info_polished_filter,
                              Ocy_to_end_tradeseqs_kda_full_zhang_info_polished_filter) %>% subset(., subset = hyper_bonf < 0.05)

all_filter_lin2_dfs <- rbind(OBP_to_OB2_tradeseqs_kda_full_zhang_info_polished_filter,
                              OB2_to_end_tradeseqs_kda_full_zhang_info_polished_filter) %>% subset(., subset = hyper_bonf < 0.05)

all_filter_lin3_dfs <- rbind(LMP_to_MALP_tradeseqs_kda_full_zhang_info_polished_filter,
                              MALP_to_end_tradeseqs_kda_full_zhang_info_polished_filter) %>% subset(., subset = hyper_bonf < 0.05)

############################################################          ############################################################
############################################################  Fig 5A  ############################################################
############################################################          ############################################################

bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}

x = bn2igraph(LMP_module_P1_I22_M20_bn)
subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,3,nodes = "Fgfrl1"))))
plot(subgraph,
     vertex.label.cex=0.9,
     vertex.size=5,
     vertex.color="orange",
     vertex.label.color= "blue",
     margin=-0.2, 
     vertex.label.dist=2.25, 
     vertex.label.degree=-pi,
     edge.arrow.size = 0.25,
     edge.width=1, 
     layout = layout.fruchterman.reingold)

############################################################          ############################################################
############################################################  Fig 5B  ############################################################
############################################################          ############################################################

fgfrl1_genes <- c("Igfbp4", "Mark1","Cyp1b1", "St5","Hnmt", "Pdzrn4")
fgfrl1_colors <- c("#ff4500", "#2f4f4f", "#0000cd", "#8b0000", "#ff00ff","#00ffff")

for (i in 1:length(fgfrl1_genes)) {
  
  models = gamFIT_10knots
  counts = countsFiltered
  gene = fgfrl1_genes[i]
  nPoints = 100
  lwd = 2
  size = 2/3 
  xlab = "Pseudotime"
  ylab =  paste(gene, "Log(expression + 1)")
  border = FALSE
  alpha = 2/3 
  sample = 1 
  pointCol = NULL
  curvesCols = NULL
  plotLineages = TRUE 
  lineagesToPlot = c(3)
  linecolor = fgfrl1_colors[i]
  
  # Predicting fits ----
  # lpmatrix given X and design
  predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL){
    # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
    # INPUT:
    # lpmatrix is the linear predictor matrix of the GAM model
    # df is a data frame of values for which we want the lpmatrix
    # pseudotime is the n x l matrix of pseudotimes
    # conditions is the vector of conditions, if present.
    
    # if pseudotime is vector, make it a matrix.
    if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)
    
    condPresent <- !is.null(conditions)
    if(condPresent) nConditions <- nlevels(conditions)
    
    # for each curve, specify basis function IDs for lpmatrix
    allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")
    
    if(!condPresent){
      lineages <- sub(pattern = "s\\(", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
      lineages <- sub(pattern = "\\):.*", replacement = "",
                      x = lineages)
      nCurves <- length(unique(lineages))
      for (ii in seq_len(nCurves)) {
        assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
      }
    } else if(condPresent){
      lineages <- sub(pattern = "s\\(t", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
      lineages <- sub(pattern = "\\):.*", replacement = "",
                      x = lineages)
      nLineages <- length(unique(lineages))
      curves <- sub(pattern = ".*:l", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
      curves <- sub(pattern = "\\..*", replacement = "",
                    x = curves)
      nCurves <- length(unique(curves))
      for (ii in seq_len(nLineages)) {
        for(kk in seq_len(nConditions))
          assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
      }
    }
    
    
    # specify lineage assignment for each cell (i.e., row of lpmatrix)
    if(!condPresent){
      lineageID <- apply(lpmatrix, 1, function(x){
        for (ii in seq_len(nCurves)) {
          if (!all(x[get(paste0("id", ii))] == 0)) {
            return(ii)
          }
        }
      })
    } else if(condPresent){
      # first number is lineage, second number is condition.
      lineageID <- apply(lpmatrix, 1, function(x){
        for (ii in seq_len(nLineages)) {
          # loop over lineages
          for(kk in seq_len(nConditions)){
            # loop over conditions
            if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
              return(as.numeric(paste0(ii, kk)))
            }
          }
        }
      })
    }
    
    
    # fit splinefun for each basis function based on assigned cells
    if(!condPresent) {
      for (ii in seq_len(nCurves)) { # loop over curves
        for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
          assign(paste0("l",ii,".",jj),
                 stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                  y = lpmatrix[lineageID == ii, #only cells for lineage
                                               get(paste0("id", ii))[jj]],
                                  ties = mean)) #basis function
        }
      }
    } else if(condPresent) {
      for (ii in  seq_len(nLineages)) {
        # loop over curves
        for(kk in seq_len(nConditions)){
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
            #within curve, loop over basis functions
            assign(paste0("l",ii, "_", kk,".",jj),
                   stats::splinefun(
                     x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                     y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                  get(paste0("id", ii, "_", kk))[jj]],
                     ties = mean)) #basis function
          }
        }
      }
    }
    
    
    # use input to estimate X for each basis function
    Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
    if(!condPresent){
      for (ii in seq_len(nCurves)) { # loop over curves
        if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
            f <- get(paste0("l", ii, ".", jj))
            Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    } else if(condPresent){
      # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
      for (ii in seq_len(nLineages)) {
        # loop over curves
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
            for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
              # within curve, loop over basis functions
              f <- get(paste0("l", ii, "_", kk, ".", jj))
              Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
            }
          }
        }
      }
    }
    
    
    # add fixed covariates as in df
    dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
    dfOffsetID <- grep(x = colnames(df), pattern = "offset")
    Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]
    
    # return
    colnames(Xout) <- colnames(lpmatrix)
    return(Xout)
  }
  
  # get the first non-errored fit in models
  .getModelReference <- function(models){
    for (i in seq_len(length(models))) {
      m <- models[[i]]
      if (is(m)[1] != "try-error") return(m)
    }
    stop("All models errored")
  }
  
  .getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
    vars <- dm[1, ]
    if ("y" %in% colnames(vars)) {
      vars <- vars[!colnames(vars) %in% "y"]
      off <- 1
    } else {
      off <- 0
    }
    offsetId <- grep(x = colnames(vars), pattern = "offset")
    offsetName <- colnames(vars)[offsetId]
    offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
    names(vars)[offsetId] <- offsetName
    # set all times on 0
    vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
    # set all lineages on 0
    vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
    # duplicate to nPoints
    vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
    # set range of pseudotime for lineage of interest
    if (is.null(conditionId)) {
      lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
    } else {
      lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                          "_", conditionId, "$"))
    }
    if (length(lineageIds) == 1){
      lineageData <- dm[dm[, lineageIds + off] == 1,
                        paste0("t", lineageId)]
    } else {
      lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                        paste0("t", lineageId)]
    }
    # make sure lineage starts at zero
    if(min(lineageData) / max(lineageData) < .01) {
      lineageData[which.min(lineageData)] <- 0
    }
    vars[, lineageIds] <- 1 / length(lineageIds)
    # set lineage
    vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                          max(lineageData),
                                          length = nPoints)
    # set offset
    vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                         pattern = "offset")])
    return(vars)
  }
  
  
  #input is singleCellExperiment object.
  if (is.null(names(models))) {
    rownames(models) <- rownames(counts) <- seq_len(nrow(models))
    message(paste0(
      "The sce object has no rownames. Assuming that the counts and the sce ",
      "objects are ordered in the same way"))
  }
  if (length(gene) > 1) stop("Only provide a single gene's ID with the ",
                             "gene argument.")
  # check if all gene IDs provided are present in the models object.
  if (is(gene, "character")) {
    if (!all(gene %in% names(models))) {
      stop("The gene ID is not present in the models object.")
    }
    id <- which(names(models) %in% gene)
  } else id <- gene
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  y <- unname(counts[names(models),][id,])
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  # pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
  #                                     pattern = "pseudotime")]
  pseudotime <- pseudotimeFiltered
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id,]
  
  
  #construct time variable based on cell assignments.
  lcol <- timeAll <- rep(0, nrow(dm))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(dm))) {
      if (dm[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- dm[ii, paste0("t", jj)]
        lcol[ii] <- jj
      } else {
        next
      }
    }
  }
  
  if (!is.null(pointCol)) {
    if (length(pointCol) == 1) {
      col <- colData(models)[,pointCol]
    } else if (length(pointCol) == ncol(models)) {
      col <- pointCol
    } else {
      col <- lcol
      message(paste("pointCol should have length of either 1 or the number of cells,",
                    "reverting to default color scheme."))
    }
  } else {
    col <- lcol
  }
  
  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "pCol" = as.character(col),
                   "lineage" = as.character(lcol))
  rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
  df <- df[rows, ]
  if(!is.null(lineagesToPlot)){
    df <- df[df$lineage %in% lineagesToPlot,]
  }
  p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if(is.null(pointCol)){
    p <- p #+
    #geom_point(size = size, aes(col = lineage), colour="#ffffff") +
    #add this to the geom_point parameters if necessary: colour="#ffffff"
    #scale_color_viridis_d(alpha = alpha)
  } else {
    p <- p #+
    #geom_point(size = size, alpha = alpha, aes(col = pCol)) +
    #scale_color_discrete() +
    #labs(col = "Cell labels")
  }
  
  
  
  # predict and plot smoothers across the range
  if (plotLineages) {
    if (!is.null(curvesCols)) {
      if (length(curvesCols) != nCurves) {
        curvesCols <- viridis::viridis(nCurves)
        message("Incorrect number of lineage colors. Default to viridis")
      }
    } else {
      curvesCols <- viridis::viridis(nCurves)
    }
    if(is.null(lineagesToPlot)){
      lineagesToPlot <- seq_len(nCurves)
    }
    for (jj in lineagesToPlot) {
      df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime)
      yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
      if (border) {
        p <- p +
          geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                      "gene_count" = yhat,
                                      "lineage" = as.character(jj),
                                      "pCol" = as.character(jj)),
                    lwd = lwd + 1, colour = "white")
      }
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj),
                                    "pCol" = as.character(jj)),
                  lwd = lwd, col = linecolor)
      #curvesCols[jj])
    }
  }
  
  assign(paste(gene, "smooth_plot", sep = "_"), p)
  
  final_plot <- p +
    #ggtitle(paste(gene)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    #scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2", "Adipogenic"), values= c("#552586","#009999","#FFFF33")) +
    #scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2"), values= c("#552586","#009999")) +
    #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 1"), values= c("#552586")) +
    #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 2"), values= c("#009999")) +
    #scale_color_manual(name = "Cell Lineages", labels = c("Adipogenic"), values= c("#FFFF33")) + 
    scale_x_continuous(limits = c(0, 85)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 40, face = "italic"),
          axis.title.x = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.text.y = element_text(size = 25),
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length  = unit(0.25, "cm"),
          legend.position = "none",
          legend.title = element_text(size =15),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",color = NA)) + 
    
    geom_vline(xintercept = 0, color = 'red')  +
    annotate("label", x=0, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MPC", size = 7) +
    
    geom_vline(xintercept = 20, color = 'red')  +
    annotate("label", x=20, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "LMP", size = 7) +
    
    geom_vline(xintercept = 36, color = 'red')  +
    annotate("label", x=33, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MALP", size = 7)
  
  final_plot
  
  assign(paste(gene, "final_plot", sep = "_"), final_plot)
}

fgfrl1_network_genes <- 
  Hnmt_smooth_plot + 
  St5_smooth_plot[[2]] + 
  Igfbp4_smooth_plot[[2]] + 
  Cyp1b1_smooth_plot[[2]] + 
  Pdzrn4_smooth_plot[[2]] +
  Mark1_smooth_plot[[2]]

fgfrl1_network_ts_genes_graph  <- fgfrl1_network_genes + 
  theme(axis.title.x = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        #legend.position = "none",
        legend.title = element_text(size =15),
        legend.text = element_text(size = 15),
        #legend.background = element_rect(fill = "transparent"),
        #panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",color = NA)) + 
  labs(y = "Fgfrl1 Network TradeSeq Genes \n Log(expression + 1)") +
  geom_vline(xintercept = 0, color = 'red')  +
  annotate("label", x=1, y=5, label= "MPC", size = 6) +
  geom_vline(xintercept = 20, color = 'red')  +
  annotate("label", x=20, y=5, label= "LMP", size = 6) +
  geom_vline(xintercept = 36, color = 'red')  +
  annotate("label", x=37, y=5, label= "MALP", size = 6)

#graph
fgfrl1_network_ts_genes_graph 

#legend
names <- fgfrl1_genes
clrs <- fgfrl1_colors
ltype <- c(1, 1, 1, 1, 1, 1)
plot(NULL, xaxt='n', yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", title="Fgfrl1 Network\nTradeSeq Genes", legend = names, lty=ltype, lwd=5, cex=1.25,
       bty='n', col = clrs)

############################################################          ############################################################
############################################################  Fig 5C  ############################################################
############################################################          ############################################################

bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}

x = bn2igraph(LMP_module_P1_I22_M16_bn)
subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,3,nodes = "Tpx2"))))
plot(subgraph,
     vertex.label.cex=0.9,
     vertex.size=5,
     vertex.color="orange",
     vertex.label.color= "blue",
     margin=-0.2, 
     vertex.label.dist=2.25, 
     vertex.label.degree=-pi,
     edge.arrow.size = 0.25,
     edge.width=1, 
     layout = layout.fruchterman.reingold)

############################################################          ############################################################
############################################################  Fig 5D  ############################################################
############################################################          ############################################################

tpx2_genes <- c("Tpx2", "Top2a", "Kif4", "Iqgap3", "Prc1", "Kif11", "Ect2", "Sgo2a", "Ube2c")
tpx2_colors <- c("#8b0000", "#2f4f4f", "#228b22", "#ffd700", "#0000cd", "#ff00ff", "#00ffff", "#1e90ff", "#ff4500")

for (i in 1:length(tpx2_genes)) {

  models = gamFIT_10knots
  counts = countsFiltered
  gene = tpx2_genes[i]
  nPoints = 100
  lwd = 2
  size = 2/3 
  xlab = "Pseudotime"
  ylab =  paste(gene, "Log(expression + 1)")
  border = FALSE
  alpha = 2/3 
  sample = 1 
  pointCol = NULL
  curvesCols = NULL
  plotLineages = TRUE 
  lineagesToPlot = c(1)
  linecolor = tpx2_colors[i]
  
  # Predicting fits ----
  # lpmatrix given X and design
  predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL){
    # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
    # INPUT:
    # lpmatrix is the linear predictor matrix of the GAM model
    # df is a data frame of values for which we want the lpmatrix
    # pseudotime is the n x l matrix of pseudotimes
    # conditions is the vector of conditions, if present.
    
    # if pseudotime is vector, make it a matrix.
    if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)
    
    condPresent <- !is.null(conditions)
    if(condPresent) nConditions <- nlevels(conditions)
    
    # for each curve, specify basis function IDs for lpmatrix
    allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")
    
    if(!condPresent){
      lineages <- sub(pattern = "s\\(", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
      lineages <- sub(pattern = "\\):.*", replacement = "",
                      x = lineages)
      nCurves <- length(unique(lineages))
      for (ii in seq_len(nCurves)) {
        assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
      }
    } else if(condPresent){
      lineages <- sub(pattern = "s\\(t", replacement = "",
                      x = colnames(lpmatrix[,allBs]))
      lineages <- sub(pattern = "\\):.*", replacement = "",
                      x = lineages)
      nLineages <- length(unique(lineages))
      curves <- sub(pattern = ".*:l", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
      curves <- sub(pattern = "\\..*", replacement = "",
                    x = curves)
      nCurves <- length(unique(curves))
      for (ii in seq_len(nLineages)) {
        for(kk in seq_len(nConditions))
          assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
      }
    }
    
    
    # specify lineage assignment for each cell (i.e., row of lpmatrix)
    if(!condPresent){
      lineageID <- apply(lpmatrix, 1, function(x){
        for (ii in seq_len(nCurves)) {
          if (!all(x[get(paste0("id", ii))] == 0)) {
            return(ii)
          }
        }
      })
    } else if(condPresent){
      # first number is lineage, second number is condition.
      lineageID <- apply(lpmatrix, 1, function(x){
        for (ii in seq_len(nLineages)) {
          # loop over lineages
          for(kk in seq_len(nConditions)){
            # loop over conditions
            if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
              return(as.numeric(paste0(ii, kk)))
            }
          }
        }
      })
    }
    
    
    # fit splinefun for each basis function based on assigned cells
    if(!condPresent) {
      for (ii in seq_len(nCurves)) { # loop over curves
        for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
          assign(paste0("l",ii,".",jj),
                 stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                  y = lpmatrix[lineageID == ii, #only cells for lineage
                                               get(paste0("id", ii))[jj]],
                                  ties = mean)) #basis function
        }
      }
    } else if(condPresent) {
      for (ii in  seq_len(nLineages)) {
        # loop over curves
        for(kk in seq_len(nConditions)){
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
            #within curve, loop over basis functions
            assign(paste0("l",ii, "_", kk,".",jj),
                   stats::splinefun(
                     x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                     y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                  get(paste0("id", ii, "_", kk))[jj]],
                     ties = mean)) #basis function
          }
        }
      }
    }
    
    
    # use input to estimate X for each basis function
    Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
    if(!condPresent){
      for (ii in seq_len(nCurves)) { # loop over curves
        if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
            f <- get(paste0("l", ii, ".", jj))
            Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    } else if(condPresent){
      # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
      for (ii in seq_len(nLineages)) {
        # loop over curves
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
            for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
              # within curve, loop over basis functions
              f <- get(paste0("l", ii, "_", kk, ".", jj))
              Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
            }
          }
        }
      }
    }
    
    
    # add fixed covariates as in df
    dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
    dfOffsetID <- grep(x = colnames(df), pattern = "offset")
    Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]
    
    # return
    colnames(Xout) <- colnames(lpmatrix)
    return(Xout)
  }
  
  # get the first non-errored fit in models
  .getModelReference <- function(models){
    for (i in seq_len(length(models))) {
      m <- models[[i]]
      if (is(m)[1] != "try-error") return(m)
    }
    stop("All models errored")
  }
  
  .getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
    vars <- dm[1, ]
    if ("y" %in% colnames(vars)) {
      vars <- vars[!colnames(vars) %in% "y"]
      off <- 1
    } else {
      off <- 0
    }
    offsetId <- grep(x = colnames(vars), pattern = "offset")
    offsetName <- colnames(vars)[offsetId]
    offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
    names(vars)[offsetId] <- offsetName
    # set all times on 0
    vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
    # set all lineages on 0
    vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
    # duplicate to nPoints
    vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
    # set range of pseudotime for lineage of interest
    if (is.null(conditionId)) {
      lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
    } else {
      lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                          "_", conditionId, "$"))
    }
    if (length(lineageIds) == 1){
      lineageData <- dm[dm[, lineageIds + off] == 1,
                        paste0("t", lineageId)]
    } else {
      lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                        paste0("t", lineageId)]
    }
    # make sure lineage starts at zero
    if(min(lineageData) / max(lineageData) < .01) {
      lineageData[which.min(lineageData)] <- 0
    }
    vars[, lineageIds] <- 1 / length(lineageIds)
    # set lineage
    vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                          max(lineageData),
                                          length = nPoints)
    # set offset
    vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                         pattern = "offset")])
    return(vars)
  }
  
  
  #input is singleCellExperiment object.
  if (is.null(names(models))) {
    rownames(models) <- rownames(counts) <- seq_len(nrow(models))
    message(paste0(
      "The sce object has no rownames. Assuming that the counts and the sce ",
      "objects are ordered in the same way"))
  }
  if (length(gene) > 1) stop("Only provide a single gene's ID with the ",
                             "gene argument.")
  # check if all gene IDs provided are present in the models object.
  if (is(gene, "character")) {
    if (!all(gene %in% names(models))) {
      stop("The gene ID is not present in the models object.")
    }
    id <- which(names(models) %in% gene)
  } else id <- gene
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  y <- unname(counts[names(models),][id,])
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  # pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
  #                                     pattern = "pseudotime")]
  pseudotime <- pseudotimeFiltered
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id,]
  
  
  #construct time variable based on cell assignments.
  lcol <- timeAll <- rep(0, nrow(dm))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(dm))) {
      if (dm[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- dm[ii, paste0("t", jj)]
        lcol[ii] <- jj
      } else {
        next
      }
    }
  }
  
  if (!is.null(pointCol)) {
    if (length(pointCol) == 1) {
      col <- colData(models)[,pointCol]
    } else if (length(pointCol) == ncol(models)) {
      col <- pointCol
    } else {
      col <- lcol
      message(paste("pointCol should have length of either 1 or the number of cells,",
                    "reverting to default color scheme."))
    }
  } else {
    col <- lcol
  }
  
  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "pCol" = as.character(col),
                   "lineage" = as.character(lcol))
  rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
  df <- df[rows, ]
  if(!is.null(lineagesToPlot)){
    df <- df[df$lineage %in% lineagesToPlot,]
  }
  p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if(is.null(pointCol)){
    p <- p #+
    #geom_point(size = size, aes(col = lineage), colour="#ffffff") +
    #add this to the geom_point parameters if necessary: colour="#ffffff"
    #scale_color_viridis_d(alpha = alpha)
  } else {
    p <- p #+
    #geom_point(size = size, alpha = alpha, aes(col = pCol)) +
    #scale_color_discrete() +
    #labs(col = "Cell labels")
  }
  
  
  
  # predict and plot smoothers across the range
  if (plotLineages) {
    if (!is.null(curvesCols)) {
      if (length(curvesCols) != nCurves) {
        curvesCols <- viridis::viridis(nCurves)
        message("Incorrect number of lineage colors. Default to viridis")
      }
    } else {
      curvesCols <- viridis::viridis(nCurves)
    }
    if(is.null(lineagesToPlot)){
      lineagesToPlot <- seq_len(nCurves)
    }
    for (jj in lineagesToPlot) {
      df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime)
      yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
      if (border) {
        p <- p +
          geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                      "gene_count" = yhat,
                                      "lineage" = as.character(jj),
                                      "pCol" = as.character(jj)),
                    lwd = lwd + 1, colour = "white")
      }
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj),
                                    "pCol" = as.character(jj)),
                  lwd = lwd, col = linecolor)
      #curvesCols[jj])
    }
  }
  
  assign(paste(gene, "smooth_plot", sep = "_"), p)
  
  final_plot <- p +
    #ggtitle(paste(gene)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    #scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2", "Adipogenic"), values= c("#552586","#009999","#FFFF33")) +
    #scale_color_manual(name = "Pseudotime \n Trajectories", labels = c("Osteogenic 1", "Osteogenic 2"), values= c("#552586","#009999")) +
    #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 1"), values= c("#552586")) +
    #scale_color_manual(name = "Cell Lineages", labels = c("Osteogenic 2"), values= c("#009999")) +
    #scale_color_manual(name = "Cell Lineages", labels = c("Adipogenic"), values= c("#FFFF33")) + 
    scale_x_continuous(limits = c(0, 85)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 40, face = "italic"),
          axis.title.x = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.text.y = element_text(size = 25),
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length  = unit(0.25, "cm"),
          legend.position = "none",
          legend.title = element_text(size =15),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",color = NA)) + 
    
    geom_vline(xintercept = 0, color = 'red')  +
    annotate("label", x=0, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "MPC", size = 7) +
    
    geom_vline(xintercept = 20, color = 'red')  +
    annotate("label", x=20, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "LMP", size = 7) +
    
    geom_vline(xintercept=42, color = 'red') + 
    annotate("label", x=42, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "OBP", size = 7) +
    
    geom_vline(xintercept=55, color = 'red') +
    annotate("label", x=54, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.35), label= "OB1/\nOB2", size = 7) +
    
    geom_vline(xintercept=62, color = 'red') +
    annotate("label", x=62, y=(round(max(log1p(countsFiltered[gene,])), digits = 0) + 0.5), label= "Ocy", size = 7)
  
  final_plot
  
  assign(paste(gene, "final_plot", sep = "_"), final_plot)
}

tpx2_network_genes <- 
  Tpx2_smooth_plot + 
  Top2a_smooth_plot[[2]] + 
  Kif4_smooth_plot[[2]] + 
  Iqgap3_smooth_plot[[2]] + 
  Prc1_smooth_plot[[2]] +
  Kif11_smooth_plot[[2]] +
  Ect2_smooth_plot[[2]] +
  Sgo2a_smooth_plot[[2]] +
  Ube2c_smooth_plot[[2]]

tpx2_network_ts_genes_graph  <- tpx2_network_genes + 
  theme(axis.title.x = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length  = unit(0.25, "cm"),
        #legend.position = "none",
        legend.title = element_text(size =15),
        legend.text = element_text(size = 15),
        #legend.background = element_rect(fill = "transparent"),
        #panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",color = NA)) + 
  labs(y = "Tpx2 Network TradeSeq Genes \n Log(expression + 1)") +
  geom_vline(xintercept = 0, color = 'red')  +
  annotate("label", x=0, y=2, label= "MPC", size = 6) +
  geom_vline(xintercept = 20, color = 'red')  +
  annotate("label", x=20, y=2, label= "LMP", size = 6) +
  geom_vline(xintercept=42, color = 'red') + 
  annotate("label", x=42, y=2, label= "OBP", size = 6) +
  geom_vline(xintercept=55, color = 'red') +
  annotate("label", x=53, y=1.95, label= "OB1/\nOB2", size = 6) +
  geom_vline(xintercept=62, color = 'red') +
  annotate("label", x=62, y=2, label= "Ocy", size = 6)


#graph
tpx2_network_ts_genes_graph 

#legend
names <- tpx2_genes
clrs <- tpx2_colors
ltype <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
plot(NULL, xaxt='n', yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", title="Tpx2 Network\nTradeSeq Genes", legend = names, lty=ltype, lwd=5, cex=1.25,
       bty='n', col = clrs)

############################################################            ############################################################
############################################################  Fig 5BD   ############################################################
############################################################            ############################################################

grid.arrange(fgfrl1_network_ts_genes_graph, tpx2_network_ts_genes_graph, nrow = 1)

############################################################          ############################################################
############################################################  Fig 5E  ############################################################
############################################################          ############################################################

IMPC_data$sex <- recode(IMPC_data$sex, "male" = "Male", "female" = "Female")

p <- ggplot(IMPC_data, aes(x=sex, y=data_point, fill=biological_sample_group)) +
  geom_boxplot(position = position_dodge(1), outlier.colour = NA) +
  geom_point(size = 0.1,
             position = position_jitterdodge(dodge.width = 1)) +
  labs(title = expression(paste(italic("Tpx2 "), "IMPC Mutant")),
       y = expression("Bone Mineral Density (excluding skull) (g/"*cm^{2}*")"),
       x = "Sex")

p + scale_fill_manual(values=c("#999999","#56B4E9"),
                      name="Genotype",
                      breaks=c("control", "experimental"),
                      labels=c("Control", expression("Tpx2"^"em1(IMPC)J"))) +
  theme(legend.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.length  = unit(0.25, "cm"),
        legend.key.size = unit(0.25, 'cm'),
        legend.text=element_text(size=12, hjust = 0),
        legend.position = "right")

############################################################
############################################################ 
############################################################