library(data.table)

##

loom2matrix <- function(x){
  require(loomR)
  lfile <- connect(x)
  out <- t(lfile$matrix[,])
  colnames(out) <- lfile$col.attrs$CellID[]
  row.names(out) <- lfile$row.attrs$Gene[]
  lfile$close_all()
  return(out)
}

outlier.mad <- function(x, nmad=3, ...){
  xmad <- mad(x, ... )
  x < median(x) - (nmad*xmad)
}

##
gn <- fread("data/SS2_CellType_sc_run1_GeneNames.txt")
meta1 <- fread("data/SS2_CellType_sc_run1_SampleSheet.csv", skip=15)
meta2 <- fread("data/SS2_CellType_sc_run2_SampleSheet.csv", skip=15)
meta1[, sample_bc := paste0(index,index2)]
meta2[, sample_bc := paste0(index,index2)]
mat1 <- loom2matrix("data/SS2_CellType_sc_run1_Readcount_ALL.loom")
mat2 <- loom2matrix("data/SS2_CellType_sc_run2_Readcount_ALL.loom")

library(Seurat)
library(magrittr)

so.ss2.celltype <- merge(
  CreateSeuratObject(counts = mat1, meta.data = data.frame(meta1[match(colnames(mat1), sample_bc)], row.names = "sample_bc") ),
  CreateSeuratObject(counts = mat2, meta.data = data.frame(meta2[match(colnames(mat2), sample_bc)], row.names = "sample_bc") ),
  add.cell.ids = c("1", "2")
)

so.ss2.celltype <- so.ss2.celltype[,!outlier.mad(log10(so.ss2.celltype$nCount_RNA)) | !outlier.mad(log10(so.ss2.celltype$nFeature_RNA))]
so.celltype$inhibitor <- ifelse(so.celltype$inhibitor_type == "SS2", "RRI", "SEQURNA")

so.ss2.celltype %<>% NormalizeData()
so.ss2.celltype %<>% FindVariableFeatures(selection.method = "vst")
so.ss2.celltype %<>% ScaleData()
so.ss2.celltype %<>% RunPCA(verbose=F)
so.ss2.celltype %<>% RunUMAP(dims = 1:50, n.neighbors=30)
so.ss2.celltype %<>% FindNeighbors(dims = 1:50, k.param =30)
so.ss2.celltype %<>% FindClusters(resolution = 0.3)

## PLOT
library(ggplot2)
library(cowplot)
ggsave2("plots/ss2_celltype_umap.pdf", width = 12, height = 3,
  (DimPlot(so.ss2.celltype, group.by = "cell_type", cols = c("#6BAED6", "#F28B41")) + ggtitle(NULL) + theme(aspect.ratio = 1)) +
  (DimPlot(so.ss2.celltype, group.by = "inhibitor_concentration", cols = c("#A3CF99", "#31A353", "#756BAF")) + ggtitle(NULL) + theme(aspect.ratio = 1)) +
  (DimPlot(so.ss2.celltype, group.by = "seurat_clusters", cols = c("#676BB0", "#D7616B", "#B5D06B", "#C06FAB")) + ggtitle(NULL) + theme(aspect.ratio = 1))
)

ggsave2("plots/ss2_celltype_features.pdf", width = 4, height = 12,
  (FeaturePlot(so.ss2.celltype, features = gn[gene_name == "Cd79a", gene_id] ) + ggtitle("Cd79a") + scale_color_viridis_c() + theme(aspect.ratio = 1)) +
  (FeaturePlot(so.ss2.celltype, features = gn[gene_name == "Trbc1", gene_id] ) + ggtitle("Trbc1") + scale_color_viridis_c() + theme(aspect.ratio = 1)) +
  (FeaturePlot(so.ss2.celltype, features = gn[gene_name == "Ahsg", gene_id] ) + ggtitle("Ahsg") + scale_color_viridis_c() + theme(aspect.ratio = 1)) +
  (FeaturePlot(so.ss2.celltype, features = gn[gene_name == "Cd300a", gene_id] ) + ggtitle("Cd300a") + scale_color_viridis_c() + theme(aspect.ratio = 1))
)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

markers <- as.data.table(FindAllMarkers(so.ss2.celltype, only.pos=T))
top.gn <- markers[,head(gene,30), by="cluster"][, V1]
names(top.gn) <- gn[match(top.gn, gene_id), gene_name]
mat.z <- t(apply(as.matrix(so.ss2.celltype@assays$RNA@data[top.gn,]), 1, scale))

cha <- HeatmapAnnotation(Cluster = factor(so.ss2.celltype$seurat_clusters), Tissue = so.ss2.celltype$cell_type, SEQURNA = so.ss2.celltype$inhibitor_concentration, col= list(Cluster = setNames(c("#676BB0", "#D7616B", "#B5D06B", "#C06FAB"), 0:3), Tissue = setNames(c("#6BAED6", "#F28B41"), c("liver", "spleen")), SEQURNA = c("3" ="#A3CF99", "4.5" = "#31A353", "RRI" = "#756BAF")), border=T)
rha <- rowAnnotation(gn = anno_mark(at = c(1:10, 31:40, 61:70, 91:100), labels = names(top.gn[c(1:10, 31:40, 61:70, 91:100)])))

p.heat <-
  Heatmap(
    matrix = mat.z,
    name = "Z-score",
    border = "black",
    cluster_rows = T,
    show_column_names = F,
    show_row_names = F,
    col = colorRamp2(seq(-2,2,length.out = 11), rev(brewer.pal(11, "RdBu"))),
    column_split = so.ss2.celltype$seurat_clusters, 
    cluster_column_slices = F,
    show_row_dend = F,
    #column_gap = unit(0,"mm"),
    top_annotation = cha,
    right_annotation = rha,
    use_raster = T
  )

pdf("plots/ss2_celltype_heatmap.pdf", width = 6, height = 8)
  p.heat
dev.off()
