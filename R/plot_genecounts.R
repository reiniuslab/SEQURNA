library(data.table)
library(ggplot2)
library(cowplot)
library(SingleCellExperiment)

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

##

## SS2 bulk

meta.ss2.bulk <- fread("dataSS2_RNA_bulk_BarcodeAnnotations.csv")

mat.ss2.bulk <- loom2matrix("data/SS2_RNA_bulk_Readcount_Downsampled100k.loom")

sce.ss2.bulk <- SingleCellExperiment(assays = list(counts = mat.ss2.bulk), colData=meta.ss2.bulk[match(colnames(mat.ss2.bulk), sample_barcode)])

ggplot(data.table(y=colSums(assay(sce.ss2.bulk)>0), x= factor(sce.ss2.bulk$inhibitor_concentration, levels=c(0,0.075,0.15,0.3,0.75,1.5,3,4.5,6,9,12,18,24,30,"RRI")), time = sce.ss2.bulk$rt_incubation_length_hr, rna = sce.ss2.bulk$rna_concentration)[time < 2 & rna > 50], aes(x=x,y=y)) + geom_boxplot() + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS2 bulk") + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## SS2 single-cell

meta.ss2.hek1 <- fread("data/SS2_HEK_sc_run1_BarcodeAnnotations.csv")
meta.ss2.hek2 <- fread("data/SS2_HEK_sc_run2_BarcodeAnnotations.csv")

mat.ss2.hek1 <- loom2matrix("data/SS2_HEK_sc_run1_Readcounts_Downsampled100k.loom")
mat.ss2.hek2 <- loom2matrix("data/SS2_HEK_sc_run2_Readcount_Downsampled100k.loom")

idx.row.ss2.hek <- Reduce(intersect, list(row.names(mat.ss2.hek1), row.names(mat.ss2.hek2)))

sce.ss2.hek <- SingleCellExperiment(assays = list(counts = cbind(mat.ss2.hek1[idx.row.ss2.hek, ], mat.ss2.hek2[idx.row.ss2.hek,])), colData=rbindlist(list(meta.ss2.hek1[match(colnames(mat.ss2.hek1), sample_barcode)],meta.ss2.hek2[match(colnames(mat.ss2.hek2), sample_barcode)])))

ggplot(data.frame(y=colSums(assay(sce.ss2.hek)>0), x= factor(sce.ss2.hek$inhibitor_concentration, levels=c(0,0.15,0.3,0.75,3,4.5,6,9,12,18,24, "RRI"))), aes(x=x,y=y)) + geom_boxplot() + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS2 single-cell HEK") + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## SS2 cell type

meta.ss2.celltype1 <- fread("data/SS2_CellType_sc_run1_SampleSheet.csv", skip=15)
meta.ss2.celltype1[, sample_bc := paste0(index,index2)]
meta.ss2.celltype2 <- fread("data/SS2_CellType_sc_run2_SampleSheet.csv", skip=15)
meta.ss2.celltype2[, sample_bc := paste0(index,index2)]

mat.ss2.celltype1 <- loom2matrix("data/SS2_CellType_sc_run1_Readcount_Downsampled100k.loom")
mat.ss2.celltype2 <- loom2matrix("data/SS2_CellType_sc_run2_Readcount_Downsampled100k.loom")

idx.row.ss2.celltype <- Reduce(intersect, list(row.names(mat.ss2.celltype1), row.names(mat.ss2.celltype2)))

sce.ss2.celltype <- SingleCellExperiment(assays = list(counts = cbind(mat.ss2.celltype1[idx.row.ss2.celltype, ], mat.ss2.celltype2[idx.row.ss2.celltype,])), colData=rbindlist(list(meta.ss2.celltype1[match(colnames(mat.ss2.celltype1), sample_bc)],meta.ss2.celltype2[match(colnames(mat.ss2.celltype2), sample_bc)])))

ggplot(data.frame(y=colSums(assay(sce.ss2.celltype)>0), x= factor(sce.ss2.celltype$inhibitor_concentration, levels=c(3,4.5,"RRI")), cell_type = sce.ss2.celltype$cell_type), aes(x=x,y=y)) + geom_boxplot() + facet_wrap(~cell_type) + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS2 celltype") + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## SS3 bulk + single-cell

meta.ss3 <- fread("data/SS3_HEK_scbulk_BarcodeAnnotations.csv")

mat.ss3.bulk <- loom2matrix("data/SS3_RNA_bulk_Readcount_Downsampled100k.loom")
mat2.ss3.bulk <- loom2matrix("data/SS3_RNA_bulk_UMIcount_Downsampled100k.loom")
mat.ss3.hek <- loom2matrix("data/SS3_HEK_sc_Readcount_Downsampled100k.loom")
mat2.ss3.hek <- loom2matrix("data/SS3_HEK_sc_UMIcount_Downsampled100k.loom")

idx.col.ss3.bulk <- intersect(meta.ss3[sample_source == "bulk", sample_barcode], colnames(mat.ss3.bulk))
idx.col.ss3.hek <- intersect(meta.ss3[sample_source == "HEK", sample_barcode], colnames(mat.ss3.hek))

sce.ss3.bulk <- SingleCellExperiment(assays = list(counts = mat.ss3.bulk[, idx.col.ss3.bulk], umi = mat2.ss3.bulk[, idx.col.ss3.bulk]), colData=meta.ss3[match(idx.col.ss3.bulk, sample_barcode)])
sce.ss3.hek <- SingleCellExperiment(assays = list(counts = mat.ss3.hek[, idx.col.ss3.hek], umi = mat2.ss3.hek[, idx.col.ss3.hek]), colData=meta.ss3[match(idx.col.ss3.hek, sample_barcode)])

ggplot(data.table(y=colSums(assay(sce.ss3.bulk)>0), x= factor(sce.ss3.bulk$inhibitor_concentration, levels = c(0,0.006, 0.015, 0.03,0.06,0.15,0.3,0.6,1.5,2.25,3,6,9,"RRI")), sample_source = sce.ss3.bulk$sample_source), aes(x=x,y=y)) + geom_boxplot() + facet_wrap(~sample_source) + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS3") + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(data.table(y=colSums(assay(sce.ss3.hek)>0), x= factor(sce.ss3.hek$inhibitor_concentration, levels = c(0.15,0.3,0.6,"RRI")), sample_source = sce.ss3.hek$sample_source), aes(x=x,y=y)) + geom_boxplot() + facet_wrap(~sample_source) + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS3") + theme_cowplot()  + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggplot(data.table(y=colSums(assay(sce.ss3.bulk,"umi")), x= factor(sce.ss3.bulk$inhibitor_concentration, levels = c(0,0.006, 0.015, 0.03,0.06,0.15,0.3,0.6,1.5,2.25,3,6,9,"RRI")), sample_source = sce.ss3.bulk$sample_source), aes(x=x,y=y)) + geom_boxplot() + facet_wrap(~sample_source) + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS3") + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(data.table(y=colSums(assay(sce.ss3.hek,"umi")), x= factor(sce.ss3.hek$inhibitor_concentration, levels = c(0.15,0.3,0.6,"RRI")), sample_source = sce.ss3.hek$sample_source), aes(x=x,y=y)) + geom_boxplot() + facet_wrap(~sample_source) + expand_limits(y=0) + labs(y="Genes", x=NULL, title="SS3") + theme_cowplot()  + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

## Commercial kits

meta.kit.bulk <- fread("data/Kit_RNA_bulk_BarcodeAnnotations.csv")

mat.kit.bulk <- loom2matrix("data/Kit_RNA_bulk.readcount.exon.downsampled_100000.loom")

sce.kit.bulk <- SingleCellExperiment(assays = list(counts = mat.kit.bulk), colData=meta.kit.bulk[match(colnames(mat.kit.bulk), sample_barcode)])

ggplot(data.table(y=colSums(assay(sce.kit.bulk)>0), x= factor(sce.kit.bulk$inhibitor_concentration, levels=c("RRI", 0, 0.12,0.24,0.3,0.6,1.2,3,12,150,1500,15000)), kit = sce.kit.bulk$kit), aes(x=x,y=y)) + geom_boxplot() + expand_limits(y=0) + facet_wrap(~kit, scales="free_x") + labs(y="Genes", x=NULL, title="Kit bulk") + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())

## SS3xpress single-cell

loom2matrix_new <- function(x){
  require(loomR)
  lfile <- connect(x)
  out <- t(lfile$matrix[,])
  colnames(out) <- lfile$col.attrs$cell_names[]
  row.names(out) <- lfile$row.attrs$gene_names[]
  lfile$close_all()
  return(out)
}

meta.ss3xpress <- fread("data/SS3xpress_HEK_sc_BarcodeAnnotations.csv")

mat.ss3xpress.all <- loom2matrix_new("data/SS3xpress_HEK_sc_Readcount_ALL.loom")
mat.ss3xpress.ds <- loom2matrix_new("data/SS3xpress_HEK_sc_Readcount_Downsampled100k.loom")
mat2.ss3xpress.all <- loom2matrix_new("data/SS3xpress_HEK_sc_UMIcount_ALL.loom")
mat2.ss3xpress.ds <- loom2matrix_new("data/SS3xpress_HEK_sc_UMIcount_Downsampled100k.loom")

sce.ss3xpress.all <- SingleCellExperiment(assays = list(counts = mat.ss3xpress.all, umi = mat2.ss3xpress.all), colData=meta.ss3xpress[match(colnames(mat.ss3xpress.all), sample_barcode)])
sce.ss3xpress.ds <- SingleCellExperiment(assays = list(counts = mat.ss3xpress.ds, umi = mat2.ss3xpress.ds), colData=meta.ss3xpress[match(colnames(mat.ss3xpress.ds), sample_barcode)])

stats.ss3xpress.ds <- fread("SS3xpress/downsampled_stats.txt")

ggplot(data.table(y=colSums(assay(sce.ss3xpress.all)>0), x= factor(sce.ss3xpress.all$inhibitor_concentration, levels=c(0, 0.06, 0.15, 0.3, 0.6, 1.5, 3, "RRI")), day=sce.ss3xpress.all$day, temp=sce.ss3xpress.all$temperature)[day == "day0"], aes(x=x,y=y)) + geom_boxplot() + expand_limits(y=0) + labs(y="Genes", x=NULL) + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())
ggplot(data.table(y=colSums(assay(sce.ss3xpress.all, "umi")), x= factor(sce.ss3xpress.all$inhibitor_concentration, levels=c(0, 0.06, 0.15, 0.3, 0.6, 1.5, 3, "RRI")), day=sce.ss3xpress.all$day, temp=sce.ss3xpress.all$temperature)[day == "day0"], aes(x=x,y=y)) + geom_boxplot() + expand_limits(y=0) + labs(y="Genes", x=NULL) + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())

ggplot(data.table(y=colSums(assay(sce.ss3xpress.ds)>0), x= factor(sce.ss3xpress.ds$inhibitor_concentration, levels=c(0, 0.06, 0.15, 0.3, 0.6, 1.5, 3, "RRI")), day=sce.ss3xpress.ds$day, temp=sce.ss3xpress.ds$temperature)[(temp == "RT" & day != "day14") | day == "day0"], aes(x=x,y=y)) + geom_boxplot() + facet_grid(~day) + labs(y="Genes", x=NULL) + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())
ggplot(data.table(y=colSums(assay(sce.ss3xpress.ds)>0), x= factor(sce.ss3xpress.ds$inhibitor_concentration, levels=c(0, 0.06, 0.15, 0.3, 0.6, 1.5, 3, "RRI")), day=sce.ss3xpress.ds$day, temp=sce.ss3xpress.ds$temperature)[temp == "4C" | day == "day0"], aes(x=x,y=y)) + geom_boxplot() + facet_grid(~day) + labs(y="Genes", x=NULL) + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())

ggplot(stats.ss3xpress.ds[(grepl("RT", condition) | grepl("day0", condition)) & d == 1e5 & !grepl("day14", day)], aes(x=factor(gsub("ug.*","",condition3),levels=c("no RRI", "0.6", "1.5", "3.0", "6.0", "RRI")), y=nUMIs_2/1e3)) + geom_boxplot() + facet_grid(~day) + coord_cartesian(ylim=c(30,60)) + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())
ggplot(stats.ss3xpress.ds[(grepl("4C", condition) | grepl("day0", condition)) & d == 1e5], aes(x=factor(gsub("ug.*","",condition3),levels=c("no RRI", "0.6", "1.5", "3.0", "6.0", "RRI")), y=nUMIs_2/1e3)) + geom_boxplot() + facet_grid(~day) + coord_cartesian(ylim=c(30,60)) + theme_cowplot() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank())

##
