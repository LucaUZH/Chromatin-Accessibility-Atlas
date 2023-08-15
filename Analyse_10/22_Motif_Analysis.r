## first line of code is to clear R's memory
rm(list=ls())

if (FALSE){
    install.packages('ggseqlogo')
    BiocManager::install("motifmatchr")
    BiocManager::install("chromVAR")
}

#libraries
library(Signac)
library(Seurat)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(qs)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cicero)




#####################
#Read processed Seurat Object
###################
seurat_object = qread("complete/GSE184462_merged_1_10_CoAc.rds")

seurat_object
p1 <- DimPlot(seurat_object, label = TRUE, pt.size = 0.1) + NoLegend()
p1

library(BSgenome.Hsapiens.UCSC.hg38)
seurat_object[["peaks"]]
#Removign scaffold debris
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(seurat_object)) %in% main.chroms)
seurat_object[["peaks"]] <- subset(seurat_object[["peaks"]], features = rownames(seurat_object[["peaks"]])[keep.peaks])


###################
# extract gene annotations from EnsDb
###################
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

genome(seurat_object) <- genome(annotations)

#Get list of motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
seurat_object <- AddMotifs(
  object = seurat_object,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

Idents(seurat_object) <- seurat_object$cell.type
ident_lyst <- unique(Idents(seurat_object))
ident_lyst

#find markers
DefaultAssay(seurat_object) <- "peaks"

da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = "Endothelial Cell (General) 1",
  ident.2 = "Adipocyte",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

head(da_peaks)

enriched.motifs <- FindMotifs(
  object = seurat_object,
  features = top.da.peak,
)

DefaultAssay(seurat_object) <- "peaks"
mp <- MotifPlot(
  object = seurat_object,
  motifs = head(rownames(enriched.motifs))
)

mp


if (TRUE){
ggsave(mp, 
       filename = "exports/motif_plot.png",
       device = "png",
       height = 4, width = 6, units = "in")    
}


enriched.motifs

###################
#Computing motif activities
###################
seurat_object <- RunChromVAR(
  object = seurat_object,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(seurat_object) <- 'chromvar'



# look at the activity of MA0079.4
p2 <- FeaturePlot(
  object = seurat_object,
  features = "MA1513.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

if (TRUE){
ggsave(p1, 
       filename = "exports/cluster_plot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in") 
    
ggsave(p2, 
       filename = "exports/featur_MA1513.1_plot.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in") 
}

differential.activity <- FindMarkers(
  object = seurat_object,
  ident.1 = "Adipocyte",
  ident.2 = "Endothelial Cell (General) 1",
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

mp <- MotifPlot(
  object = seurat_object,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)
mp

differential.activity

if (TRUE){
ggsave(mp, 
       filename = "exports/motif_plot.pdf",
       device = "pdf",
       height = 6, width = 8, units = "in") 
    

}

#qsave(seurat_object, "complete/GSE184462_merged_1_10_Motif.rds")

#qsave(seurat_object, "complete/GSE184462_merged_1_10_Motif.rds")







enriched.motifs

dplyr::filter(enriched.motifs, motif.name == "CEBPG")
