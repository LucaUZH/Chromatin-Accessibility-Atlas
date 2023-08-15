## first line of code is to clear R's memory
rm(list=ls())

#libraries
library(Signac)
library(Seurat)
library(rtracklayer)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(qs)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)




#####################
#Read processed Seurat Object
###################
seurat_object = qread("complete/GSE184462_merged_1_10_Motif.rds")

Idents(seurat_object)  <- seurat_object$cell.type

unique(Idents(seurat_object))

#find markers
DefaultAssay(seurat_object) <- "peaks"

da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = "Mesothelial Cell",
  ident.2 = "Endothelial Cell (General) 1", 
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

enriched.motifs <- FindMotifs(
  object = seurat_object,
  features = top.da.peak,
)

MotifPlot(
  object = seurat_object,
  motifs = head(rownames(enriched.motifs))
)

enriched.motifs

# gather the footprinting information for sets of motifs
DefaultAssay(seurat_object) <- "peaks"
seurat_object <- Footprint(
  object = seurat_object,
  motif.name = c("CEBPG"), #, "Wt1", "KLF5"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)



# plot the footprint data for each group of cells
p3 <- PlotFootprint(seurat_object, features = c("CEBPG"), label = FALSE,
                    idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Macrophage (General)", "Adipocyte"))
p4 <- p3+ patchwork::plot_layout(ncol = 1)
p4



if (TRUE){
ggsave(p4, 
       filename = "exports/footprint_plot.pdf",
       device = "pdf",
       height = 6, width = 15, units = "in") 

}


