## first line of code is to clear R's memory
rm(list=ls())

#libraries
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
  BiocManager::install("EnsDb.Hsapiens.v86")
library(Signac)
library(Seurat)
library(rtracklayer)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(readr)
library(qs)
library(dplyr)
library(stringr)
library(Azimuth)
library(harmony)
library(EnsDb.Hsapiens.v86)

#####################
#Read Initial Seurat Object and filter it
###################
seurat_object = readRDS("complete/GSE184462_merged_1_10_test.rds")

peaks <- CallPeaks(
  object = seurat_object,
  group.by = "cell.type"
)

metadata <- seurat_object[[]]
frags <- Fragments(seurat_object)


peak_matrix <- FeatureMatrix(
  fragments = frags,
  features = peaks
)

chrom_assay <- CreateChromatinAssay(
  counts = peak_matrix,
  sep = c(":", "-"),
  fragments = frags
)

SeruatObj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

SeruatObj

#seurat_object <- SeruatObj
#saveRDS(SeruatObj, "complete/GSE184462_merged_1_10_newpeaks.rds")

seurat_object <- readRDS("complete/GSE184462_merged_1_10_newpeaks.rds")
 
#Add barcode to metadata
seurat_object$barcode <- Cells(seurat_object)

#####################
#Annotation
###################
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(seurat_object) <- annotations
    


#####################
#QC Metrics
###################
# compute nucleosome signal score per cell
seurat_object <- NucleosomeSignal(object = seurat_object)

# compute TSS enrichment score per cell
seurat_object <- TSSEnrichment(object = seurat_object, fast = FALSE)

DensityScatter(seurat_object, x = 'nCount_peaks', y = 'tsse', log_x = TRUE, quantiles = TRUE)

VlnPlot(
  object = seurat_object,
  features = c('nCount_peaks', 'tsse', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

DensityScatter(seurat_object, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)


#filtering out bad reads
seurat_object <- subset(
    x = seurat_object,
    subset = nCount_peaks < 4277 &
    nCount_peaks > 302 &
    TSS.enrichment < 6.88
)

DensityScatter(seurat_object, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)



#qsave(seurat_object, "exports/GSE184462_merged_1_10_QC.rds")



##########################
#Normalization and linear dimensional reduction
#########################
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q50')
seurat_object <- RunSVD(seurat_object)

DepthCor(seurat_object)





##########################
#Comparison to non-batchcontroled
#########################
if (FALSE){
unintegrated <- RunUMAP(seurat_object, dims = 2:50, reduction = 'lsi')
unintegrated <- FindNeighbors(object = unintegrated, reduction = 'lsi', dims = 2:30, k.param = 50)
unintegrated <- FindClusters(object = unintegrated, verbose = FALSE, algorithm = 3)
DimPlot(object = unintegrated, group.by = 'tissue') + NoLegend()
    
}



##########################
#Batch Correction (harmony)
#########################
if (FALSE){
seurat_object <- RunHarmony(seurat_object, group.by.vars = "tissue", reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30)
DimPlot(seurat_object, group.by = 'tissue', pt.size = 0.1)
}


##########################
#Non-linear dimension reduction and clustering
#########################
seurat_object <- RunHarmony(seurat_object, group.by.vars = "tissue", reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 2:30)
seurat_object <- FindNeighbors(object = seurat_object, reduction = 'harmony', dims = 2:30)
seurat_object <- FindClusters(object = seurat_object, verbose = FALSE, algorithm = 3)
DimPlot(object = seurat_object, group.by = 'tissue') + NoLegend()




length(unique(seurat_object$seurat_clusters))
length(unique(seurat_object$cell.type))

DimPlot(object = seurat_object, label = TRUE)
DimPlot(object = seurat_object, group.by = 'cell.type') + NoLegend()

unique(seurat_object$cell.type)

ct <- unique(seurat_object$cell.type)[5]
st <- subset( x = seurat_object, subset = cell.type == ct)

DimPlot(object = st, group.by = "seurat_clusters", label = TRUE) + NoLegend() + ggtitle(ct)





##########################
#Create a gene activity matrix
#########################
gene.activities <- GeneActivity(seurat_object)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
seurat_object[['RNA']] <- CreateAssayObject(counts = gene.activities)
seurat_object <- NormalizeData(
  object = seurat_object,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat_object$nCount_RNA)
)



DefaultAssay(seurat_object) <- 'RNA'

FeaturePlot(
  object = seurat_object,
  features = c('ACACB', 'ADIPOQ', 'ANO2', 'VWF'),
  pt.size = 0.05,
  max.cutoff = 'q95',
  ncol = 2
)
#Fibroblasts: VIM, PDGFRA
#Adipocyte: LEP, HOXC8
#Endothel: LDB2, JAM2, ANO2, VWF
#Adipocytes: GPAM, ADIPOQ, ACACB, PCDH9, GHR, TRHDE, SORBS1, PLIN1, PDE3B, WDPCP

DefaultAssay(seurat_object) <- 'peaks'

length(Features(seurat_object[["RNA"]]))
#qsave(seurat_object, "complete/GSE184462_merged_1_10_processed.rds")

unique(seurat_object$cell.type)

# change back to working with peaks instead of gene activities
DefaultAssay(seurat_object) <- 'peaks'
Idents(seurat_object) <- seurat_object$cell.type

da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = "Adipocyte",
  ident.2 = "Endothelial Cell (General) 1",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = seurat_object,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("Mesothelial Cell","Endothelial Cell (General) 1")
)
plot2 <- FeaturePlot(
  object = seurat_object,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

da_peaks

#seurat_object <- qread("complete/GSE184462_merged_1_10_processed.rds")

seurat_object_unnoramlized<- seurat_object
seurat_object_unnoramlized[['RNA']] <- NULL
seurat_object_unnoramlized[['RNA']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(seurat_object_unnoramlized) <- "RNA"

az <- RunAzimuth(seurat_object_unnoramlized,  reference = "adiposeref", assay = "RNA")

DimPlot(az, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
DimPlot(az, group.by = "predicted.celltype.l2")

az[[]][c(9,21,23)]

colnames(az[[]])

#Check quality
unique(az$cell.type)
unique(az$predicted.celltype.l1)
unique(az$predicted.celltype.l2)

az$correctcall <- FALSE
#Mesothel
az$correctcall[az$cell.type=="Mesothelial Cell" & az$predicted.celltype.l1=="Mesothelial"] <- TRUE
#Adipocytes
az$correctcall[az$cell.type=="Adipocyte" & az$predicted.celltype.l1 %in% c("Adipocyte", "ASPC")] <- TRUE
#Endothel
az$correctcall[az$cell.type=="Endothelial Cell (General) 1" & az$predicted.celltype.l1=="Endothelial"] <- TRUE
az$correctcall[az$cell.type=="Endothelial Cell (General) 2" & az$predicted.celltype.l1=="Endothelial"] <- TRUE
az$correctcall[az$cell.type=="Endothelial Cell (General) 3" & az$predicted.celltype.l1=="Endothelial"] <- TRUE
az$correctcall[az$cell.type=="Endothelial Cell (Myocardial)" & az$predicted.celltype.l1=="Endothelial"] <- TRUE
az$correctcall[az$cell.type=="Alveolar Capillary Endothelial Cell" & az$predicted.celltype.l1=="Endothelial"] <- TRUE
az$correctcall[az$cell.type=="Endothelial (Exocrine Tissues)" & az$predicted.celltype.l1=="Endothelial"] <- TRUE
az$correctcall[az$cell.type=="Lymphatic Endothelial Cell" & az$predicted.celltype.l1=="Lymphatic Endothelial"] <- TRUE
#Macrophages
az$correctcall[az$cell.type=="Macrophage (General)" & az$predicted.celltype.l1=="Macrophage"] <- TRUE
az$correctcall[az$cell.type=="Macrophage (General,Alveolar)" & az$predicted.celltype.l1=="Macrophage"] <- TRUE
#Pericytes
az$correctcall[az$cell.type=="Pericyte (General) 1" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (General) 2" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (General) 3" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (General) 4" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 1" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 2" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 3" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 4" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (Esophageal Muscularis)" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
#Smooth Muscle
az$correctcall[az$cell.type=="Smooth Muscle (General)" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Vascular Smooth Muscle 1" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Vascular Smooth Muscle 2" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Smooth Muscle (Esophageal Muscularis) 3" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Smooth Muscle (GE Junction)" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
#Immune cells
az$correctcall[az$cell.type=="Memory B Cell" & az$predicted.celltype.l1=="B"] <- TRUE
az$correctcall[az$cell.type=="Naive T cell" & az$predicted.celltype.l1=="T"] <- TRUE
az$correctcall[az$cell.type=="Natural Killer T Cell" & az$predicted.celltype.l1 %in% c("T", "NK")] <- TRUE
az$correctcall[az$cell.type=="Mast Cell" & az$predicted.celltype.l1=="Mast"] <- TRUE



az_types <- c("Mesothelial Cell", "Adipocyte", "Endothelial Cell (General) 1", "Endothelial Cell (General) 2",
             "Endothelial Cell (General) 3", "Endothelial Cell (Myocardial)", "Endothelial (Exocrine Tissues)",
             "Lymphatic Endothelial Cell", "Pericyte (General) 1", "Pericyte (General) 2", "Pericyte (General) 3", 
             "Pericyte (General) 4", "Cardiac Pericyte 1", "Cardiac Pericyte 2", "Cardiac Pericyte 3", "Cardiac Pericyte 4",
             "Pericyte (Esophageal Muscularis)", "Smooth Muscle (General)", "Vascular Smooth Muscle 1", "Vascular Smooth Muscle 2",
             "Smooth Muscle (Esophageal Muscularis) 3")
az_sub <- subset(az, cell.type %in% az_types)

az_sub[[]]

az_quality <- data.frame(row.names = az_types)


for (i in seq(1,length(az_types))){
    temp <- subset(az_sub, cell.type == az_types[i])
    total <- length(temp$cell.type)
    correct <-  sum(temp$correctcall, na.rm = TRUE)
    az_quality[i,1] <- total
    az_quality[i,2] <- correct
    
}


az_quality <- mutate(az_quality, V3 = V2/V1)
az_quality
