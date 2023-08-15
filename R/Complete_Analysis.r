## first line of code is to clear R's memory
rm(list=ls())

#libraries
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
  BiocManager::install("EnsDb.Hsapiens.v86")
library(Signac)
library(Seurat)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(readr)
library(qs)
library(dplyr)
library(stringr)
library(Azimuth)
library(harmony)
library(cicero)
library(SeuratWrappers)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)


#####################
#Read Initial Seurat Object and filter it
###################
path <- paste("/", file, "_newpeaks.rds", sep = "")
seurat_object = qread(path)
 
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



path <- paste("exports/", file, "_QC.rds", sep = "")
#qsave(seurat_object, path)





path <- paste("exports/", file, "_QC.rds", sep = "")
#seurat_object <- qread(path)

##########################
#Normalization and linear dimensional reduction
#########################
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q50')
seurat_object <- RunSVD(seurat_object)

DepthCor(seurat_object)

##########################
#Non-linear dimension reduction and clustering
#########################
seurat_object <- RunHarmony(seurat_object, group.by.vars = "tissue", reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
seurat_object <- RunUMAP(seurat_object, reduction = "lsi", dims = 2:30)
seurat_object <- FindNeighbors(object = seurat_object, reduction = 'lsi', dims = 2:30)
seurat_object <- FindClusters(object = seurat_object, verbose = FALSE, algorithm = 3)
DimPlot(object = seurat_object, group.by = 'tissue') + NoLegend()



#Dimension plots
DimPlot(object = seurat_object, group.by = 'seurat_clusters', label = TRUE)
DimPlot(object = seurat_object, group.by = 'cell.type') + NoLegend()

#Dimension plot filtered for single cell type
for (i in seq(1,50)){

ct <- unique(seurat_object$cell.type)[i]
st <- subset( x = seurat_object, subset = cell.type == ct)
p1 <- DimPlot(object = st, group.by = "seurat_clusters", label = TRUE) + NoLegend() + ggtitle(ct)

    
print(p1)
    
}

##########################
#Create a gene activity matrix
#########################
#Calculate gene activity
gene.activities <- GeneActivity(seurat_object)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
seurat_object[['RNA']] <- CreateAssayObject(counts = gene.activities)

path <- paste("report/", file, "_nonnormalizedRNA.rds", sep = "")
#qsave(seurat_object, path)

#normalize gen activity data
seurat_object <- NormalizeData(
  object = seurat_object,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat_object$nCount_RNA)
)



#Feature plot
DefaultAssay(seurat_object) <- 'RNA'

fp <-FeaturePlot(
  object = seurat_object,
  features = c('KCTD8', 'PKHD1L1'),
  pt.size = 1.3,
  max.cutoff = 'q95',
  ncol = 2
)
fp
#Fibroblasts: VIM, PDGFRA
#Adipocyte: LEP, HOXC8
#Endothel: ANO2, VWF
#Mesothel: KCTD8, PKHD1L1

DefaultAssay(seurat_object) <- 'peaks'

Idents(seurat_object) <- seurat_object$cell.type
path <- paste("report/", file, "_processed.rds", sep = "")
#qsave(seurat_object, path)

##########################
#Find differentially accessible peaks
#########################
# change back to working with peaks instead of gene activities
DefaultAssay(seurat_object) <- 'peaks'
id <- unique(Idents(seurat_object))

#Calculate differentialy accessible peaks between two cell types
da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = id[1], #"Mesothelial Cell",
  ident.2 = id[-1], #"Endothelial Cell (General) 1",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

head(da_peaks)

#Violin plot
plot1 <- VlnPlot(
  object = seurat_object,
  features = rownames(da_peaks)[1],
  pt.size = 0.05,
  idents = c("Mesothelial Cell","Endothelial Cell (General) 1")
)

#feature plot
plot2 <- FeaturePlot(
  object = seurat_object,
  features = rownames(da_peaks)[2],
  pt.size = 1.3
)

plot1 | plot2

##########################
#Run Azimuth
#########################
DefaultAssay(seurat_object) <- "RNA"

az <- RunAzimuth(seurat_object,  reference = "adiposeref", assay = "RNA")

DimPlot(az, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
DimPlot(az, group.by = "predicted.celltype.l2")

az[[]][c(9,21,23)]

#####################
#Find co-accessible peaks with cicero
###################
#seurat_object = qread("complete/GSE184462_merged_1_10_processed.rds")

# convert to CellDataSet format and make the cicero object
DefaultAssay(seurat_object) <- "peaks"
seurat.cds <- as.cell_data_set(x = seurat_object)
cicero_object <- make_cicero_cds(seurat.cds, reduced_coordinates = reducedDims(seurat.cds)$UMAP)


########################
#Find Cicero connections
########################
# get the chromosome sizes from hg19
data(human.hg19.genome)
genome.df <- human.hg19.genome


# use chromosome 1 to save some time
# omit this step to run on the whole genome
#genome.df <- subset(human.hg19.genome, V1 == "chr18")

# run cicero
conns <- run_cicero(cicero_object, genomic_coords = genome.df, sample_num = 100)

########################
#Find cis-co-accessible networks (CCANs)
########################
ccans <- generate_ccans(conns)

########################
#Add links to a Seurat object
########################
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seurat_object) <- links

CoveragePlot(seurat_object, region = "VWF",
            idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Macrophage (General)"))

path <- paste("report/", file, "_CoAc.rds", sep = "")
#qsave(seurat_object, path)

#####################
#Motif analysis
###################
path <- paste("report/", file, "_CoAc.rds", sep = "")
#seurat_object = qread(path)

#Removing scaffold debris
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(seurat_object)) %in% main.chroms)
seurat_object[["peaks"]] <- subset(seurat_object[["peaks"]], features = rownames(seurat_object[["peaks"]])[keep.peaks])

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



ident_lyst <- unique(Idents(seurat_object))
ident_lyst

#find markers
DefaultAssay(seurat_object) <- "peaks"

da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = ident_lyst[1],
  ident.2 = ident_lyst[-1],
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

#Search for enriched motifs
enriched.motifs <- FindMotifs(
  object = seurat_object,
  features = top.da.peak,
)

DefaultAssay(seurat_object) <- "peaks"
MotifPlot(
  object = seurat_object,
  motifs = head(rownames(enriched.motifs))
)

###################
#Computing motif activities
###################
DefaultAssay(seurat_object) <- 'peaks'
seurat_object <- RunChromVAR(
  object = seurat_object,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(seurat_object) <- 'chromvar'

# look at the activity of Mef2c
#Plot nomral dimension plot with clusters
p1 <- DimPlot(seurat_object, label = TRUE, pt.size = 0.1) + NoLegend()

p2 <- FeaturePlot(
  object = seurat_object,
  features = "MA0079.4",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

differential.activity <- FindMarkers(
  object = seurat_object,
  ident.1 = 0,
  ident.2 = 1,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = seurat_object,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)

path <- paste("report/", file, "_Motif.rds", sep = "")
#qsave(seurat_object, path)

#####################
#TF Footprinting
###################
path <- paste("report/", file, "_Motif.rds", sep = "")
#seurat_object = qread(path)

Idents(seurat_object)  <- seurat_object$cell.type
unique(Idents(seurat_object))

# gather the footprinting information for sets of motifs
seurat_object <- Footprint(
  object = seurat_object,
  motif.name = c("EGR1", "Wt1", "KLF5"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p3 <- PlotFootprint(seurat_object, features = c("EGR1", "Wt1", "KLF5"), label = FALSE,
                    idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Macrophage (General)"))
p4 <- p3+ patchwork::plot_layout(ncol = 1)
p4

path <- paste("report/", file, "_TF_Footprinted.rds", sep = "")
#qsave(seurat_object, path)
