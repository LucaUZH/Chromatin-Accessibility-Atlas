## first line of code is to clear R's memory
rm(list=ls())

# Install Cicero
if (FALSE){
    BiocManager::install("cicero")
    remotes::install_github('satijalab/seurat-wrappers')
    remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")  
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    
}

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
library(cicero)




#####################
#Read processed Seurat Object
###################
#seurat_object <- qread("complete/GSE184462_merged_1_10_processed.rds")
seurat_object <- qread("complete/GSE149683_placenta_processed.rds")

seurat_object

# convert to CellDataSet format and make the cicero object
seurat.cds <- as.cell_data_set(x = seurat_object)
cicero_object <- make_cicero_cds(seurat.cds, reduced_coordinates = reducedDims(seurat.cds)$UMAP)

cicero_object

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

head(conns)

########################
#Find cis-co-accessible networks (CCANs)
########################
ccans <- generate_ccans(conns)

head(ccans)

########################
#Add links to a Seurat object
########################
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seurat_object) <- links




Idents(seurat_object) <- seurat_object$seurat_clusters
CoveragePlot(seurat_object, region = "VWF")

CoveragePlot(seurat_object, region = "chr12-6100000-6150000")

Links(seurat_object)

#qsave(seurat_object, "complete/GSE184462_merged_1_10_CoAc.rds")

seurat_object <- qread("complete/GSE149683_placenta_CoAc.rds")



Idents(seurat_object) <- seurat_object$cell.type
CoveragePlot(seurat_object, region = "VWF", 
            idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Pericyte (General) 1", "Adipocyte"))

CoveragePlot(seurat_object, region = "chr12-6100000-6150000",
            idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Pericyte (General) 1", "Adipocyte"))






