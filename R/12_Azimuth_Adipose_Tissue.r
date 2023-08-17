## first line of code is to clear R's memory
rm(list=ls())

#remotes::install_github('satijalab/azimuth', ref = 'master')
#remotes::install_github("satijalab/seurat-data", "seurat5", force = TRUE)
#remotes::install_github("satijalab/azimuth", "seurat5", force = TRUE)

#libraries
if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE))
  BiocManager::install("EnsDb.Hsapiens.v75")
library(Signac)
library(Seurat)
library(rtracklayer)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(readr)
library(qs)
library(dplyr)
library(stringr)
library(harmony)
library(Azimuth)
library(EnsDb.Hsapiens.v86)

adipose_tissue <- readRDS("complete/GSE184462_merged_1_3.rds")

peaks <- CallPeaks(
  object = adipose_tissue,
  group.by = "cell.type"
)

metadata <- adipose_tissue[[]]
frags <- Fragments(adipose_tissue)


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

adipose_tissue <- SeruatObj
saveRDS(SeruatObj, "complete/GSE184462_merged_1_3_newpeaks.rds")

#####################
#Annotation
###################
adipose_tissue <- readRDS("complete/GSE184462_merged_1_3_newpeaks.rds")

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(adipose_tissue) <- annotations
    

#####################
#QC Metrics
###################
# compute nucleosome signal score per cell
adipose_tissue <- NucleosomeSignal(object = adipose_tissue)

# compute TSS enrichment score per cell
adipose_tissue <- TSSEnrichment(object = adipose_tissue, fast = FALSE)

##########################
#Normalization and linear dimensional reduction
#########################
adipose_tissue <- RunTFIDF(adipose_tissue)
adipose_tissue <- FindTopFeatures(adipose_tissue, min.cutoff = 'q0')
adipose_tissue <- RunSVD(adipose_tissue)


##########################
#Non-linear dimension reduction and clustering
##########################
adipose_tissue <- RunHarmony(adipose_tissue, group.by.vars = "tissue", reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
adipose_tissue <- RunUMAP(adipose_tissue, reduction = "harmony", dims = 2:30)
adipose_tissue <- FindNeighbors(object = adipose_tissue, reduction = 'harmony', dims = 2:30)
adipose_tissue <- FindClusters(object = adipose_tissue, verbose = FALSE, algorithm = 3)
DimPlot(object = adipose_tissue, group.by = 'tissue') + NoLegend()

DimPlot(object = adipose_tissue, group.by = 'cell.type') + NoLegend()
DimPlot(object = adipose_tissue, group.by = 'cell.type')

##########################
#Create a gene activity matrix
#########################
gene.activities <- GeneActivity(adipose_tissue)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
adipose_tissue[['RNA']] <- CreateAssayObject(counts = gene.activities)

#adipose_tissue <- NormalizeData(
#  object = adipose_tissue,
#  assay = 'RNA',
#  normalization.method = 'LogNormalize',
#  scale.factor = median(adipose_tissue$nCount_RNA)
#)

#saveRDS(adipose_tissue, "complete/GSE184462_merged_1_3_processed.rds")

adipose_tissue <- readRDS("complete/GSE184462_merged_1_3_processed.rds")

DefaultAssay(adipose_tissue) <- "RNA"

az <- RunAzimuth(adipose_tissue, reference = "adiposeref", assay = "RNA")
DimPlot(az, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()

az[[]][c(9,20,22)]


az_plot <- DimPlot(az, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3)
az_plot
if (FALSE){
ggsave(az_plot, 
       filename = "exports/az_plot.pdf",
       device = "pdf",
       height = 6, width = 8, units = "in") 
    

}

#Azimuth with different reference
az <- RunAzimuth(adipose_tissue, reference = "pbmcref", assay = "RNA")
DimPlot(az, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
az[[]][c(9,18,20)]

#Mesothelium
subset(az, cell.type == unique(az$cell.type)[1] )[[]][c(9,18,20)]

#Adipocytes
subset(az, cell.type == unique(az$cell.type)[10] )[[]][c(9,18,20)]

#Endothelium
subset(az, cell.type == unique(az$cell.type)[c(5, 9, 12, 18, 20, 34)] )[[]][c(9,18,20)]

#NK T-cells
subset(az, cell.type == unique(az$cell.type)[41] )[[]][c(9,18,20)]

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
az$correctcall[az$cell.type=="Lymphatic Endothelial Cell" & az$predicted.celltype.l1=="Lymphatic Endothelial"] <- TRUE
#Macrophages
az$correctcall[az$cell.type=="Macrophage (General)" & az$predicted.celltype.l1=="Macrophage"] <- TRUE
az$correctcall[az$cell.type=="Macrophage (General,Alveolar)" & az$predicted.celltype.l1=="Macrophage"] <- TRUE
#Pericytes
az$correctcall[az$cell.type=="Pericyte (General) 1" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (General) 2" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (General) 3" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (General) 4" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 2" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 3" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Cardiac Pericyte 4" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
az$correctcall[az$cell.type=="Pericyte (Esophageal Muscularis)" & az$predicted.celltype.l1=="Pericyte"] <- TRUE
#Smooth Muscle
az$correctcall[az$cell.type=="Smooth Muscle (General)" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Vascular Smooth Muscle 1" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Vascular Smooth Muscle 2" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
az$correctcall[az$cell.type=="Smooth Muscle (GE Junction)" & az$predicted.celltype.l1=="Smooth Muscle"] <- TRUE
#Immune cells
az$correctcall[az$cell.type=="Naive T cell" & az$predicted.celltype.l1=="T"] <- TRUE
az$correctcall[az$cell.type=="Natural Killer T Cell" & az$predicted.celltype.l1 %in% c("T", "NK")] <- TRUE
az$correctcall[az$cell.type=="Mast Cell" & az$predicted.celltype.l1=="Mast"] <- TRUE



az_types <- c("Mesothelial Cell", "Adipocyte", "Endothelial Cell (General) 1", "Endothelial Cell (General) 2",
             "Endothelial Cell (General) 3", "Endothelial Cell (Myocardial)","Lymphatic Endothelial Cell", 
              "Pericyte (General) 1", "Pericyte (General) 2", "Pericyte (General) 3", "Pericyte (General) 4", 
              "Cardiac Pericyte 2", "Cardiac Pericyte 3", "Cardiac Pericyte 4", "Pericyte (Esophageal Muscularis)", 
              "Smooth Muscle (General)", "Vascular Smooth Muscle 1", "Vascular Smooth Muscle 2",
              "Naive T cell", "Natural Killer T Cell", "Mast Cell")
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

colnames(az_quality) <- c("total", "correct")
az_quality <- mutate(az_quality, false = total - correct)
az_quality <- mutate(az_quality, perc.correct = correct/total)
az_quality <- mutate(az_quality, perc.false = false/total)
az_quality$id <- rownames(az_quality)


az_quality$id <- c("Mesothelial Cell", "Adipocyte", "Endothelial Cell", "Endothelial Cell",
             "Endothelial Cell", "Endothelial Cell","Lymphatic Endothelial Cell", 
              "Pericyte", "Pericyte", "Pericyte", "Pericyte", 
              "Pericyte", "Pericyte", "Pericyte", "Pericyte", 
              "Smooth Muscle", "Smooth Muscle", "Smooth Muscle",
              "Naive T cell", "Natural Killer T Cell", "Mast Cell")



az_quality

df <- az_quality
df$id <- NULL

#Combine enodthelial cells
df[3,]<-df[3,]+df[4,]+df[5,]+df[6,]
df <-df[c(-4,-5,-6),]

#combine pericytes
df[5,]<-df[5,]+df[6,]+df[7,]+df[8,]+df[9,]+df[10,]+df[11,]+df[12,]
df <-df[c(-6,-7,-8,-9,-10,-11,-12),]
df[6,]<-df[6,]+df[7,]+df[8,]
df <-df[c(-7,-8),]

df$id <- c("Mesothelial Cell", "Adipocyte", "Endothelial Cell","Lymphatic Endothelial Cell", 
              "Pericyte","Smooth Muscle", "Naive T cell", "Natural Killer T Cell", "Mast Cell")


ggplot(df, aes(x= id, y = correct/total, color = id, fill = id))+
geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
ggtitle("Correct predictions")

ggplot(df, aes(x= id, y = correct/total, color = id, fill = id))+
geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
ggtitle("Correct predictions")+
geom_text(aes(label=total), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)









