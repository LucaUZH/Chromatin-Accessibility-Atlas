## first line of code is to clear R's memory
rm(list=ls())

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
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)

#####################################
#Cluster and cell types
#####################################
seurat_object <- qread("report/GSE184462_merged_1_20_Motif.rds")

id <- seurat_object$cell.type
id <- as.data.frame(id)

#Assing summariesed id's
id[id =="Adipocyte"] <- "Adipocyte"
id[id =="Alveolar Capillary Endothelial Cell"] <- "Endothelial"
id[id =="Alveolar Type 2 (AT2) Cell"] <- "Alveolar"
id[id =="Alverolar Type 2,Immune"] <- "Alveolar"
id[id =="Cardiac Fibroblasts"] <- "Fibroblast"
id[id =="Cardiac Pericyte 1"] <- "Pericyte"
id[id =="Cardiac Pericyte 2"] <- "Pericyte"
id[id =="Cardiac Pericyte 3"] <- "Pericyte"
id[id =="Cardiac Pericyte 4"] <- "Pericyte"
id[id =="Chief Cell"] <- "Chief Cell"
id[id =="CNS,Enteric Neuron"] <- "Neuron"
id[id =="Cortical Epithelial-like"] <- "Epithelial"
id[id =="Ductal Cell (Pancreatic)"] <- "Ductal Cell"
id[id =="Endocardial Cell"] <- "Endocardial"
id[id =="Endothelial (Exocrine Tissues)"] <- "Endothelial"
id[id =="Endothelial Cell (General) 1"] <- "Endothelial"
id[id =="Endothelial Cell (General) 2"] <- "Endothelial"
id[id =="Endothelial Cell (General) 3"] <- "Endothelial"
id[id =="Endothelial Cell (Myocardial)"] <- "Endothelial"
id[id =="Fibroblast (Epithelial)"] <- "Fibroblast"
id[id =="Fibroblast (General)"] <- "Fibroblast"
id[id =="Fibroblast (Liver Adrenal)"] <- "Fibroblast"
id[id =="Fibroblast (Peripheral Nerve)"] <- "Fibroblast"
id[id =="Foveolar Cell"] <- "Foveolar Cell"
id[id =="Keratinocyte 1"] <- "Keratinocyte"
id[id =="Luteal Cell (Ovarian)"] <- "Luteal Cell"
id[id =="Lymphatic Endothelial Cell"] <- "Endothelial"
id[id =="Macrophage (General,Alveolar)"] <- "Macrophage"
id[id =="Macrophage (General)"] <- "Macrophage"
id[id =="Mast Cell"] <- "Mast Cell"
id[id =="Memory B Cell"] <- "B"
id[id =="Mesothelial Cell"] <- "Mesothelial"
id[id =="Naive T cell"] <- "T"
id[id =="Natural Killer T Cell"] <- "T"
id[id =="Pancreatic Acinar Cell"] <- "Acinar Cell"
id[id =="Pericyte (Esophageal Muscularis)"] <- "Pericyte"
id[id =="Pericyte (General) 1"] <- "Pericyte"
id[id =="Pericyte (General) 2"] <- "Pericyte"
id[id =="Pericyte (General) 3"] <- "Pericyte"
id[id =="Pericyte (General) 4"] <- "Pericyte"
id[id =="Peripheral Nerve Stromal"] <- "Neuron"
id[id =="Plasma Cell"] <- "Plasma Cell"
id[id =="Schwann Cell (General)"] <- "Schwann Cell"
id[id =="Small Intestinal Enterocyte"] <- "Enterocyte"
id[id =="Smooth Muscle (Esophageal Muscularis) 3"] <- "Smooth Muscle"
id[id =="Smooth Muscle (GE Junction)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (General)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Vaginal)"] <- "Smooth Muscle"
id[id =="T Lymphocyte 1 (CD8+)"] <- "T"
id[id =="T lymphocyte 2 (CD4+)"] <- "T"
id[id =="Transitional Zone Cortical Cell"] <- "Cortical Cell"
id[id =="Vascular Smooth Muscle 1"] <- "Smooth Muscle"
id[id =="Vascular Smooth Muscle 2"] <- "Smooth Muscle"
id[id =="Zona Fasciculata Cortical Cell"] <- "Cortical Cell"
id[id =="Zona Glomerulosa Cortical Cell"] <- "Cortical Cell"

#Additional Cellt types (11-20)
id[id =="Airway Goblet Cell"] <- "Goblet Cell"
id[id =="Astrocyte 2"] <- "Astrocyte"
id[id =="Basal Epidermal (Skin)"] <- "Epidermal"
id[id =="Basal Epithelial (Mammary)"] <- "Epithelial"
id[id =="Blood Brain Barrier Endothelial Cell"] <- "Endothelial"
id[id =="Club Cell"] <- "Club Cell"
id[id =="Colon Epithelial Cell 1"] <- "Epithelial"
id[id =="Colon Epithelial Cell 2"] <- "Epithelial"
id[id =="Colon Epithelial Cell 3"] <- "Epithelial"
id[id =="Colonic Goblet Cell"] <- "Goblet Cell"
id[id =="Enterochromaffin Cell"] <- "Enterochromaffin Cell"
id[id =="Esophageal Epithelial Cell"] <- "Epithelial"
id[id =="Fibroblast (Gastrointestinal)"] <- "Fibroblast"
id[id =="Fibroblast (Sk Muscle Associated)"] <- "Fibroblast"
id[id =="GABAergic Neuron 1"] <- "Neuron"
id[id =="Gastric Neuroendocrine Cell"] <- "Neuron"
id[id =="Granular Epidermal (Skin)"] <- "Epidermal"
id[id =="Mammary Luminal Epithelial Cell 2"] <- "Epithelial"
id[id =="Melanocyte"] <- "Melanocyte"
id[id =="Myoepithelial (Skin)"] <- "Myoepithelial"
id[id =="Oligodendrocyte"] <- "Oligodendrocyte"
id[id =="Oligodendrocyte Precursor"] <- "Oligodendrocyte"
id[id =="Paneth Cell"] <- "Paneth Cell"
id[id =="Small Intestinal Goblet Cell"] <- "Adipocyte"
id[id =="Smooth Muscle (Colon) 1"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Colon) 2"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Mucosal)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Muscularis) 1"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Muscularis) 2"] <- "Smooth Muscle"
id[id =="Smooth Muscle (General Gastrointestinal)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Uterine)"] <- "Smooth Muscle"
id[id =="Tuft Cell"] <- "Tuft Cell"


seurat_object$id <- id

#Create data frame for each summarised id and how many cells they contain
id <- unique(seurat_object$id)
df <- data.frame(row.names = id)
meta <- as.data.frame(seurat_object[[]])
for (i in id){
    temp <- filter(meta, id == i)
    df[i,1] <- length(temp$id)
    
}

colnames(df) <- c("cell_count")

#Filter for cell types (id) with more than 100 cells
filter_id = rownames(filter(df, cell_count > 100 ))
temp <- subset(seurat_object, id %in% filter_id)

#Plot dim plot with summarised and filtered cell types
p1 <- DimPlot(object = seurat_object, label = TRUE)
p2 <- DimPlot(object = seurat_object, group.by = 'id', label = TRUE) +ggtitle("Dimension Plot")

p3 <- DimPlot(object = temp, label = TRUE)
p4 <- DimPlot(object = temp, group.by = 'id', label = TRUE)+ ggtitle("Dimension Plot, cell types")
p1
p2
p3
p4

ggsave(p1, 
       filename = "report/dimplot_summarised_legend.pdf",
       device = "pdf",
       height = 6, width = 20, units = "in") 
ggsave(p2, 
       filename = "report/dimplot_summarised_label.pdf",
       device = "pdf",
       height = 6, width = 13, units = "in") 
ggsave(p3, 
       filename = "report/dimplot_summarised_legend_red.pdf",
       device = "pdf",
       height = 6, width = 20, units = "in") 
ggsave(p4, 
       filename = "report/dimplot_summarised_label_red.pdf",
       device = "pdf",
       height = 6, width = 13, units = "in") 

for (i in seq(1,19)){

ct <- unique(temp$id)[i]
st <- subset( x = temp, subset = id == ct)
ct
p1 <- DimPlot(object = st) + NoLegend() + ggtitle(ct)

    
print(p1)
    
}

#Make dataframe with clusters as rows and cell types as columnes
cols = sort(unique(temp$id))
rows = sort(unique(temp$seurat_clusters))
meta <- as.data.frame(temp[[]])
df <- data.frame(row.names = rows)

#Count how many cells per cell types there are in each cluster
for (r in rows){
    for (l in seq(1,length(cols))){
        m <- dplyr::filter(meta, id == cols[l])
        m <- dplyr::filter(m, seurat_clusters == r)
        df[as.numeric(r)+1, l] <- length(m$id)        
    }  
}

colnames(df) <- cols

#Add the toal number of cells per cluster
df$total <- rowSums(df)



#Calculate purrity of each cluster
df$purrity <- 0
for (r in rows){
    df[r,20] <- max(df[-19][r,])/df[r,19]
}

#Add predominant cell typ of cluster to df
df$cluster_celltype <- colnames(df[-19])[max.col(df[-19],ties.method="first")]
#Add cluster as columne
df$cluster <- rows



df[c(19,20,21)]

ggplot(df, aes(x= cluster, y = purrity, color = cluster_celltype, fill = cluster_celltype))+
geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
ggtitle("Cluster purities")+
geom_text(aes(label=cluster_celltype), position = position_stack(vjust = 0.5), angle = 90, color="black", size=3.5)


#####################################
#DA peaks
#####################################
seurat_object <- qread("report/GSE184462_merged_1_20_processed.rds")

# change back to working with peaks instead of gene activities
DefaultAssay(seurat_object) <- 'peaks'
Idents(seurat_object) <- seurat_object$id
ident <- unique(Idents(seurat_object))

da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = "Endothelial", #"Mesothelial Cell",
  ident.2 = "Adipocyte", #rest,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
head(da_peaks)

plot1 <- VlnPlot(
  object = seurat_object,
  features = "chr12-6124336-6125604",
  pt.size = 0.05,
  idents = c("Endothelial", "Adipocyte")
)
plot2 <- FeaturePlot(
  object = seurat_object,
  features = "chr12-6124336-6125604",
  pt.size = 1.3
)

plot1 | plot2

# change back to working with peaks instead of gene activities
DefaultAssay(seurat_object) <- 'RNA'
Idents(seurat_object) <- seurat_object$id
ident <- unique(Idents(seurat_object))

da_peaks <- FindMarkers(
  object = seurat_object,
  ident.1 = "Endothelial", #"Mesothelial Cell",
  ident.2 = "Adipocyte", #rest,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
head(da_peaks)

#Filter for peaks more active in ident1, statistical significance and activity over 50%
da_peaks <- filter(da_peaks, pct.1 > pct.2)
da_peaks <- filter(da_peaks, pct.1 > 0.5)
da_peaks <- filter(da_peaks, p_val_adj < 0.05)
da_peaks

plot1 <- VlnPlot(
  object = seurat_object,
  features = rownames(da_peaks)[1],
  pt.size = 0,
  idents = c("Endothelial", "Adipocyte")
)
plot2 <- FeaturePlot(
  object = seurat_object,
  features = rownames(da_peaks)[1],
  pt.size = 1.3
)

plot1 | plot2

ggsave(plot1, 
       filename = "report/Violin_plot_Endo_Adipo.pdf",
       device = "pdf",
       height = 10, width = 10, units = "in") 

ggsave(plot2, 
       filename = "report/Feature_plot_Endo_Adipo.pdf",
       device = "pdf",
       height = 10, width = 10, units = "in") 

#####################################
#CRE
#####################################

links <- Links(seurat_object)
keep <- !is.na(mcols(links)$group)
links <- links[keep,]
links <- plyranges::filter(links, score > 0.22)
seurat_test <- seurat_object
Links(seurat_test) <- links


Idents(seurat_test) <- seurat_test$seurat_clusters
CoveragePlot(seurat_test, region = "chr1-1100000-1130000") + ggtitle("Coverage Plot C1orf159")
CoveragePlot(seurat_test, region = "chr12-6100000-6150000") + ggtitle("Coverage Plot VWF")

Idents(seurat_test) <- seurat_test$id
CoveragePlot(seurat_test, region = "chr1-1100000-1130000",
            idents = c("Mesothelial", "Endothelial", "Pericyte", "Adipocyte", "T")) + ggtitle("Coverage Plot C1orf159")
CoveragePlot(seurat_test, region = "chr12-6100000-6150000",
            idents = c("Mesothelial", "Endothelial", "Pericyte", "Adipocyte", "T")) + ggtitle("Coverage Plot VWF")

Idents(seurat_test) <- seurat_test$cell.type
CoveragePlot(seurat_test, region = "chr1-1100000-1130000",
            idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Pericyte (General) 1", "Adipocyte", "T Lymphocyte 1 (CD8+)"))+ ggtitle("Coverage Plot C1orf159")
CoveragePlot(seurat_test, region = "chr12-6100000-6150000",
            idents = c("Mesothelial Cell", "Endothelial Cell (General) 1", "Pericyte (General) 1", "Adipocyte", "T Lymphocyte 1 (CD8+)")) + ggtitle("Coverage Plot VWF")

##########################3
#motif & heatmap
#############################3
seurat_object <- qread("report/GSE184462_merged_1_20_Motif.rds")

id <- seurat_object$cell.type
id <- as.data.frame(id)
#Assing summariesed id's
id[id =="Adipocyte"] <- "Adipocyte"
id[id =="Alveolar Capillary Endothelial Cell"] <- "Endothelial"
id[id =="Alveolar Type 2 (AT2) Cell"] <- "Alveolar"
id[id =="Alverolar Type 2,Immune"] <- "Alveolar"
id[id =="Cardiac Fibroblasts"] <- "Fibroblast"
id[id =="Cardiac Pericyte 1"] <- "Pericyte"
id[id =="Cardiac Pericyte 2"] <- "Pericyte"
id[id =="Cardiac Pericyte 3"] <- "Pericyte"
id[id =="Cardiac Pericyte 4"] <- "Pericyte"
id[id =="Chief Cell"] <- "Chief Cell"
id[id =="CNS,Enteric Neuron"] <- "Neuron"
id[id =="Cortical Epithelial-like"] <- "Epithelial"
id[id =="Ductal Cell (Pancreatic)"] <- "Ductal Cell"
id[id =="Endocardial Cell"] <- "Endocardial"
id[id =="Endothelial (Exocrine Tissues)"] <- "Endothelial"
id[id =="Endothelial Cell (General) 1"] <- "Endothelial"
id[id =="Endothelial Cell (General) 2"] <- "Endothelial"
id[id =="Endothelial Cell (General) 3"] <- "Endothelial"
id[id =="Endothelial Cell (Myocardial)"] <- "Endothelial"
id[id =="Fibroblast (Epithelial)"] <- "Fibroblast"
id[id =="Fibroblast (General)"] <- "Fibroblast"
id[id =="Fibroblast (Liver Adrenal)"] <- "Fibroblast"
id[id =="Fibroblast (Peripheral Nerve)"] <- "Fibroblast"
id[id =="Foveolar Cell"] <- "Foveolar Cell"
id[id =="Keratinocyte 1"] <- "Keratinocyte"
id[id =="Luteal Cell (Ovarian)"] <- "Luteal Cell"
id[id =="Lymphatic Endothelial Cell"] <- "Endothelial"
id[id =="Macrophage (General,Alveolar)"] <- "Macrophage"
id[id =="Macrophage (General)"] <- "Macrophage"
id[id =="Mast Cell"] <- "Mast Cell"
id[id =="Memory B Cell"] <- "B"
id[id =="Mesothelial Cell"] <- "Mesothelial"
id[id =="Naive T cell"] <- "T"
id[id =="Natural Killer T Cell"] <- "T"
id[id =="Pancreatic Acinar Cell"] <- "Acinar Cell"
id[id =="Pericyte (Esophageal Muscularis)"] <- "Pericyte"
id[id =="Pericyte (General) 1"] <- "Pericyte"
id[id =="Pericyte (General) 2"] <- "Pericyte"
id[id =="Pericyte (General) 3"] <- "Pericyte"
id[id =="Pericyte (General) 4"] <- "Pericyte"
id[id =="Peripheral Nerve Stromal"] <- "Neuron"
id[id =="Plasma Cell"] <- "Plasma Cell"
id[id =="Schwann Cell (General)"] <- "Schwann Cell"
id[id =="Small Intestinal Enterocyte"] <- "Enterocyte"
id[id =="Smooth Muscle (Esophageal Muscularis) 3"] <- "Smooth Muscle"
id[id =="Smooth Muscle (GE Junction)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (General)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Vaginal)"] <- "Smooth Muscle"
id[id =="T Lymphocyte 1 (CD8+)"] <- "T"
id[id =="T lymphocyte 2 (CD4+)"] <- "T"
id[id =="Transitional Zone Cortical Cell"] <- "Cortical Cell"
id[id =="Vascular Smooth Muscle 1"] <- "Smooth Muscle"
id[id =="Vascular Smooth Muscle 2"] <- "Smooth Muscle"
id[id =="Zona Fasciculata Cortical Cell"] <- "Cortical Cell"
id[id =="Zona Glomerulosa Cortical Cell"] <- "Cortical Cell"

#Additional Cellt types (11-20)
id[id =="Airway Goblet Cell"] <- "Goblet Cell"
id[id =="Astrocyte 2"] <- "Astrocyte"
id[id =="Basal Epidermal (Skin)"] <- "Epidermal"
id[id =="Basal Epithelial (Mammary)"] <- "Epithelial"
id[id =="Blood Brain Barrier Endothelial Cell"] <- "Endothelial"
id[id =="Club Cell"] <- "Club Cell"
id[id =="Colon Epithelial Cell 1"] <- "Epithelial"
id[id =="Colon Epithelial Cell 2"] <- "Epithelial"
id[id =="Colon Epithelial Cell 3"] <- "Epithelial"
id[id =="Colonic Goblet Cell"] <- "Goblet Cell"
id[id =="Enterochromaffin Cell"] <- "Enterochromaffin Cell"
id[id =="Esophageal Epithelial Cell"] <- "Epithelial"
id[id =="Fibroblast (Gastrointestinal)"] <- "Fibroblast"
id[id =="Fibroblast (Sk Muscle Associated)"] <- "Fibroblast"
id[id =="GABAergic Neuron 1"] <- "Neuron"
id[id =="Gastric Neuroendocrine Cell"] <- "Neuron"
id[id =="Granular Epidermal (Skin)"] <- "Epidermal"
id[id =="Mammary Luminal Epithelial Cell 2"] <- "Epithelial"
id[id =="Melanocyte"] <- "Melanocyte"
id[id =="Myoepithelial (Skin)"] <- "Myoepithelial"
id[id =="Oligodendrocyte"] <- "Oligodendrocyte"
id[id =="Oligodendrocyte Precursor"] <- "Oligodendrocyte"
id[id =="Paneth Cell"] <- "Paneth Cell"
id[id =="Small Intestinal Goblet Cell"] <- "Adipocyte"
id[id =="Smooth Muscle (Colon) 1"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Colon) 2"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Mucosal)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Muscularis) 1"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Muscularis) 2"] <- "Smooth Muscle"
id[id =="Smooth Muscle (General Gastrointestinal)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Uterine)"] <- "Smooth Muscle"
id[id =="Tuft Cell"] <- "Tuft Cell"


seurat_object$id <- id

#Create data frame for each summarised id and how many cells they contain
id <- unique(seurat_object$id)
df <- data.frame(row.names = id)
meta <- as.data.frame(seurat_object[[]])
t <- 1


for (i in id){
    temp <- filter(meta, id == i)
    df[i,1] <- length(temp$id)
    
}

colnames(df) <- c("cell_count")

#Filter for cell types (id) with more than 100 cells
filter_id = rownames(filter(df, cell_count > 100 ))

temp <- subset(seurat_object, id %in% filter_id)

Idents(temp) <- seurat_object$id 
DefaultAssay(temp) <- 'chromvar'
ident <- unique(temp$id)
df <- NULL
df <- data.frame(row.names = ident)

for (i in ident){
differential.activity <- FindMarkers(
  object = temp,
  ident.1 = i,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

    
df[t,1] <- rownames(differential.activity)[1]
df[t,2] <- rownames(differential.activity)[2]
df[t,3] <- rownames(differential.activity)[3]
df[t,4] <- rownames(differential.activity)[4]
df[t,5] <- rownames(differential.activity)[5]
t <- t + 1
   
if (FALSE){
mp <- MotifPlot(
  object = seurat_object,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)

path <- paste("report/motif_plot", i, ".pdf", sep = "")
ggsave(mp, 
       filename = path,
       device = "pdf",
       height = 6, width = 8, units = "in") 
    
}
} 

df

hp_featurs <- NULL
for (r in seq(1, nrow(df))){
    hp_featurs <- append(hp_featurs, df[r,1])
    hp_featurs <- append(hp_featurs, df[r,2])
    hp_featurs <- append(hp_featurs, df[r,3])
    hp_featurs <- append(hp_featurs, df[r,4])
    hp_featurs <- append(hp_featurs, df[r,5])    
}

DefaultAssay(temp) <- "chromvar"
Idents(temp) <- temp$id
hp <- DoHeatmap(
  temp,
  features = hp_featurs,
  group.by = "ident",
  slot = "data",
  size = 3,
    
  combine = TRUE
)

hp <- hp + ggtitle("Heatmap Motifs")
ggsave(hp, 
       filename = "report/HeatMap_motifs.pdf",
       device = "pdf",
       height = 10, width = 10, units = "in") 

DefaultAssay(temp) <- "RNA"
df2 <- NULL
df2 <- data.frame(row.names = ident)
t <- 1

for (i in ident){
differential.activity <- FindMarkers(
  object = temp,
  ident.1 = i,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

    
df2[t,1] <- rownames(differential.activity)[1]
df2[t,2] <- rownames(differential.activity)[2]
df2[t,3] <- rownames(differential.activity)[3]
df2[t,4] <- rownames(differential.activity)[4]
df2[t,5] <- rownames(differential.activity)[5]
t <- t + 1
    
}

df2

hp_featurs <- NULL
for (r in seq(1, nrow(df))){
    hp_featurs <- append(hp_featurs, df2[r,1])
    hp_featurs <- append(hp_featurs, df2[r,2])
    hp_featurs <- append(hp_featurs, df2[r,3])
    hp_featurs <- append(hp_featurs, df2[r,4])
    hp_featurs <- append(hp_featurs, df2[r,5])    
}

DefaultAssay(temp) <- "RNA"
Idents(temp) <- temp$id
hp <- DoHeatmap(
  temp,
  features = hp_featurs,
  group.by = "ident",
  slot = "data",
  size = 3,
  disp.max = 2,    
  combine = TRUE
)

ggsave(hp, 
       filename = "report/HeatMap_genes.pdf",
       device = "pdf",
       height = 10, width = 10, units = "in") 

temp_f <- subset(temp, id == "Fibroblast")
temp_f <- subset(temp_f, cell.type %in% c("Fibroblast (Gastrointestinal)", "Fibroblast (Liver Adrenal)", "Fibroblast (General)", "Fibroblast (Peripheral Nerve)")


Idents(temp_f) <- temp_f$cell.type 
DefaultAssay(temp_f) <- 'chromvar'
ident <- unique(temp_f$cell.type)

df3 <- NULL
df3 <- data.frame(row.names = ident)
t <- 1

for (i in ident){
differential.activity <- FindMarkers(
  object = temp_f,
  ident.1 = i,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

    
df3[t,1] <- rownames(differential.activity)[1]
df3[t,2] <- rownames(differential.activity)[2]
df3[t,3] <- rownames(differential.activity)[3]
df3[t,4] <- rownames(differential.activity)[4]
df3[t,5] <- rownames(differential.activity)[5]
t <- t + 1
    
}

df3

hp_featurs <- NULL
for (r in seq(1, nrow(df))){
    hp_featurs <- append(hp_featurs, df3[r,1])
    hp_featurs <- append(hp_featurs, df3[r,2])
    hp_featurs <- append(hp_featurs, df3[r,3])
    hp_featurs <- append(hp_featurs, df3[r,4])
    hp_featurs <- append(hp_featurs, df3[r,5])    
}

DefaultAssay(temp_f) <- "chromvar"
Idents(temp_f) <- temp_f$cell.type
hp <- DoHeatmap(
  temp_f,
  features = hp_featurs,
  group.by = "ident",
  slot = "data",
  size = 3,    
  combine = TRUE
)

ggsave(hp, 
       filename = "report/HeatMap_motifs_fibro.pdf",
       device = "pdf",
       height = 10, width = 10, units = "in") 

fp <- FeaturePlot(
  object = temp,
  features = "MA0071.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.5
)

ggsave(fp, 
       filename = "report/Featureplot_MA0071.1.pdf",
       device = "pdf",
       height = 10, width = 10, units = "in")

###################################
#Do a heat map
###################################
seurat_object <- qread("report/GSE184462_merged_1_20_Motif.rds")

id <- seurat_object$cell.type
id <- as.data.frame(id)
#Assing summariesed id's
id[id =="Adipocyte"] <- "Adipocyte"
id[id =="Alveolar Capillary Endothelial Cell"] <- "Endothelial"
id[id =="Alveolar Type 2 (AT2) Cell"] <- "Alveolar"
id[id =="Alverolar Type 2,Immune"] <- "Alveolar"
id[id =="Cardiac Fibroblasts"] <- "Fibroblast"
id[id =="Cardiac Pericyte 1"] <- "Pericyte"
id[id =="Cardiac Pericyte 2"] <- "Pericyte"
id[id =="Cardiac Pericyte 3"] <- "Pericyte"
id[id =="Cardiac Pericyte 4"] <- "Pericyte"
id[id =="Chief Cell"] <- "Chief Cell"
id[id =="CNS,Enteric Neuron"] <- "Neuron"
id[id =="Cortical Epithelial-like"] <- "Epithelial"
id[id =="Ductal Cell (Pancreatic)"] <- "Ductal Cell"
id[id =="Endocardial Cell"] <- "Endocardial"
id[id =="Endothelial (Exocrine Tissues)"] <- "Endothelial"
id[id =="Endothelial Cell (General) 1"] <- "Endothelial"
id[id =="Endothelial Cell (General) 2"] <- "Endothelial"
id[id =="Endothelial Cell (General) 3"] <- "Endothelial"
id[id =="Endothelial Cell (Myocardial)"] <- "Endothelial"
id[id =="Fibroblast (Epithelial)"] <- "Fibroblast"
id[id =="Fibroblast (General)"] <- "Fibroblast"
id[id =="Fibroblast (Liver Adrenal)"] <- "Fibroblast"
id[id =="Fibroblast (Peripheral Nerve)"] <- "Fibroblast"
id[id =="Foveolar Cell"] <- "Foveolar Cell"
id[id =="Keratinocyte 1"] <- "Keratinocyte"
id[id =="Luteal Cell (Ovarian)"] <- "Luteal Cell"
id[id =="Lymphatic Endothelial Cell"] <- "Endothelial"
id[id =="Macrophage (General,Alveolar)"] <- "Macrophage"
id[id =="Macrophage (General)"] <- "Macrophage"
id[id =="Mast Cell"] <- "Mast Cell"
id[id =="Memory B Cell"] <- "B"
id[id =="Mesothelial Cell"] <- "Mesothelial"
id[id =="Naive T cell"] <- "T"
id[id =="Natural Killer T Cell"] <- "T"
id[id =="Pancreatic Acinar Cell"] <- "Acinar Cell"
id[id =="Pericyte (Esophageal Muscularis)"] <- "Pericyte"
id[id =="Pericyte (General) 1"] <- "Pericyte"
id[id =="Pericyte (General) 2"] <- "Pericyte"
id[id =="Pericyte (General) 3"] <- "Pericyte"
id[id =="Pericyte (General) 4"] <- "Pericyte"
id[id =="Peripheral Nerve Stromal"] <- "Neuron"
id[id =="Plasma Cell"] <- "Plasma Cell"
id[id =="Schwann Cell (General)"] <- "Schwann Cell"
id[id =="Small Intestinal Enterocyte"] <- "Enterocyte"
id[id =="Smooth Muscle (Esophageal Muscularis) 3"] <- "Smooth Muscle"
id[id =="Smooth Muscle (GE Junction)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (General)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Vaginal)"] <- "Smooth Muscle"
id[id =="T Lymphocyte 1 (CD8+)"] <- "T"
id[id =="T lymphocyte 2 (CD4+)"] <- "T"
id[id =="Transitional Zone Cortical Cell"] <- "Cortical Cell"
id[id =="Vascular Smooth Muscle 1"] <- "Smooth Muscle"
id[id =="Vascular Smooth Muscle 2"] <- "Smooth Muscle"
id[id =="Zona Fasciculata Cortical Cell"] <- "Cortical Cell"
id[id =="Zona Glomerulosa Cortical Cell"] <- "Cortical Cell"

#Additional Cellt types (11-20)
id[id =="Airway Goblet Cell"] <- "Goblet Cell"
id[id =="Astrocyte 2"] <- "Astrocyte"
id[id =="Basal Epidermal (Skin)"] <- "Epidermal"
id[id =="Basal Epithelial (Mammary)"] <- "Epithelial"
id[id =="Blood Brain Barrier Endothelial Cell"] <- "Endothelial"
id[id =="Club Cell"] <- "Club Cell"
id[id =="Colon Epithelial Cell 1"] <- "Epithelial"
id[id =="Colon Epithelial Cell 2"] <- "Epithelial"
id[id =="Colon Epithelial Cell 3"] <- "Epithelial"
id[id =="Colonic Goblet Cell"] <- "Goblet Cell"
id[id =="Enterochromaffin Cell"] <- "Enterochromaffin Cell"
id[id =="Esophageal Epithelial Cell"] <- "Epithelial"
id[id =="Fibroblast (Gastrointestinal)"] <- "Fibroblast"
id[id =="Fibroblast (Sk Muscle Associated)"] <- "Fibroblast"
id[id =="GABAergic Neuron 1"] <- "Neuron"
id[id =="Gastric Neuroendocrine Cell"] <- "Neuron"
id[id =="Granular Epidermal (Skin)"] <- "Epidermal"
id[id =="Mammary Luminal Epithelial Cell 2"] <- "Epithelial"
id[id =="Melanocyte"] <- "Melanocyte"
id[id =="Myoepithelial (Skin)"] <- "Myoepithelial"
id[id =="Oligodendrocyte"] <- "Oligodendrocyte"
id[id =="Oligodendrocyte Precursor"] <- "Oligodendrocyte"
id[id =="Paneth Cell"] <- "Paneth Cell"
id[id =="Small Intestinal Goblet Cell"] <- "Adipocyte"
id[id =="Smooth Muscle (Colon) 1"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Colon) 2"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Mucosal)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Muscularis) 1"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Esophageal Muscularis) 2"] <- "Smooth Muscle"
id[id =="Smooth Muscle (General Gastrointestinal)"] <- "Smooth Muscle"
id[id =="Smooth Muscle (Uterine)"] <- "Smooth Muscle"
id[id =="Tuft Cell"] <- "Tuft Cell"


seurat_object$id <- id

Idents(seurat_object) <- seurat_object$id
id <- Idents(seurat_object)

d <- unique(seurat_object$id)
df <- data.frame(row.names = id)
meta <- as.data.frame(seurat_object[[]])
df <- data.frame(row.names = ident)
t <- 1


for (i in id){
    temp <- filter(meta, id == i)
    df[i,1] <- length(temp$id)
    
}

colnames(df) <- c("cell_count")

#Filter for cell types (id) with more than 100 cells
filter_id = rownames(filter(df, cell_count > 100 ))

temp <- subset(seurat_object, id %in% filter_id)
Idents(temp) <- temp$id

# gather the footprinting information for sets of motifs
DefaultAssay(temp) <- "peaks"
seurat_object <- Footprint(
  object = temp,
  motif.name = c("RORA"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

Idents(seurat_object) <- seurat_object$id
# plot the footprint data for each group of cells
fp <- PlotFootprint(seurat_object, features = c("RORA"), label = FALSE,
                    idents = c("Mesothelial", "Fibroblast", "Endothelial", "T", "Cortical Cell" ))


ggsave(fp, 
       filename = "report/TF_Footprint_RORA.pdf",
       device = "pdf",
       height = 5, width = 10, units = "in")
