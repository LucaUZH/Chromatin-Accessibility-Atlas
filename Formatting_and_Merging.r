## first line of code is to clear R's memory
rm(list=ls())

#libraries
library(Signac)
library(Seurat)
library(rtracklayer)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(qs)
library(GenomicRanges)
library(future)


#Adult datasets
file1 = "GSM5589344_adipose_omentum_SM-ADYHB_rep1_fragments"
file2 = "GSM5589345_adipose_omentum_SM-CHZRM_rep1_fragments"
file3 = "GSM5589346_adipose_omentum_SM-CSSD4_rep1_fragments"
file4 = "GSM5589347_adipose_omentum_SM-IOBHJ_rep1_fragments"
file5 = "GSM5589348_adrenal_gland_SM-A8WNO_rep1_fragments"
file6 = "GSM5589349_artery_aorta_SM-C1MLC_rep1_fragments"
file7 = "GSM5589350_artery_aorta_SM-C1PX3_rep1_fragments"
file8 = "GSM5589351_artery_aorta_SM-CR89M_rep1_fragments"
file9 = "GSM5589352_artery_aorta_SM-JF1NU_rep1_fragments"
file10 = "GSM5589353_artery_tibial_SM-CHLWW_rep1_fragments"
file11 = "GSM5589355_colon_sigmoid_SM-AZPYO_rep1_fragments"
file12 = "GSM5589356_colon_sigmoid_SM-JF1O8_rep1_fragments"
file13 = "GSM5589357_colon_transverse_SM-A9HOW_rep1_fragments"
file14 = "GSM5589358_colon_transverse_SM-A9VP4_rep1_fragments"
file15 = "GSM5589359_colon_transverse_SM-ACCQ1_rep1_fragments"
file16 = "GSM5589360_colon_transverse_SM-BZ2ZS_rep1_fragments"
file17 = "GSM5589361_colon_transverse_SM-CSSDA_rep1_fragments"
file18 = "GSM5589362_esophagus_ge_junction_SM-CTD24_rep1_fragments"
file19 = "GSM5589363_esophagus_ge_junction_SM-IOERG_rep1_fragments"
file20 = "GSM5589364_esophagus_mucosa_SM-A9HOR_rep1_fragments"
file21 = "GSM5589365_esophagus_mucosa_SM-A9VPA_rep1_fragments"
file22 = "GSM5589366_esophagus_mucosa_SM-AZPYJ_rep1_fragments"
file23 = "GSM5589367_esophagus_muscularis_SM-A8CPH_rep1_fragments"
file24 = "GSM5589368_esophagus_muscularis_SM-CSSCV_rep1_fragments"
file25 = "GSM5589369_esophagus_muscularis_SM-IOBHM_rep1_fragments"
file26 = "GSM5589370_esophagus_muscularis_SM-IQYD1_rep1_fragments"
file27 = "GSM5589371_heart_atrial_appendage_SM-IOBHN_rep1_fragments"
file28 = "GSM5589372_heart_atrial_appendage_SM-JF1NX_rep1_fragments"
file29 = "GSM5589373_heart_lv_SM-IOBHO_rep1_fragments"
file30 = "GSM5589374_heart_lv_SM-JF1NY_rep1_fragments"
file31 = "GSM5589375_liver_SM-A8WNZ_rep1_fragments"
file32 = "GSM5589376_lung_SM-A62E9_rep1_fragments"
file33 = "GSM5589377_lung_SM-A8WNH_rep1_fragments"
file34 = "GSM5589378_lung_SM-ACCPU_rep1_fragments"
file35 = "GSM5589379_lung_SM-JF1NZ_rep1_fragments"
file36 = "GSM5589380_mammary_tissue_SM-IOBHL_rep1_fragments"
file37 = "GSM5589381_mammary_tissue_SM-JF1NV_rep1_fragments"
file38 = "GSM5589382_muscle_SM-ADA6L_rep1_fragments"
file39 = "GSM5589383_muscle_SM-C1MKW_rep1_fragments"
file40 = "GSM5589384_muscle_SM-C1PWV_rep1_fragments"
file41 = "GSM5589385_muscle_SM-IOBHP_rep1_fragments"
file42 = "GSM5589386_muscle_SM-JF1O9_rep1_fragments"
file43 = "GSM5589387_nerve_tibial_SM-CHLWU_rep1_fragments"
file44 = "GSM5589388_nerve_tibial_SM-CP2V6_rep1_fragments"
file45 = "GSM5589389_nerve_tibial_SM-IOBHQ_rep1_fragments"
file46 = "GSM5589390_ovary_SM-IOBHR_rep1_fragments"
file47 = "GSM5589391_pancreas_SM-ADRUQ_rep1_fragments"
file48 = "GSM5589392_pancreas_SM-IOBHS_rep1_fragments"
file49 = "GSM5589393_pancreas_SM-JF1NS_rep1_fragments"
file50 = "GSM5589394_pancreas_SM-JF1O6_rep1_fragments"
file51 = "GSM5589395_skin_SM-IOBHT_rep1_fragments"
file52 = "GSM5589396_skin_SM-JF1O1_rep1_fragments"
file53 = "GSM5589397_skin_sun_exposed_SM-ADYHK_rep1_fragments"
file54 = "GSM5589398_skin_sun_exposed_SM-IOBHU_rep1_fragments"
file55 = "GSM5589399_skin_sun_exposed_SM-IQYCP_rep1_fragments"
file56 = "GSM5589400_skin_sun_exposed_SM-JF1NT_rep1_fragments"
file57 = "GSM5589401_small_intestine_SM-A62GO_rep1_fragments"
file58 = "GSM5589402_small_intestine_SM-ADA5F_rep1_fragments"
file59 = "GSM5589403_small_intestine_SM-JF1O2_rep1_fragments"
file60 = "GSM5589404_stomach_SM-CHLWL_rep1_fragments"
file61 = "GSM5589405_stomach_SM-IOBHV_rep1_fragments"
file62 = "GSM5589406_stomach_SM-JF1NP_rep1_fragments"
file63 = "GSM5589407_stomach_SM-JF1O3_rep1_fragments"
file64 = "GSM5589408_thyroid_SM-C1MKY_rep1_fragments"
file65 = "GSM5589409_thyroid_SM-IOBHW_rep1_fragments"
file66 = "GSM5589410_thyroid_SM-JF1O4_rep1_fragments"
file67 = "GSM5589411_uterus_SM-A87A2_rep1_fragments"
file68 = "GSM5589412_uterus_SM-IOBHX_rep1_fragments"
file69 = "GSM5589413_vagina_SM-A9HOS_rep1_fragments"
file70 = "GSM5589414_UMB4540_snATAC_frontal_cortex_rep1_fragments"
file71 = "GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments"

lyst_adult = list(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12, file13,file14, file15, file16,
           file17, file18, file19, file20, file21, file22, file23, file24, file25, file26, file27, file28, file29, file30, 
           file31, file32, file33, file34, file35, file36, file37, file38, file39, file40, file41, file42, file43, file44, file45,
           file46, file47, file48, file49, file50, file51, file52, file53, file54, file55, file56, file57, file58, file59, file60,
           file61, file62, file63, file64, file65, file66, file67, file68, file69, file70, file71)

#Fetal datasets
file80 = "GSM4508928_adrenal_filtered.seurat"
file81 = "GSM4508929_cerebellum_filtered.seurat"
file82 = "GSM4508930_cerebrum_filtered.seurat"
file83 = "GSM4508931_eye_filtered.seurat"
file84 = "GSM4508932_heart_filtered.seurat"
file85 = "GSM4508933_intestine_filtered.seurat"
file86 = "GSM4508934_kidney_filtered.seurat"
file87 = "GSM4508935_liver_filtered.seurat"
file88 = "GSM4508936_lung_filtered.seurat"
file89 = "GSM4508937_muscle_filtered.seurat"
file90 = "GSM4508938_pancreas_filtered.seurat"
file91 = "GSM4508939_placenta_filtered.seurat"
file92 = "GSM4508940_spleen_filtered.seurat"
file93 = "GSM4508941_stomach_filtered.seurat"
file94 = "GSM4508942_thymus_filtered.seurat"

lyst_fetal = list(file80, file81, file82, file83, file84, file85, file86, file87, file88, file89, file90, file91, file92, 
                  file93, file94)


################################################
#Set cores and memory
################################################
plan("multicore", workers = 16)
options(future.globals.maxSize = 150000 * 1024^2) # for 150 Gb RAM

################################################
#BED to TSV (Fragment) files
################################################
for (f in lyst_adult) {  
    #Get file path for the bed file
    peak_file = paste("GSE184462_extr/",  f, ".bed", sep = "")
    #set file path for saving the fragments file
    tsv_file = paste("GSE184462_extr/fragments/", f, ".tsv", sep = "")

    #Read bed files
    peaks <- import(peak_file)

    #Rearraanging columns V1= chrom, V2 = frag start, V3 = end frag, V4 = barcode, V5 = # of reads 
    chrom <- as.vector(seqnames(peaks))
    start <- start(ranges(peaks))
    end   <- end(ranges(peaks))
    barcode  <- peaks$name 
    readCount<- peaks$score 
    
    #Create dataframe and save as TSV file
    dd <- data.frame(chrom, start, end, barcode, readCount)
    readr::write_tsv(dd,tsv_file, col_names=FALSE)
} 


################################################
#TSV files have to be sorted, gzipped and indexed in Ubuntu
################################################

################################################
#Create seurat Objects from fragment & metadata
################################################

#Read metadata
metadata <- read_tsv(
  file = "GSE184462_metadata.tsv"
)

for (f in lyst_adult) {
    #Get file pahts
    peak_file = paste("GSE184462_extr/",  f, ".bed", sep = "")
    tsv_file = paste("GSE184462_extr/fragments/", f, ".sort.tsv.gz", sep = "")
    #Set file path for new seurat object
    seurat_file = paste("GSE184462_extr/seurat/meta/", f, "_meta", ".rds", sep = "")

    #Read bed file to get cells (barcodes)
    peaks <- import(peak_file)    
    cells_from_peaks <- unique(peaks$name)

    #Create a Fragment object 
    fragments <- CreateFragmentObject(
      path = tsv_file,
      cells = cells_from_peaks, 
      validate.fragments = TRUE
    )

    #Call peaks using MACS2 using the Fragments
    peak_range <- CallPeaks(
      object = fragments,
      group.by = "predicted.id"
    )

    #Create peak-by-cell matrix
    peak_matrix <- FeatureMatrix(
      fragments = fragments,
      features = peak_range
    )
 

    #Convert peak-by-cell matrix to chromatin assay object
    chrom_assay <- CreateChromatinAssay(
      counts = peak_matrix,
      sep = c(":", "-"),
      fragments = tsv_file
    )

    #Hard coding tissue_name in metadata for brain datasets
    tissue_name <- gsub("\\_rep.*", "", substr(f, 12, 100))
    if (f == file70){
        tissue_name <- "Human_brain_1"
    }

    if (f == file71){
        tissue_name <- "Human_brain_2"
    }
    
    #Filte metadata for cells relevant for the dataset by looking at the tissue_name
    df <- as.data.frame(metadata, , row.names = 1)
    df2 <- filter(df, tissue == tissue_name)
    df2 <- mutate(df2, cellID = str_remove(cellID, sample))
    df2 <- mutate(df2, cellID = str_remove(cellID, fixed("+")))
    rownames(df2) <- df2$cellID
    df2 <- select(df2, -1)
    
    #Create Seurat objects
    SeruatObj <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "peaks",
      meta.data = df2
    )
    
    #Save it to designated path
    saveRDS(SeruatObj, seurat_file)
    
} 

################################################
#Merging adult datasets
################################################

#Lists containing the seurat objects, & their id
seurat_list <- NULL
cell_id <- NULL

for (f in lyst_adult){
    #Get file path for seurat objects
    seurat_file = paste("GSE184462_extr/seurat/meta/", f, "_meta", ".rds", sep = "")
 
    #Read seurat files & filter for cells present in metadata
    seurat_object <- readRDS(seurat_file)
    seurat_object <- subset( x = seurat_object, subset = tsse > 0)

    #Add every seurat object & its id to the lists
    if (f == file1){
        seurat_list <- list(seurat_object)
        cell_id <- c(f)
    }
    else{
        seurat_list <- append(seurat_list, seurat_object)
        cell_id <- append(cell_id, f)
    }
}

#combine the seurat object form the list
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = cell_id
)

#Save the combined files
qsave(combined, "complete/GSE184462_merged.rds")


################################################
#Change metadata of fetal seurat objects to adult metadata
################################################

#Read metadata
metadata_184462 <- read_tsv(
  file = "GSE184462_metadata.tsv"
)

#Get tissues_names in from the metadata
tissue_meta <- unique(metadata_184462$tissue)

for (f in lyst_fetal) {
    #Get file paths for old seurat object
    seurat_file_fetal <- paste("GSE149683_extr/", f, ".RDS", sep = "")  
    #Set file path for new seurat object
    seurat_file_fetal_save <- paste("GSE149683_extr/", f, "_meta.rds", sep = "")
    
    #Read old seurat object
    seurat_object = readRDS(seurat_file_fetal) 

    #Extract tissues form file names
    tissue_name <- gsub("\\_filtered.*", "", substr(f, 12, 100))
    tissue_name_list <- NULL

    #Search for tisse_names matching with the tissue extracted from the file name
    for (t in tissue_meta){
        if (str_detect(t, tissue_name) && str_detect(t, "sample")){
            #Save all tissue_names matching with dataset
            tissue_name_list <- append(tissue_name_list, t)
        }
    }

    #Filter metadata for relevant cells based on tissue_names
    df <- as.data.frame(metadata_184462, , row.names = 1)
    df2 <- filter(df, tissue %in% tissue_name_list)    
    df2 <- mutate(df2, cellID = str_remove(cellID, sample))
    df2 <- mutate(df2, cellID = str_remove(cellID, fixed("+")))
    rownames(df2) <- df2$cellID
    df2 <- select(df2, -1)

    #Extract data from old seurat object
    DefaultAssay(seurat_object) <- "peaks"
    data <- GetAssayData(seurat_object, slot = "data")
    chromatinassay <- CreateChromatinAssay(counts = data, genome = "hg38")
    
    #Generate new seurat object with new metadata
    new_seurat <- CreateSeuratObject(counts = chromatinassay, assay = "peaks",  meta.data = df2)    

    #Save new seurat object
    saveRDS(new_seurat, seurat_file_fetal_save)
}

################################################
#Merging fetal datasets
################################################

#Lists for seurat objects & their ids
new_seurat_list <- NULL
cell_id <- NULL

for (f in lyst_fetal){
    #Get file path for seurat object
    seurat_file_fetal <- paste("GSE149683_extr/", f, "_meta.rds", sep = "")
    
    #Read seurat object
    seurat_object = readRDS(seurat_file_fetal) 
    
    #Add all seuruat object & ids to the lists
    if (f == file80){
        new_seurat_list <- list(new_seurat)
        cell_id <- c(f)
    }
    else{
        new_seurat_list <- append(new_seurat_list, new_seurat)
        cell_id <- append(cell_id, f)
    }
    
}

#Combine the seurat object form the list
combined <- merge(
  x = new_seurat_list[[1]],
  y = new_seurat_list[-1],
  add.cell.ids = cell_id
)

#Save the combined objects
qsave(combined, "complete/GSE149683_merged_.rds")

################################################
#Merging complete datasets
################################################
seurat_adult <- qread("complete/GSE184462_merged.rds")
seurat_fetal <- qread("complete/GSE149683_merged_.rds")

#Combine both seurat object 
combined <- merge(
  x = seurat_adult,
  y = seurat_fetal,
  add.cell.ids = c("adult", "fetal")
)

#Save the combined objects
qsave(combined, "complete/GSE149683_merged_.rds")

f <- file80
seurat_file_fetal <- paste("GSE149683_extr/", f, ".RDS", sep = "")
seurat_object = readRDS(seurat_file_fetal) 


seurat_object
