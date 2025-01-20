library(Seurat)
library(tidyverse)
library(harmony)
library(Matrix)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(dplyr)
theme_set(theme_cowplot())
set.seed(12345)

# Setting up the file directory
sample_dirs <- ('./Isoform_Seurat/Combined_ensemble_gtf/')
names(sample_dirs) <- c('VR03')

# Import expression count in Matrix Market format
iso_X_list <- lapply(sample_dirs, function(x){Seurat::ReadMtx(paste0(x, 'isoforms_seurat/matrix.mtx'), 
                                                              paste0(x, 'isoforms_seurat/barcodes_RC.tsv'), 
                                                              paste0(x, 'isoforms_seurat/genes.tsv'))})
gene_X_list <- lapply(sample_dirs, function(x){Seurat::ReadMtx(paste0(x, 'genes_seurat/matrix.mtx'), 
                                                               paste0(x, 'genes_seurat/barcodes_RC.tsv'), 
                                                               paste0(x, 'genes_seurat/genes.tsv'))})
seurat_list <- lapply(1:length(sample_dirs), function(i){
  sobj <- Seurat::CreateSeuratObject(gene_X_list[[i]], meta.data = VR03all_metadata, min.cells = 1)
  sobj$barcode <- colnames(sobj)
  sobj$Sample <- names(sample_dirs)[[i]]
  sobj[["iso"]] <- Seurat::CreateAssayObject(counts = iso_X_list[[i]], min.cells = 1)
  sobj
})

# Combine both gene-level and isoform-level count into single seurat object
VR03all_seurat <- Reduce(merge, seurat_list)

# remove spatial spots that are on the irrelevent tissue domain 
shoot_list <- rownames(subset(VR03all_seurat@meta.data, domains=="shoot"))
VR03all_seurat <- VR03all_seurat[, !colnames(VR03all_seurat) %in% shoot_list]
root_list <- rownames(subset(VR03all_seurat@meta.data, domains=="root"))
VR03all_seurat <- VR03all_seurat[, !colnames(VR03all_seurat) %in% root_list]
features <- c('nCount_RNA', 'nFeature_RNA', 'nCount_iso', 'nFeature_iso')

# Processing gene-level expression using standard Seurat workflow
DefaultAssay(VR03all_seurat) <- 'RNA'
VR03all_seurat <- VR03all_seurat %>%
  NormalizeData() %>%
  FindVariableFeature(selection.method = 'vst', nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, reduction.name = 'pca_RNA')
# Batch-correction using the four tissue slices as variable
VR03all_seurat <- RunHarmony(VR03all_seurat, group.by.cars = 'Section', reduction = 'pca_RNA', reduction.save = 'harmony_RNA', assay.use = 'RNA')
# Unsupervised clustering based on the batch-corrected UMAP
VR03all_seurat <- VR03all_seurat %>%
  RunUMAP(min.dist=0.3, reduction = 'harmony_RNA', reduction.name = 'umap_RNA') %>%
  FindNeighbors(dims = 1:30, reduction = 'harmony_RNA') %>%
  FindClusters(cluster.name = 'suerat_cluster_RNA')

# Processing isoform-level expression using standard Seurat workflow
DefaultAssay(VR03all_seurat) <- 'iso'
VR03all_seurat <- VR03all_seurat %>%
  NormalizeData() %>%
  FindVariableFeature(selection.method = 'vst', nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, reduction.name = 'pca_iso')
# Batch-correction using the four tissue slices as variable
VR03all_seurat <- RunHarmony(VR03all_seurat, group.by.cars = 'Section', reduction = 'pca_iso', reduction.save = 'harmony_iso', assay.use = 'iso')
# Unsupervised clustering based on the batch-corrected UMAP
VR03all_seurat <- VR03all_seurat %>%
  RunUMAP(min.dist=0.3, reduction = 'harmony_iso', reduction.name = 'umap_iso') %>%
  FindNeighbors(dims = 1:30, reduction = 'harmony_iso') %>%
  FindClusters(cluster.name = 'seurat_cluster_iso')

# Import tissue histology image into seurat
library(stringr)
library(Biostrings)
library(STutility)

# Histology image informations were retrieved from short-read Space Ranger output
VR03_samples <- "./Seurat_rds/XGE21_VR03/outs/spatial/filtered_feature_bc_matrix.h5"
VR03_imgs <- "./Seurat_rds/XGE21_VR03/outs/spatial/tissue_hires_image.png"
VR03_spotfiles <- "./Seurat_rds/XGE21_VR03/outs/spatial/tissue_positions_list.csv"
VR03_json <- "./Seurat_rds/XGE21_VR03/outs/spatial/scale_factors_json.json"
# process the image information, so the information can be transfer to long-read dataset
VR03_infoTable <- data.frame(VR03_samples, VR03_imgs, VR03_spotfiles, VR03_json)
VR03_SE <- InputFromTable(VR03_infoTable)
VR03_SE <- VR03_SE %>%
  LoadImages(time.resolve = FALSE) %>%
  ManualAnnotation()
transform <- list("1"  = list("shift.x" = -20, "shift.y" = -10))
VR03_SE <- VR03_SE %>%
  MaskImages() %>%
  WatpImages(transform)
VR03all_seurat@tools$Staffli <- GetStaffli(VR03_SE)
VR03all_seurat@tools$Staffli@meta.data$barcode <- rownames(VR03all_seurat@tools$Staffli@meta.data)
# Plotting the cluster domains from unsupervised clustering on histology image
FeatureOverlay(VR03all_seurat, features = 'seurat_cluster_RNA' pt.size = 3, show.sb = FALSE, pt.alpha = TRUE) | 
  FeatureOverlay(VR03all_seurat, features = 'seurat_cluster_iso', pt.size = 3, show.sb = FALSE, pt.alpha = TRUE)

