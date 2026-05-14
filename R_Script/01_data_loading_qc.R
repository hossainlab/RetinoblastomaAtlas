# Script 01: Data Loading and Quality Control (Updated)
# Project: RetinoblastomaAtlas
# Author: Gemini CLI (Interactive Agent)
# Date: 2026-05-07

library(Seurat)
options(Seurat.warn.umap.method = FALSE)
library(tidyverse)
library(patchwork)
library(DoubletFinder)

# 1. Setup Paths and Metadata
metadata   <- read_csv("results/tables/sample_metadata.csv")
data_dir   <- "data/raw"
output_dir <- "results/figures"
dir.create("data/processed", showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)

# 2. Sequential Loading and Per-Sample QC
seurat_list <- list()

for (i in 1:nrow(metadata)) {
  sample_id <- metadata$sample_id[[i]]
  dataset   <- metadata$dataset[[i]]
  message("\n>>> Processing sample [", i, "/", nrow(metadata), "]: ", sample_id)
  
  # Load 10x data
  path <- file.path(data_dir, dataset, sample_id)
  
  if (!dir.exists(path)) {
    message("!!! Directory not found at ", path, " - Skipping.")
    next
  }
  
  counts <- Read10X(data.dir = path)
  obj    <- CreateSeuratObject(counts = counts, project = sample_id)
  
  # Add ALL metadata from harmonized table
  obj$dataset       <- dataset
  obj$disease_stage <- metadata$disease_stage[[i]]
  obj$patient_id    <- metadata$patient_id[[i]]
  obj$tissue        <- metadata$tissue[[i]]
  obj$sex           <- metadata$sex[[i]]
  obj$age           <- metadata$age[[i]]
  obj$rb1_mutation  <- metadata$rb1_mutation[[i]]
  obj$geo_accession <- metadata$geo_accession[[i]]
  obj$replicate     <- metadata$replicate[[i]]
  
  # Calculate MT percentage
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Save raw QC plot before filtering
  p_raw <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
                   ncol = 3, pt.size = 0.1) + plot_annotation(title = paste(sample_id, "Raw QC"))
  ggsave(file.path(output_dir, paste0(sample_id, "_qc_vln_raw.png")), p_raw, width = 10, height = 5)

  # Adaptive Filtering (as per GEMINI.md core mandates)
  # Lower bounds: nFeature > 300, nCount > 500
  # Upper bounds: nFeature < 7500, MT < 20%
  obj <- subset(obj, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 &
                  percent_mt < 20 & nCount_RNA > 500)
  
  if (ncol(obj) < 100) {
    message("!!! Too few cells remaining after QC (", ncol(obj), ") - Skipping.")
    next
  }

  # 3. Preprocessing for DoubletFinder
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  # Join layers for DoubletFinder compatibility with Seurat v5
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  
  # pK Identification
  sweep_res   <- paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn       <- find.pK(sweep_stats)
  pk          <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Homotypic Doublet Proportion Estimate
  annotations    <- obj@active.ident
  homotypic_prop <- modelHomotypic(annotations)
  nExp_poi       <- round(0.075 * ncol(obj)) # Assuming 7.5% doublet rate
  nExp_poi_adj   <- round(nExp_poi * (1 - homotypic_prop))
  
  if (nExp_poi_adj == 0) nExp_poi_adj <- 1
  
  # Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = pk,
                       nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)
  
  # Remove Doublets
  doublet_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  obj$doublet_status <- obj@meta.data[[doublet_col]]
  obj <- subset(obj, subset = doublet_status == "Singlet")
  
  message("--- Remaining cells after QC and DoubletFinder: ", ncol(obj))
  
  seurat_list[[sample_id]] <- obj
  
  # Cleanup per iteration to save RAM (32GB limit)
  rm(counts, sweep_res, sweep_stats, bcmvn, obj)
  gc()
}

# 4. Merge and Save
if (length(seurat_list) > 1) {
  message("\n>>> Merging all samples...")
  merged_obj <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
                      add.cell.ids = names(seurat_list), project = "RBAtlas")
} else {
  merged_obj <- seurat_list[[1]]
}

# 5. Final Global QC Plots
p1 <- VlnPlot(merged_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
              group.by = "disease_stage", ncol = 3, pt.size = 0)
ggsave(file.path(output_dir, "qc_vlnplot_final_by_stage.png"),
       p1, width = 12, height = 6, dpi = 300)

p2 <- VlnPlot(merged_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
              group.by = "dataset", ncol = 3, pt.size = 0)
ggsave(file.path(output_dir, "qc_vlnplot_final_by_dataset.png"),
       p2, width = 12, height = 6, dpi = 300)

# 6. Save Processed Object
message("\n>>> Saving integrated object...")
saveRDS(merged_obj, "data/processed/01_qc_filtered.rds")
message("Done. Saved to data/processed/01_qc_filtered.rds")
