# Script 02: Normalization, Integration, and Clustering
# Project: RetinoblastomaAtlas
# Author: Md. Jubayer Hossain 
# Date: 2026-05-05

library(Seurat)
options(Seurat.warn.umap.method = FALSE)
library(harmony)
library(tidyverse)
library(scCustom)

# 1. Load Data
merged_obj <- readRDS("data/processed/qc_filtered.rds")
message("Atlas loaded with ", ncol(merged_obj), " cells.")

# 2. Standard Workflow
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 3000)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = FALSE)

# 3. Harmony Integration
message("Running Harmony integration...")
merged_obj <- RunHarmony(merged_obj, group.by.vars = c("dataset", "patient_id"), 
                         plot_convergence = FALSE, dims.use = 1:30)

# 4. Dimensionality Reduction and Clustering
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:30, umap.method = 'uwot')
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = c(0.2, 0.5, 0.8))

# Rename resolution column for underscore consistency
merged_obj$RNA_snn_res_0_5 <- merged_obj$RNA_snn_res.0.5

# 5. Visualization
p1 <- DimPlot(merged_obj, reduction = "umap", group.by = "dataset") + 
  theme_minimal() + ggtitle("Integration by Dataset")
p2 <- DimPlot(merged_obj, reduction = "umap", group.by = "disease_stage") + 
  theme_minimal() + ggtitle("Integration by Disease Stage")

ggsave("results/figures/umap_integration_check.png", p1 + p2, width = 14, height = 6, dpi = 300)

SetIdent(merged_obj) <- "RNA_snn_res_0_5"
p3 <- DimPlot_scCustom(merged_obj, reduction = "umap", label = TRUE) + 
  theme_minimal() + ggtitle("Leiden Clusters (Res 0.5)")
ggsave("results/figures/umap_clusters_res0.5.png", p3, width = 8, height = 7, dpi = 300)

# 6. Save Integrated Object
saveRDS(merged_obj, "data/processed/integrated.rds")