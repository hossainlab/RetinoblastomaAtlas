# Retinoblastoma Atlas: R Pipeline Documentation

This directory contains a clean, concise, and high-performance R pipeline for single-cell RNA-seq analysis of Retinoblastoma, optimized for a **32GB RAM** environment.

## 🚀 Workflow Overview

The pipeline is split into 6 modular scripts. Each script saves an intermediate RDS file to `data/processed/` to ensure you can resume at any point.

| Order | Script | Description | Key Tools |
| :--- | :--- | :--- | :--- |
| 1 | `01_data_loading_qc.R` | Loads 140k cells, performs adaptive QC & doublet detection. | Seurat, DoubletFinder |
| 2 | `02_normalization_integration.R` | Normalizes data and removes batch effects between datasets. | Harmony, Seurat v5 |
| 3 | `03_annotation_markers.R` | Identifies cluster markers and annotates cell types. | Seurat, scCustom |
| 4 | `04_scoring_pathway.R` | Scores tumor subtypes and TGF-B pathway activity. | decoupleR, PROGENy |
| 5 | `05_cell_communication.R` | Infers ligand-receptor networks comparing Intra vs Extraocular. | CellChat |
| 6 | `06_trajectory_analysis.R` | Maps cone precursor evolution toward invasive states. | Slingshot |

## 🛠️ Performance Tips for 32GB RAM

- **Garbage Collection:** Scripts explicitly call `gc()` and delete large intermediate objects.
- **Harmony Integration:** Used instead of Seurat's `IntegrateData` to minimize memory footprint.
- **Parallelization:** You can enable multi-threading for certain steps by running `plan("multisession", workers = 4)` at the start of scripts (requires `future` package).

## 📊 Publication-Ready Outputs

- **Figures:** Saved in `results/figures/` as high-DPI PNGs/PDFs using `scCustom` and `ggplot2`.
- **Tables:** Marker lists and scoring metrics saved in `results/tables/`.

## 📦 Required Packages

Ensure you have the following installed in R:
```R
install.packages(c("Seurat", "tidyverse", "patchwork", "harmony", "grDevices", "BiocManager"))
BiocManager::install(c("DoubletFinder", "decoupleR", "CellChat", "slingshot", "SingleCellExperiment"))
remotes::install_github("samuel-marsh/scCustom")
```
