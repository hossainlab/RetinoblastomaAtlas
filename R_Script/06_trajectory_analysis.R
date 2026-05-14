# Script 06: Trajectory Inference (Slingshot)
# Project: RetinoblastomaAtlas
# Author: Md. Jubayer Hossain
# Date: 2026-05-05

library(Seurat)
library(slingshot)
library(tidyverse)
library(grDevices)

# 1. Load Data
merged_obj <- readRDS("data/processed/final_scored.rds")
cone_obj <- subset(merged_obj, idents = "Cone_Precursor")

# 2. Run Slingshot
message("Running Slingshot trajectory on Cone Precursors...")
sce <- as.SingleCellExperiment(cone_obj)
sce <- slingshot(sce, clusterLabels = 'ident', reducedDim = 'UMAP')

# 3. Visualization
png("results/figures/slingshot_trajectory_umap.png", width = 800, height = 800, res = 150)
plot(reducedDims(sce)$UMAP, col = rainbow(10)[as.numeric(as.factor(sce$ident))], 
     pch = 16, cex = 0.5, main = "Cone Precursor Trajectory")
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
dev.off()

# 4. Pseudotime Driver Genes
pseudotime <- colData(sce)$slingPseudotime_1
gene_expr <- as.matrix(logcounts(sce))
cor_results <- apply(gene_expr, 1, function(x) cor(x, pseudotime, use = "complete.obs"))
drivers <- sort(cor_results, decreasing = TRUE)

write_csv(data.frame(gene = names(drivers), correlation = drivers), 
          "results/tables/trajectory_driver_genes.csv")
