# Script 04: Subtype Scoring and Pathway Analysis
# Project: RetinoblastomaAtlas
# Author: Md. Jubayer Hossain
# Date: 2026-05-05

library(Seurat)
library(tidyverse)
library(decoupleR)
library(scCustom)

# 1. Load Data
merged_obj <- readRDS("data/processed/annotated.rds")

# 2. Subtype Scoring
subtype_genes <- list(
  Subtype1 = c("ARR3", "OPN1LW", "OPN1MW", "GNGT2", "GUCA1C"),
  Subtype2 = c("SOX4", "MYCN", "SOX9", "POU4F2", "NRL")
)

merged_obj <- AddModuleScore(merged_obj, features = subtype_genes, name = "RB_Subtype")
merged_obj$Subtype1_Score <- merged_obj$RB_Subtype1
merged_obj$Subtype2_Score <- merged_obj$RB_Subtype2

p1 <- FeaturePlot_scCustom(merged_obj, features = c("Subtype1_Score", "Subtype2_Score"))
ggsave("results/figures/featureplot_subtype_scores.png", p1, width = 12, height = 6, dpi = 300)

# 3. TGF-B Pathway Scoring (decoupleR)
message("Running decoupleR for TGF-B pathway...")
net <- get_progeny(organism = 'human', top = 500)
acts <- run_mlm(mat = as.matrix(merged_obj@assays$RNA$data), net = net, .source='source', .target='target', .weight='weight')

tgfb_act <- acts %>% filter(source == "TGFb") %>% select(condition, score)
merged_obj$TGFb_Activity <- tgfb_act$score[match(colnames(merged_obj), tgfb_act$condition)]

p2 <- VlnPlot(merged_obj, features = "TGFb_Activity", group.by = "disease_stage", pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white") +
  theme_classic() + ggtitle("TGF-B Activity (PROGENy)")

ggsave("results/figures/vlnplot_tgfb_activity.png", p2, width = 7, height = 6, dpi = 300)

# 4. Save Final Object
saveRDS(merged_obj, "data/processed/final_scored.rds")
message("Subtype and pathway analysis completed.")
