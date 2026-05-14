# Script 05: Cell-Cell Communication (CellChat)
# Project: RetinoblastomaAtlas
# Author: Md. Jubayer Hossain
# Date: 2026-05-05

library(CellChat)
library(Seurat)
library(tidyverse)
library(patchwork)

# 1. Load Data
merged_obj <- readRDS("data/processed/final_scored.rds")

# 2. Split by Disease Stage
obj_list <- SplitObject(merged_obj, split.by = "disease_stage")

run_cellchat <- function(seurat_obj) {
  cellchat <- createCellChat(object = seurat_obj, group.by = "ident", assay = "RNA")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

message("Running CellChat on Intraocular...")
chat_intra <- run_cellchat(obj_list$intraocular)

message("Running CellChat on Extraocular...")
chat_extra <- run_cellchat(obj_list$extraocular)

# 3. Comparative Analysis
chat_list <- list(Intra = chat_intra, Extra = chat_extra)
chat_merged <- mergeCellChat(chat_list, add.names = names(chat_list))

p1 <- compareInteractions(chat_merged, show.legend = F, group = c(1,2))
p2 <- compareInteractions(chat_merged, show.legend = F, group = c(1,2), measure = "weight")
ggsave("results/figures/cellchat_interaction_comparison.png", p1 + p2, width = 10, height = 5, dpi = 300)

# TGF-B Pathway Comparison
pathway <- "TGFb"
p3 <- netVisual_aggregate(chat_intra, signaling = pathway, layout = "circle")
p4 <- netVisual_aggregate(chat_extra, signaling = pathway, layout = "circle")

png("results/figures/cellchat_tgfb_circle.png", width = 1200, height = 600, res = 150)
par(mfrow = c(1,2))
p3; p4
dev.off()

# 4. Save CellChat Objects
saveRDS(chat_merged, "data/processed/cellchat_merged.rds")
message("CellChat analysis completed.")
