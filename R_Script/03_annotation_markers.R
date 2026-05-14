# Script 03: Cell Type Annotation and Marker Identification
# Project: RetinoblastomaAtlas
# Author: Md. Jubayer Hossain
# Date: 2026-05-05

library(Seurat)
library(tidyverse)
library(scCustom)

# 1. Load Data
merged_obj <- readRDS("data/processed/integrated.rds")

# 2. Marker Identification
message("Identifying cluster markers...")
markers <- FindAllMarkers(merged_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(markers, "results/tables/cluster_markers.csv")

# 3. Canonical Markers
marker_list <- list(
  "Cone_Precursor" = c("ARR3", "RXRG", "THRB", "PRDM1", "CRX"),
  "Muller_Glia"    = c("GFAP", "VIM", "SOX9", "RLBP1"),
  "TAM_Microglia"  = c("IBA1", "CD68", "AIF1", "CX3CR1"),
  "Endothelium"    = c("PECAM1", "CDH5", "VWF"),
  "Pericyte"       = c("PDGFRB", "ACTA2", "RGS5"),
  "Fibroblast"     = c("DCN", "COL1A1", "LUM"),
  "Immune_T_NK"    = c("CD3D", "CD8A", "GNLY", "NKG7"),
  "B_Cell"         = c("CD79A", "MS4A1")
)

p1 <- DotPlot_scCustom(merged_obj, features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/dotplot_canonical_markers.png", p1, width = 12, height = 6, dpi = 300)

# 4. Proportion Analysis
prop_df <- merged_obj@meta.data |>
  group_by(disease_stage, RNA_snn_res_0_5) |>
  tally() |>
  mutate(pct = n / sum(n))

p2 <- ggplot(prop_df, aes(x = disease_stage, y = pct, fill = RNA_snn_res_0_5)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  labs(y = "Proportion of Cells", x = "Disease Stage", fill = "Cluster")

ggsave("results/figures/cell_type_proportions.png", p2, width = 6, height = 6, dpi = 300)

# 5. Save Annotated Object
saveRDS(merged_obj, "data/processed/annotated.rds")

