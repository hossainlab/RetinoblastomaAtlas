library(Seurat)
library(tidyverse)

# 1. Setup Paths and Metadata
metadata   <- read_csv("results/tables/sample_metadata.csv")
data_dir   <- "data/raw"

# Test first 2 samples
test_samples <- metadata[c(1, 11), ] # One from each dataset

for (i in 1:nrow(test_samples)) {
  sample_id <- test_samples$sample_id[[i]]
  dataset   <- test_samples$dataset[[i]]
  message("\n>>> Testing sample: ", sample_id)
  
  path <- file.path(data_dir, dataset, sample_id)
  message("Looking for data at: ", path)
  
  if (!dir.exists(path)) {
    stop("!!! Directory not found: ", path)
  }
  
  # Check if files exist
  required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  for (f in required_files) {
    if (!file.exists(file.path(path, f))) {
      stop("!!! Missing file: ", file.path(path, f))
    }
  }
  
  # Try loading
  counts <- Read10X(data.dir = path)
  obj    <- CreateSeuratObject(counts = counts, project = sample_id)
  message("Successfully loaded ", ncol(obj), " cells.")
  
  # Check metadata
  obj$dataset       <- dataset
  obj$disease_stage <- test_samples$disease_stage[[i]]
  message("Sample metadata: Dataset=", obj$dataset[1], ", Stage=", obj$disease_stage[1])
}

message("\n>>> All tests passed! Pipeline is ready.")
