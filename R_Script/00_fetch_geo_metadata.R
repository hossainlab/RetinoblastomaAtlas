# Script 00: Robust GEO Metadata Extraction and Harmonization
# Project: RetinoblastomaAtlas
# Author: Gemini CLI (Interactive Agent)
# Date: 2026-05-07

library(GEOquery)
library(tidyverse)

# Create directories
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("--- Fetching GSE168434 ---")
gse168434 <- getGEO("GSE168434", GSEMatrix = TRUE, getGPL = FALSE)
# Combine multiple platforms
p168434 <- bind_rows(lapply(gse168434, pData)) %>%
  as_tibble() %>%
  mutate(across(everything(), as.character))

message("--- Fetching GSE249995 ---")
gse249995 <- getGEO("GSE249995", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
p249995 <- pData(gse249995) %>%
  as_tibble() %>%
  mutate(across(everything(), as.character))

# --- Process GSE168434 ---
meta_168434 <- p168434 %>%
  transmute(
    original_title = title,
    geo_accession = geo_accession,
    dataset = "GSE168434",
    patient_id = str_extract(title, "RB\\d+"),
    # Standardize disease stage
    disease_stage = case_when(
      str_detect(`invasivenes:ch1`, "Non") ~ "intraocular",
      `invasivenes:ch1` == "Invasive" ~ "extraocular",
      TRUE ~ NA_character_
    ),
    tissue = `tissue:ch1`,
    sex = `Sex:ch1`,
    age = `age:ch1`,
    rb1_mutation = `blood mut (rb1):ch1`,
    # Technical metadata
    instrument = instrument_model,
    library_strategy = library_strategy,
    replicate = ifelse(str_detect(title, "rep1"), "rep1", 
                       ifelse(str_detect(title, "rep2"), "rep2", "rep1")),
    # Match directory names like GSM5139852_RB01_rep1
    sample_id = paste0(geo_accession, "_", str_replace(title, "_mRNA", ""))
  )

# --- Process GSE249995 ---
meta_249995 <- p249995 %>%
  transmute(
    original_title = title,
    geo_accession = geo_accession,
    dataset = "GSE249995",
    patient_id = paste0("P", str_extract(title, "\\d+")),
    disease_stage = ifelse(str_detect(title, "Intraocular"), "intraocular", "extraocular"),
    tissue = ifelse(str_detect(title, "Intraocular"), "primary_tumour", "optic_nerve"),
    sex = NA_character_,
    age = NA_character_,
    rb1_mutation = NA_character_,
    instrument = instrument_model,
    library_strategy = library_strategy,
    replicate = "rep1",
    # We moved these to folders S1_in1, etc.
    sample_id = case_when(
      str_detect(title, "Intraocular RB 1") ~ "S1_in1",
      str_detect(title, "Intraocular RB 2") ~ "S2_in2",
      str_detect(title, "Extraocular RB 1") ~ "S3_ex1",
      str_detect(title, "Extraocular RB 2") ~ "S4_ex2"
    )
  )

# --- Harmonize and Save ---
sample_metadata <- bind_rows(meta_168434, meta_249995)

# Final cleanup
sample_metadata <- sample_metadata %>%
  mutate(across(where(is.character), str_trim)) %>%
  # Add metadata from the original file if possible (like n_cells)
  mutate(n_cells = case_when(
    str_detect(sample_id, "RB01_rep1") ~ 6990,
    str_detect(sample_id, "RB01_rep2") ~ 6826,
    str_detect(sample_id, "RB02_rep1") ~ 4185,
    str_detect(sample_id, "RB02_rep2") ~ 2140,
    str_detect(sample_id, "RB03_rep1") ~ 7596,
    str_detect(sample_id, "RB03_rep2") ~ 7638,
    str_detect(sample_id, "RB04") ~ 14093,
    str_detect(sample_id, "RB05") ~ 13016,
    str_detect(sample_id, "RB06") ~ 14407,
    str_detect(sample_id, "RB07") ~ 14681,
    sample_id == "S1_in1" ~ 10453,
    sample_id == "S2_in2" ~ 15386,
    sample_id == "S3_ex1" ~ 11099,
    sample_id == "S4_ex2" ~ 13370,
    TRUE ~ NA_real_
  ))

write_csv(sample_metadata, "data/metadata/sample_metadata_full.csv")
write_csv(sample_metadata, "results/tables/sample_metadata.csv")

message("--- Harmonized metadata saved to results/tables/sample_metadata.csv ---")
print(sample_metadata %>% select(sample_id, dataset, disease_stage, patient_id))
