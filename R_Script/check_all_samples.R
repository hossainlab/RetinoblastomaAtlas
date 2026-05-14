library(GEOquery)
library(tidyverse)

# GSE168434 - handle multiple platforms
gse168434 <- getGEO("GSE168434", GSEMatrix = TRUE, getGPL = FALSE)
p168434_list <- lapply(gse168434, pData)
p168434 <- bind_rows(p168434_list)

cat("--- GSE168434 All Samples ---\n")
print(p168434 %>% select(title, geo_accession, `invasivenes:ch1`))

# GSE249995
gse249995 <- getGEO("GSE249995", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
p249995 <- pData(gse249995)
cat("\n--- GSE249995 All Samples ---\n")
print(p249995 %>% select(title, geo_accession))
