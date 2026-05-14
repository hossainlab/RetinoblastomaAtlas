library(GEOquery)
library(tidyverse)

# Fetch GSE168434
gse168434 <- getGEO("GSE168434", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
p168434 <- pData(gse168434)
cat("--- GSE168434 pData Columns ---\n")
print(colnames(p168434))
cat("\n--- GSE168434 First 2 rows ---\n")
print(head(p168434, 2))

# Fetch GSE249995
gse249995 <- getGEO("GSE249995", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
p249995 <- pData(gse249995)
cat("\n--- GSE249995 pData Columns ---\n")
print(colnames(p249995))
cat("\n--- GSE249995 First 2 rows ---\n")
print(head(p249995, 2))
