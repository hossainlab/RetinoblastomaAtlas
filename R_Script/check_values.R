library(GEOquery)
library(tidyverse)

# GSE168434
gse168434 <- getGEO("GSE168434", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
p168434 <- pData(gse168434)
cat("--- GSE168434 invasivenes:ch1 counts ---\n")
print(table(p168434$`invasivenes:ch1`))
cat("\n--- GSE168434 titles ---\n")
print(p168434$title)

# GSE249995
gse249995 <- getGEO("GSE249995", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
p249995 <- pData(gse249995)
cat("\n--- GSE249995 titles ---\n")
print(p249995$title)
