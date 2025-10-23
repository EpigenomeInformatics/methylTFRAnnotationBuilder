# methylTFRAnnotationBuilder


<!-- badges: start -->
[![GitHub issues](https://img.shields.io/github/issues/EpigenomeInformatics/methylTFRAnnotationBuilder)](https://github.com/EpigenomeInformatics/methylTFRAnnotationCreator/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/EpigenomeInformatics/methylTFRAnnotationBuilder)](https://github.com/EpigenomeInformatics/methylTFRAnnotationCreator/pulls)
<!-- badges: end -->

`methylTFRAnnotationBuilder` is an R package for creating methylTFR annotation packages.

## Installation instructions

The latest version from [GitHub](https://github.com/EpigenomeInformatics/methylTFRAnnotationCreator) with:


```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("EpigenomeInformatics/methylTFRAnnotationBuilder")
```

## Example usage

```r
# Load necessary packages
suppressPackageStartupMessages({
  library(TFBSTools)
  library(motifmatchr)
  library(Biostrings)
  library(data.table)
  library(BiocParallel)
  library(dplyr)
  library(stringr)
  library(parallel)
  library(logger)
  library(methylTFRAnnotationBuilder)})

# Set the path to package
pkg.base.dir <- getwd() 

# Create the package scaffold
createMethylTFRPackageScaffold("Hg38",dest = pkg.base.dir, motifSets = c("JASPAR2020"))

# Update the package path
pkg.base.dir <- paste0(pkg.base.dir,"methylTFRAnnotationHg38")

# Create annotations
build_annotations("JASPAR2020", pkg.base.dir, chunk_size = 2, "Hg38", cores = 10, enhancer = NULL)
```

```bash
cd methylTFRAnnotationHg38

# Install the package
R CMD INSTALL .
```