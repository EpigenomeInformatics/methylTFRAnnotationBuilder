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
  library(BSgenome.Hsapiens.UCSC.hg38)
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

# Restrict the genome-wide GC table to CpG positions. methylTFR reads
# that table only through findOverlaps() against methylation calls, so
# windows that overlap no CpG are never used. For CpG methylomes this is
# lossless and cuts the table by about two orders of magnitude.
cpg <- cpgSites(BSgenome.Hsapiens.UCSC.hg38)

# Create annotations
build_annotations(annotations = "JASPAR2020",
                  pkg.base.dir = pkg.base.dir,
                  chunk_size = 2,
                  genome = BSgenome.Hsapiens.UCSC.hg38,
                  cores = 24,
                  enhancer = NULL,
                  gc_sites = cpg,      # omit for the full genome-wide table
                  keep_score = FALSE)  # methylTFR never reads the score
```

### Memory

Binding-site discovery parallelises over chromosomes and matches every
motif in one pass per chromosome, so peak memory scales with the number
of workers times one chromosome, not one chromosome per motif. Setting
`cores` above the number of chromosomes gains nothing.

The genome-wide GC scan runs in tiles, so the intermediate
nucleotide-frequency matrix stays at tile size rather than reaching
several gigabytes per chromosome. Lower `tile_size` if memory is tight.

```bash
cd methylTFRAnnotationHg38

# Install the package
R CMD INSTALL .
```