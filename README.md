# methylTFRAnnotationBuilder

<!-- badges: start -->
[![GitHub issues](https://img.shields.io/github/issues/EpigenomeInformatics/methylTFRAnnotationBuilder)](https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/EpigenomeInformatics/methylTFRAnnotationBuilder)](https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder/pulls)
<!-- badges: end -->

`methylTFRAnnotationBuilder` builds the genome annotations used by
[methylTFR](https://github.com/EpigenomeInformatics/methylTFR) and
scaffolds them into an installable annotation package. It is the tool
behind
[methylTFRAnnotationHg38](https://github.com/EpigenomeInformatics/methylTFRAnnotationHg38)
and
[methylTFRAnnotationMm10](https://github.com/EpigenomeInformatics/methylTFRAnnotationMm10).

For a genome (`BSgenome`) and a set of motifs it produces three files:

| File | Content |
|---|---|
| `<set>_tf_bindsites.rds` | `GRangesList` of genome-wide motif matches, one element per motif, each extended by 200 bp on either side |
| `<set>_motif_gcfreq.rds` | list of 5 x n matrices: for each position along a motif's footprint, the fraction of binding sites in each genome-wide GC quintile |
| `genomewide_GC_<assembly>.rds` | `GRanges` of non-overlapping 30 nt windows with their GC fraction (`GC_bias`) and GC quintile (`GC_bin`) |

methylTFR uses the GC tables to correct TF deviation scores for sequence
composition.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("EpigenomeInformatics/methylTFRAnnotationBuilder")
```

You also need the `BSgenome` package for your assembly (for example
`BSgenome.Hsapiens.UCSC.hg38`) and the package that provides your motif
set (`JASPAR2020`, or `chromVARmotifs` for `cisbpv2`, `homer` and
`encode`).

## Example

```r
library(methylTFRAnnotationBuilder)
library(BSgenome.Hsapiens.UCSC.hg38)

genome <- BSgenome.Hsapiens.UCSC.hg38
dest <- getwd()

# 1. Create the package skeleton
createMethylTFRPackageScaffold("Hg38", dest = dest, motifSets = "jaspar2020")
pkg_dir <- file.path(dest, "methylTFRAnnotationHg38")

# 2. Compute binding sites, the genome-wide GC table and the motif GC
#    frequency tables into pkg_dir/inst/extdata
build_annotations(
    annotations = "jaspar2020",
    pkg.base.dir = pkg_dir,
    genome = genome,
    cores = 24,
    chunk_size = 15,
    keep_score = FALSE,       # methylTFR never reads the match score
    chromosomes = standardChrs(genome)
)

# Optional: GC frequency tables restricted to distal regulatory regions,
# written as jaspar2020_distal_motif_gcfreq.rds
# build_annotations("jaspar2020", pkg.base.dir = pkg_dir, genome = genome,
#                   enhancer = distal_regions)
```

Then install the result:

```bash
R CMD INSTALL methylTFRAnnotationHg38
```

Supported motif set names are `jaspar2020`, `jaspar2018`, `jaspar_vert`,
`jaspar2016`, `cisbp`, `cisbpv2`, `homer` and `encode`; see
`?prepareMotifmatchr`. You can also pass your own `PWMatrixList` to
`findTFBindSites()`, or a `GRangesList` of binding sites to
`build_annotations()`.

The scripts that produced the published annotations are in
`inst/scripts/make-data.R` of
[methylTFRAnnotationHg38](https://github.com/EpigenomeInformatics/methylTFRAnnotationHg38)
and
[methylTFRAnnotationMm10](https://github.com/EpigenomeInformatics/methylTFRAnnotationMm10).

### Memory

Binding-site discovery runs in parallel over chromosomes and matches every
motif in one pass per chromosome, so peak memory is roughly the number of
workers times one chromosome. Setting `cores` above the number of
chromosomes gains nothing.

The genome-wide GC scan runs in tiles, so the intermediate
nucleotide-frequency matrix stays at tile size. Lower `tile_size` if
memory is tight. The motif GC tables are computed in chunks of
`chunk_size` motifs; each finished chunk is saved, so an interrupted run
resumes where it stopped.

## Citation

If you use methylTFRAnnotationBuilder, please cite it together with
methylTFR:

> Gunduz IB, Mueller F (2026). methylTFRAnnotationBuilder: Build annotation packages for methylTFR. R package version 0.99.2. https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder

> Gunduz IB, Murugan SK, Mueller F (2026). methylTFR: Quantification of DNA methylation signatures in TFBS. https://github.com/EpigenomeInformatics/methylTFR

## License

MIT, see [LICENSE.md](LICENSE.md).
