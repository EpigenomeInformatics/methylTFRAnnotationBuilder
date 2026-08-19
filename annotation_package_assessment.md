# methylTFR annotation: assessment and plan

Review of `methylTFRAnnotationBuilder` and the generated
`methylTFRAnnotationHg38`, with a route to a Bioconductor-acceptable
structure that leaves `methylTFR` itself unchanged.

---

## 1. The size problem

`methylTFRAnnotationHg38/inst/extdata` currently holds:

| File | Size |
|---|---|
| `cisbpv2_tf_bindsites.rds` | 2.35 GB |
| `JASPAR2020_tf_bindsites.rds` | 1.75 GB |
| `altius_tf_bindsites.rds` | 999 MB |
| `genomewide_GC_hg38.rds` | 327 MB |
| `cisbpv2_motif_gcfreq.rds` | 15.6 MB |
| `JASPAR2020_motif_gcfreq.rds` | 11.9 MB |
| `cisbpv2_distal_motif_gcfreq.rds` | 9.5 MB |
| `JASPAR2020_DISTAL_motif_gcfreq.rds` | 7.2 MB |
| `altius_motif_gcfreq.rds` | 5.2 MB |
| `distal_regions.RDS` | 0.4 MB |
| **Total** | **~5.4 GB** |

A Bioconductor software package is capped at 5 MB. Annotation packages
are allowed to be larger but are expected to be modest, and this is
roughly a thousand times over. Shipping these files inside a package is
not an option in any form. That is the real blocker, not a formatting
issue.

---

## 2. Why the builder crashes

Two independent memory problems. Either alone is enough to kill a large
run; the README example (`cores = 30`) triggers both.

### 2.1 `findTFBindSites()` parallelises over the wrong axis

```r
tf_binding_list <- BiocParallel::bplapply(names(motifs), function(mo) {
  for (chr in seqNames) {
    chr_seq <- Biostrings::getSeq(genome, chr)
    motif_ix <- motifmatchr::matchMotifs(motifs[mo], chr_seq, out = "positions")
```

The outer loop is over motifs and the chromosome loop is inside. Two
consequences:

**Redundant sequence loading.** Every motif re-reads all 24
chromosomes. For a 600-motif set that is ~14,400 `getSeq()` calls where
24 would do.

**Peak memory scales with worker count.** Each of the `cores` workers
holds a full chromosome `DNAString` at the same time. `chr1` alone is
~250 MB before `matchMotifs` allocates anything. At `cores = 30` that is
7.5 GB of sequence, plus motifmatchr's internals, plus the copy-on-write
footprint of the forked `genome` object.

**Fix:** invert the loops. Parallelise over chromosomes, and inside each
chromosome call `matchMotifs()` once with the entire `PWMatrixList`
rather than once per motif. motifmatchr is built for that. This gives
24 sequence loads instead of 14,400, caps workers at 24, and each worker
holds exactly one chromosome.

### 2.2 `computeGCgenome_helper()` allocates whole-chromosome matrices

```r
nucfreqs <- Biostrings::letterFrequencyInSlidingView(
    seq_set, view.width = 30, letters = c("A", "C", "G", "T"))
```

For `chr1` this is a 248,956,393 x 4 integer matrix, about **4 GB**.
`gc_tmp` adds another ~2 GB, and the resulting `GRanges` carries a start,
a width, a `GC_bias` double and a `GC_bin` integer for every one of those
windows. Multiply by the number of workers.

`do.call(c, bplapply(...))` then holds all 24 chromosome results in
memory simultaneously before concatenating them.

**Fix:** process each chromosome in tiles (10 Mb blocks, say) so the
sliding-view matrix never exceeds a few hundred MB, and accumulate
results incrementally rather than concatenating at the end.

### 2.3 A smaller bug in the same function

```r
param <- bpparam()
param$workers <- cores
```

Assigning `$workers` on a `BiocParallelParam` does not reliably change
the worker count. Use `BiocParallel::MulticoreParam(workers = cores)`
explicitly, as `build_annotations()` already does elsewhere.

---

## 3. Most of the genome-wide GC table is never used

This is the cheapest large win available, and it needs no change to
`methylTFR`.

`methylTFR` touches `gc_dist` in exactly one place,
`addGCBintoMethylome()`:

```r
hits <- findOverlaps(msites, gcdist, ignore.strand = ignoreStrand)
gcmap <- data.table(mscore = msites[hits@from]$score,
                    gcbin  = gcdist[hits@to]$GC_bin)
```

Two facts follow.

**Only windows that overlap a methylation call are ever read.** The
current table stores every 30-bp window at 1-bp step across the genome,
about 2.9 billion windows. A CpG methylome queries roughly 28 million
positions. Restricting the shipped table to CpG contexts is a reduction
of roughly two orders of magnitude and is lossless for CpG analysis.

**Only `GC_bin` is read, never `GC_bias`.** The `GC_bias` column is used
once, inside `build_annotations()`, to derive the quantile breaks:

```r
gc_dist <- gc_genome$GC_bias
gc_bin  <- quantile(gc_dist, probs = seq(0, 1, 1/5))
```

That is builder-side. Dropping `GC_bias` from the distributed object
halves what remains; store the five quantile breaks alongside it so the
builder can still recompute bins without the full vector.

Caveat: restricting to CpG contexts assumes CpG methylation. Non-CpG
(CHG/CHH) analysis would need the fuller table, so this should be a
build-time option rather than a hard-coded assumption.

---

## 4. The TFBS objects can be trimmed too

`findTFBindSites()` attaches a `score` column to every binding site:

```r
mcols(curr_motif)$score <- mcols(hits)$score
```

`methylTFR` never reads it. `computeDeviation()` uses `tf_bindsites`
only for ranges and the derived `mid_point`; the `score` it does read
belongs to `msites`, not to the binding sites. Dropping the column saves
8 bytes per site across hundreds of millions of sites.

Beyond that, the site count itself is driven by the `matchMotifs`
p-value threshold, which is left at the default. Raising it would cut
the object proportionally, at the cost of dropping weak sites.

One thing that must **not** change: sites are stored resized to
`motif_width + 400`, and `computeDeviation()` then applies
`resize(tfbs, width(tfbs)[1] + 130)`. The stored width is load-bearing.
Storing bare motif positions would silently narrow every footprint
window.

---

## 5. Bugs in the generated accessors

### 5.1 The motif-set validation never fires

Shipped `getGCfreq.R` and `getTFbindsites.R` both contain:

```r
if (!toupper(motifSet) %in% toupper(motifSet)) {
```

The argument is compared against itself, so the condition is always
`FALSE` and no input is ever rejected. The generator in
`build_package.R` writes `tolower(motifSets)` (the vector of valid
names), so the installed package does not match what the builder
produces — it has been edited by hand at some point, or generated by an
older version.

### 5.2 Case handling makes most motif sets unreachable

The accessors uppercase the argument and then build a filename:

```r
motifSet <- toupper(motifSet)
readRDS(system.file("extdata", paste0(motifSet, "_motif_gcfreq.rds"), ...))
```

Against the files actually on disk:

| Call | Looks for | On disk | Works |
|---|---|---|---|
| `getGCfreq("JASPAR2020")` | `JASPAR2020_motif_gcfreq.rds` | same | yes |
| `getGCfreq("jaspar2020_distal")` | `JASPAR2020_DISTAL_motif_gcfreq.rds` | same | yes |
| `getGCfreq("cisbpv2")` | `CISBPV2_motif_gcfreq.rds` | `cisbpv2_...` | **no** |
| `getGCfreq("altius")` | `ALTIUS_motif_gcfreq.rds` | `altius_...` | **no** |

`system.file()` returns `""` for a missing file and `readRDS("")` then
fails with an unhelpful connection error. This is almost certainly why
the analysis scripts bypass the accessor:

```r
gcfreqs <- readRDS("/icbb/.../jaspar2020_distal_motif_gcfreq_vs2.rds")
```

`build_annotations()` writes `paste0(set_name, "_motif_gcfreq.rds")`
using the name as given, while the accessor forces case. The two need to
agree; normalising both to lower case is the simpler fix.

### 5.3 `getGenomeGC()` has no default

The signature is `getGenomeGC(assembly)` with no default, but the
analysis scripts call `getGenomeGC()` with no arguments, which errors.
The `methylTFR` vignette uses `getGenomeGC("hg38")`. Since an annotation
package is built for exactly one assembly, the assembly it was built for
is the obvious default.

### 5.4 `build_package.R` duplicates a block

Lines 45–71 write `DESCRIPTION` and define the `functions` list twice.
Harmless, but it means edits have to be made in two places.

---

## 6. Route to a Bioconductor-acceptable structure

The requirement is that the returned objects keep their current classes
so `methylTFR` needs no changes. That is satisfied by every option below,
because only the retrieval mechanism changes, not the payload.

### Preferred: AnnotationHub

`methylTFRAnnotationHg38` becomes a small package (well under 5 MB)
containing only accessor code and an `inst/extdata/metadata.csv`
describing each resource. The `.rds` files are hosted in Bioconductor's
S3 bucket and fetched on first use, then cached locally by
`BiocFileCache`.

```r
getTFbindsites <- function(motifSet = "jaspar2020") {
    ah <- AnnotationHub::AnnotationHub()
    ah[[.resource_id(motifSet, "tf_bindsites")]]
}
```

Return value: the same `GRangesList`. `methylTFR` is untouched.

What this needs: `biocViews: AnnotationHub, AnnotationData`, a
`DataSource` field, and a `metadata.csv` that validates against
`AnnotationHubData::makeAnnotationHubMetadata()`. Bioconductor hosts the
files at no cost.

The realistic obstacle is size. A single 2.35 GB resource will attract
questions, so sections 3 and 4 should be done first. If the GC table
drops to CpG contexts and the TFBS objects lose the unused `score`
column, the whole set becomes far easier to justify.

### Fallback: Zenodo plus BiocFileCache

If Bioconductor declines to host the larger resources, deposit them on
Zenodo and have the accessors download and cache through
`BiocFileCache`. This is fully acceptable in a Bioconductor package, is
less negotiation, and has a side benefit for the manuscript: Zenodo
issues a DOI, which is exactly what the Availability section needs.

### Not viable

Shipping the data inside the package, in any arrangement. The AWS
tarball route currently described in the `methylTFR` vignette is also
not acceptable to Bioconductor, since there is no versioning, no cache
integration and no guarantee of persistence.

---

## 7. Suggested order of work

1. **Fix the two memory bombs** (2.1, 2.2). Needed regardless of where
   the data ends up, and until this is done the builder cannot be re-run
   to produce anything smaller.
2. **Shrink the data** (3, 4). Restrict the GC table to CpG contexts and
   drop `GC_bias`; drop the unused `score` column from the TFBS objects.
   Re-measure before deciding anything about hosting.
3. **Fix the accessor bugs** (5). Small, self-contained, and they are
   currently forcing the analysis scripts to work around the package.
4. **Choose the hosting route** (6) once the real sizes are known.
5. `methylTFR` needs no changes at any point in this.

Steps 1 and 3 are worth doing even if the annotation package is deferred,
because the builder is a submitted-package dependency in its own right
and the accessor bugs are visible to any user who tries a motif set
other than JASPAR2020.
