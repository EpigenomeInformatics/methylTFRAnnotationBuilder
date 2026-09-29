#' @title processMotifs2Matrix
#' @description Tabulate, for every offset along a motif's footprint
#' window, how its binding sites distribute across the five GC bins.
#' @details
#' Each binding site is widened by 130 bases (65 on either side) and a
#' 30 nt window is slid along it one base at a time. Every window is
#' assigned to one of the five GC bins defined by \code{gc_bin}, and for
#' each window position the fraction of binding sites falling into each
#' bin is recorded. Windows containing \code{N} produce no GC value and
#' are dropped rather than binned.
#'
#' Sites are processed in batches: each batch is unlisted into a single
#' sequence and all of its windows are counted in one
#' \code{letterFrequency} call over a \code{Views} object, then
#' accumulated with \code{tabulate}. Windows cannot straddle a site
#' boundary because each site contributes only offsets 1 to width minus
#' 29.
#' @param motif Name of the motif to process (an element name of
#' \code{tf_bindsites}).
#' @param gc_bin The six GC bin boundaries, as returned by
#' \code{\link{gcBreaks}}.
#' @param genome A \code{BSgenome} object.
#' @param tf_bindsites A \code{GRangesList} of binding sites, as returned
#' by \code{\link{findTFBindSites}}.
#' @param enhancer Optional \code{GRanges}. When given, only binding
#' sites overlapping these regions are used (for example distal
#' regulatory regions).
#' @param batch_size Number of binding sites held in memory at once.
#' Peak memory is roughly \code{batch_size} times the window count times
#' sixteen bytes.
#' @import GenomicRanges Biostrings
#' @return A numeric matrix with five rows (GC bins, lowest first) and
#' one column per window position; each column sums to one.
#' @export
processMotifs2Matrix <- function(motif, gc_bin, genome, tf_bindsites,
                                 enhancer = NULL, batch_size = 20000L) {
  tfbs <- tf_bindsites[[motif]]
  tfbs <- resize(tfbs, width(tfbs) + 130, fix = "center")
  win <- 30L
  W <- as.integer(width(tfbs)[1])
  if (is.na(W) || W < win) {
    warning(sprintf("Motif %s has no usable binding sites.", motif))
    return(NULL)
  }
  nw <- W - win + 1L

  if (!is.null(enhancer)) {
    tfbs <- IRanges::subsetByOverlaps(tfbs, enhancer, ignore.strand = TRUE)
  }
  n <- length(tfbs)
  logger::log_info(paste0(
    "Processing ", motif, ": ", format(n, big.mark = ","), " sites"
  ))
  counts <- matrix(0, nrow = 5L, ncol = nw)
  if (n == 0L) {
    warning(sprintf("Motif %s has no binding sites in the given regions.",
      motif))
    return(sweep(counts, 2, colSums(counts), FUN = "/"))
  }

  for (from in seq.int(1L, n, by = batch_size)) {
    idx <- seq.int(from, min(from + batch_size - 1L, n))
    dna <- Biostrings::getSeq(genome, tfbs[idx])
    m <- length(dna)
    u <- unlist(dna)
    st <- rep.int((seq_len(m) - 1L) * W, rep.int(nw, m)) + rep(seq_len(nw), m)
    lf <- Biostrings::letterFrequency(
      Biostrings::Views(u, start = st, width = win),
      c("A", "C", "G", "T")
    )
    b <- findInterval(
      (lf[, "C"] + lf[, "G"]) / rowSums(lf), gc_bin, rightmost.closed = TRUE
    )
    p <- rep(seq_len(nw), m)
    keep <- !is.na(b) & b >= 1L
    counts <- counts + matrix(
      tabulate((p[keep] - 1L) * 5L + b[keep], nbins = 5L * nw),
      nrow = 5L, ncol = nw
    )
    rm(dna, u, lf, b, p, keep)
  }
  logger::log_info(paste("Normalizing matrix...", motif))
  sweep(counts, 2, colSums(counts), FUN = "/")
}

#' @title computeGCgenome
#' @description Compute GC content in 30 nt windows and assign each
#' window to one of five GC bins.
#' @details
#' Windows are non-overlapping by default: one row per 30 bases of
#' genome, about 103 million for hg38. The scan is chunked so the
#' intermediate nucleotide-frequency matrix stays bounded, and only the
#' requested windows are counted. Counting every offset and discarding
#' the unwanted rows costs thirty times the work and thirty times the
#' memory for an identical result.
#'
#' Bin boundaries are quantiles of the observed GC distribution.
#' \code{bin_scope} controls whether those quantiles are taken once
#' across the genome, the default, or per chromosome. Genome scope is
#' the correct pairing, because \code{build_annotations} derives the
#' motif GC frequency tables from genome-wide quantiles; under
#' per-chromosome bins the observed and expected sides of a deviation
#' score are binned on different scales.
#' @param genome A \code{BSgenome} object.
#' @param cores Number of parallel workers.
#' @param tile_size Number of windows scanned per tile.
#' @param step Distance between consecutive window starts, in bases.
#' The default of 30 equals the window width, giving abutting
#' non-overlapping tiles: one row per 30 bases of genome. This is what
#' the published methylTFR annotations use, and it is what
#' \code{addGCBintoMethylome} assumes, since it needs each methylation
#' call to fall in exactly one window. Setting it to 1 gives a sliding
#' scan at every offset, which produces roughly thirty times as many
#' windows -- 2.9 billion for hg38 against 103 million -- and is not
#' comparable with an annotation built at the default.
#' @param bin_scope Either "genome" (default) or "chromosome". Genome
#' scope matches the quantiles \code{build_annotations} uses for the
#' motif GC frequency tables, so the observed and expected sides of a
#' deviation score are binned on the same scale. \code{build_annotations}
#' refuses a table built with "chromosome" scope.
#' @param chromosomes Character vector of sequences to scan. Defaults to
#' the primary assembled chromosomes of \code{genome}, excluding the
#' mitochondrion (see \code{\link{standardChrs}}).
#' @return A \code{GRanges} object with \code{GC_bias} and \code{GC_bin}
#' metadata columns.
#' @import GenomicRanges Biostrings
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom stats quantile
#' @importFrom S4Vectors metadata metadata<-
#' @export
computeGCgenome <- function(genome, cores = 1, tile_size = 5e6,
    step = 30L, bin_scope = c("genome", "chromosome"),
    chromosomes = NULL) {
    bin_scope <- match.arg(bin_scope)
    if (!inherits(genome, "BSgenome")) {
        stop("genome must be a BSgenome object.")
    }
    chr_len <- GenomeInfoDb::seqlengths(genome)
    if (is.null(chromosomes)) {
        chromosomes <- standardChrs(genome)
    }
    missing_chr <- setdiff(chromosomes, names(chr_len))
    if (length(missing_chr) > 0) {
        stop(
            "These sequences are not in the genome: ",
            paste(missing_chr, collapse = ", ")
        )
    }

    # One worker per chromosome. Each holds at most one tile of sequence
    # plus its own share of the output, rather than a whole chromosome.
    param <- BiocParallel::MulticoreParam(workers = max(1L, as.integer(cores)))

    # Workers return plain vectors, not GRanges. Building one GRanges per
    # tile and concatenating them costs a full copy per concatenation and
    # carries the seqnames Rle, ranges and mcols machinery through every
    # intermediate. At genome scale that dominates both the runtime and
    # the peak memory. Vectors concatenate once, cheaply, and a single
    # GRanges is constructed at the end.
    per_chr <- BiocParallel::bplapply(chromosomes, function(chr) {
        computeGCgenome_helper(
            genome = genome, chr = chr, chr_len = chr_len,
            tile_size = tile_size, step = step
        )
    }, BPPARAM = param)
    names(per_chr) <- chromosomes

    n <- vapply(per_chr, function(x) length(x$start), integer(1))
    per_chr <- per_chr[n > 0]
    n <- n[n > 0]
    if (length(per_chr) == 0) {
        stop("No GC windows were produced; check the chromosome names.")
    }

    starts <- unlist(lapply(per_chr, `[[`, "start"), use.names = FALSE)
    gc <- unlist(lapply(per_chr, `[[`, "gc"), use.names = FALSE)
    chr_of <- names(per_chr)
    rm(per_chr)
    invisible(gc())

    if (bin_scope == "chromosome") {
        bins <- integer(length(gc))
        offset <- 0L
        breaks <- list()
        for (i in seq_along(n)) {
            idx <- seq.int(offset + 1L, offset + n[i])
            b <- gcBreaks(gc[idx])
            bins[idx] <- findInterval(gc[idx], b, rightmost.closed = TRUE)
            breaks[[chr_of[i]]] <- b
            offset <- offset + n[i]
        }
    } else {
        breaks <- gcBreaks(gc)
        bins <- findInterval(gc, breaks, rightmost.closed = TRUE)
    }

    res <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(factor(chr_of, levels = chr_of), n),
        ranges = IRanges::IRanges(start = starts, width = 30L)
    )
    rm(starts)
    invisible(gc())
    res$GC_bias <- gc
    res$GC_bin <- bins

    # The boundaries travel with the object, so build_annotations() bins
    # the motif windows with the same ones rather than recomputing
    # quantiles from the stored GC content. Recomputing is correct only
    # when the stored windows are the whole genome; keeping the numbers
    # attached removes that assumption entirely.
    S4Vectors::metadata(res)$gc_breaks <- breaks
    S4Vectors::metadata(res)$bin_scope <- bin_scope
    S4Vectors::metadata(res)$step <- as.integer(step)
    res
}


#' @title gcBreaks
#' @description The five-bin GC boundaries of a reference distribution.
#' @details Kept as its own function so the genome table and the motif
#' GC frequency tables cannot drift apart: both take their boundaries
#' from here, applied to the same reference.
#' @param reference Numeric vector of GC fractions describing the
#' genome, not merely the windows that were retained.
#' @return A numeric vector of six quantiles.
#' @importFrom stats quantile
#' @export
gcBreaks <- function(reference) {
    stats::quantile(reference, probs = seq(0, 1, 1 / 5), na.rm = TRUE)
}


#' @title computeGCgenome_helper
#' @description Compute GC content in 30 nt sliding windows across one
#' chromosome, tile by tile.
#' @param genome A \code{BSgenome} object.
#' @param chr Name of the sequence to scan.
#' @param chr_len Named vector of sequence lengths.
#' @param tile_size Number of windows scanned per tile.
#' @param step Distance between consecutive window starts.
#' @return A list with an integer \code{start} vector and a numeric
#' \code{gc} vector, one element per window. Vectors rather than a
#' \code{GRanges}, because the caller concatenates across chromosomes
#' and building the object once at the end is both faster and far
#' cheaper in memory.
#' @import GenomicRanges Biostrings
#' @keywords internal
computeGCgenome_helper <- function(genome, chr, chr_len, tile_size = 5e6,
    step = 30L) {
    win <- 30L
    len <- as.integer(chr_len[chr])
    empty <- list(start = integer(0), gc = numeric(0))
    if (len < win) {
        return(empty)
    }

    # Genomic start of every window, once. With step == win these are
    # abutting tiles; with step == 1 it is a sliding scan.
    all_starts <- seq.int(1L, len - win + 1L, by = as.integer(step))
    nwin <- length(all_starts)

    per_tile <- max(1L, as.integer(tile_size))
    tile_from <- seq.int(1L, nwin, by = per_tile)

    starts <- vector("list", length(tile_from))
    gcs <- vector("list", length(tile_from))

    for (i in seq_along(tile_from)) {
        a <- tile_from[i]
        b <- min(a + per_tile - 1L, nwin)
        gs <- all_starts[a:b]

        seq_from <- gs[1L]
        seq_to <- gs[length(gs)] + win - 1L
        tile_seq <- Biostrings::getSeq(
            genome, chr,
            start = seq_from, end = min(seq_to, len)
        )

        # Count over exactly the windows wanted rather than over every
        # offset and discarding the rest.
        v <- Biostrings::Views(
            tile_seq,
            start = gs - seq_from + 1L, width = win
        )
        nucfreqs <- Biostrings::letterFrequency(
            v, letters = c("A", "C", "G", "T")
        )

        starts[[i]] <- gs
        # Denominator is the window width, not the count of called
        # bases, so a window inside an assembly gap scores 0 rather
        # than NaN and every window is retained. This is what the
        # published annotations contain: their window count equals the
        # genome divided by the step, gaps included.
        gcs[[i]] <- (nucfreqs[, "C"] + nucfreqs[, "G"]) / win
        rm(nucfreqs, v, tile_seq)
    }

    list(
        start = unlist(starts, use.names = FALSE),
        gc = unlist(gcs, use.names = FALSE)
    )
}

