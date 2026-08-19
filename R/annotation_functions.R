#' @title calculate_gcdist
#' @description  calculating GC distribution and GC frequency table
#' @param genome The genome object to use, no default value is set
#' @param threads The number of threads to use for parallel processing, default is 1
#' @return GC content values for each chromosome
#' @import GenomicRanges Biostrings parallel
#' @export
calculate_gcdist <- function(genome, threads = 1) {
  chr_len <- seqlengths(genome)
  chr_names <- names(chr_len)[1:24]
  gc_dist <- parallel::mclapply(chr_names, get_gcdist,
    genome = genome,
    mc.cores = threads
  )
  return(unlist(gc_dist))
}

#' @title get_gcdist
#' @description  get the GC distribution for a given chromosome
#' @param chr The chromosome name
#' @param genome The genome object to use
#' @import GenomicRanges Biostrings
#' @return GC content values for the chromosome
#' @export
get_gcdist <- function(chr, genome) {
  seq_set <- Biostrings::getSeq(genome, chr)
  # calculate GC content
  nucfreqs <- Biostrings::letterFrequencyInSlidingView(
    seq_set,
    view.width = 30,
    letters = c("A", "C", "G", "T")
  )
  gc_tmp <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  gc_tmp <- na.omit(gc_tmp)

  return(gc_tmp)
}

#' @title compute_gc
#' @description  compute the GC content for a given DNA sequence
#' @param X The DNA sequence to compute GC content for
#' @return GC content value
#' @import Biostrings
#' @export
compute_gc <- function(X) {
  x <- DNAString(as.character(X))
  # center around
  nucfreqs <- Biostrings::letterFrequencyInSlidingView(x,
    view.width = 30,
    letters = c("A", "C", "G", "T")
  )
  gc <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  return(gc)
}

#' @title convert_to_bins
#' @description  convert GC content values to bins
#' @param x The GC content values to convert
#' @param gc_bin The GC bins to use for conversion
#' @return A matrix of bin indices
#' @export
convert_to_bins <- function(x, gc_bin) {
  bin <- cbind(1:length(x), findInterval(x, gc_bin, rightmost.closed = TRUE))
  return(bin)
}

#' @title convert_to_matrix
#' @description convert bin indices to a matrix
#' @param bins The bin indices to convert
#' @return A matrix of bin indices
#' @export
convert_to_matrix <- function(bins) {
  mat <- matrix(0, 5, dim(bins)[1])
  bins <- na.omit(bins)
  mat[base::cbind(bins[, 2], bins[, 1])] <- 1
  return(mat)
}


#' @title processMotifs2Matrix
#' @description  process motifs to matrix
#' @param motif The motif to process
#' @param gc_bin The GC bins to use for conversion
#' @param genome The genome object to use
#' @param tf_bindsites The TF binding sites to use
#' @param enhancer The enhancer regions to use
#' @import GenomicRanges Biostrings
#' @return A matrix of GC content values
#' @export
processMotifs2Matrix <- function(motif, gc_bin, genome, tf_bindsites, enhancer = NULL) {
  tfbs <- tf_bindsites[[motif]]
  tfbs <- resize(tfbs, width(tfbs) + 130, fix = "center")
  if (!is.null(enhancer)) {
    tfbs <- subsetByOverlaps(tfbs, enhancer, ignore.strand = T)
  }
  dna_seq <- Biostrings::getSeq(genome, tfbs)
  logger::log_info(paste("Processing compute gc .. ", motif))
  motif_gc <- lapply(dna_seq, compute_gc)
  logger::log_info(paste("Processing convert to bins .. ", motif))
  gcbins <- lapply(motif_gc, convert_to_bins, gc_bin)
  logger::log_info(paste("Processing convert to matrix ... ", motif))
  gcmat <- lapply(gcbins, convert_to_matrix)
  m_gcfreq <- Reduce("+", gcmat, accumulate = FALSE)
  logger::log_info(paste("Normalizing matrix...", motif))
  normalized_matrix <- sweep(as.matrix(m_gcfreq), 2, colSums(as.matrix(m_gcfreq)), FUN = "/")
  return(normalized_matrix)
}

#' @title cpgSites
#' @description Locate every CpG dinucleotide in a genome, as a
#' \code{GRanges} of single-base positions marking the cytosine.
#' @details Supplying the result as the \code{sites} argument of
#' \code{computeGCgenome} restricts the GC table to positions a CpG
#' methylome can actually query. \code{methylTFR} reads the table only
#' through \code{findOverlaps(msites, gcdist)}, so windows that overlap
#' no methylation call are never used. For CpG data this is lossless and
#' reduces the table by roughly two orders of magnitude.
#' @param genome A \code{BSgenome} object.
#' @param chromosomes Character vector of sequences to scan. Defaults to
#' the first 24 sequences, matching the rest of the builder.
#' @return A \code{GRanges} of width-1 ranges, one per CpG.
#' @import GenomicRanges Biostrings
#' @export
cpgSites <- function(genome, chromosomes = NULL) {
    if (!inherits(genome, "BSgenome")) {
        stop("genome must be a BSgenome object.")
    }
    if (is.null(chromosomes)) {
        chromosomes <- names(GenomeInfoDb::seqlengths(genome))[seq_len(24)]
    }
    res <- lapply(chromosomes, function(chr) {
        hits <- Biostrings::matchPattern(
            "CG", Biostrings::getSeq(genome, chr)
        )
        GenomicRanges::GRanges(
            seqnames = chr,
            ranges = IRanges::IRanges(start = start(hits), width = 1L)
        )
    })
    do.call(c, res)
}


#' @title computeGCgenome
#' @description Compute GC content in 30 nt sliding windows and assign
#' each window to one of five GC bins.
#' @details
#' The sliding-window scan is performed in tiles so that the intermediate
#' nucleotide-frequency matrix never exceeds the size of one tile. The
#' previous implementation scanned a whole chromosome in one call, which
#' allocated an n-by-4 integer matrix of roughly four gigabytes for
#' chr1, in every parallel worker simultaneously.
#'
#' Bin boundaries are quantiles of the observed GC distribution.
#' \code{bin_scope} controls whether those quantiles are taken once
#' across the genome, the default, or per chromosome. Genome scope is
#' the correct pairing, because \code{build_annotations} derives the
#' motif GC frequency tables from genome-wide quantiles; under
#' per-chromosome bins the observed and expected sides of a deviation
#' score are binned on different scales. Earlier versions binned per
#' chromosome, so annotations built before this change are not
#' comparable with annotations built after it.
#' @param genome A \code{BSgenome} object.
#' @param cores Number of parallel workers.
#' @param tile_size Number of windows scanned per tile.
#' @param sites Optional \code{GRanges} restricting the output to
#' windows starting at these positions, for example the result of
#' \code{cpgSites}. \code{NULL} keeps every position, which produces a
#' very large object.
#' @param bin_scope Either "genome" (default) or "chromosome". Genome
#' scope matches the quantiles \code{build_annotations} uses for the
#' motif GC frequency tables, so the observed and expected sides of a
#' deviation score are binned on the same scale. "chromosome"
#' reproduces the behaviour of earlier versions, where the two sides
#' were binned differently.
#' @param chromosomes Character vector of sequences to scan. Defaults to
#' the first 24 sequences.
#' @return A \code{GRanges} object with \code{GC_bias} and \code{GC_bin}
#' metadata columns.
#' @import GenomicRanges Biostrings
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom stats quantile
#' @export
computeGCgenome <- function(genome, cores = 1, tile_size = 5e6,
    sites = NULL, bin_scope = c("genome", "chromosome"),
    chromosomes = NULL) {
    bin_scope <- match.arg(bin_scope)
    if (!inherits(genome, "BSgenome")) {
        stop("genome must be a BSgenome object.")
    }
    if (!is.null(sites) && !inherits(sites, "GRanges")) {
        stop("sites must be NULL or a GRanges object.")
    }
    chr_len <- GenomeInfoDb::seqlengths(genome)
    if (is.null(chromosomes)) {
        chromosomes <- names(chr_len)[seq_len(24)]
    }

    # One worker per chromosome. Each holds at most one tile of sequence
    # plus its own share of the output, rather than a whole chromosome.
    param <- BiocParallel::MulticoreParam(workers = max(1L, as.integer(cores)))

    per_chr <- BiocParallel::bplapply(chromosomes, function(chr) {
        computeGCgenome_helper(
            genome = genome, chr = chr, chr_len = chr_len,
            tile_size = tile_size, sites = sites
        )
    }, BPPARAM = param)

    per_chr <- per_chr[vapply(per_chr, length, integer(1)) > 0]
    if (length(per_chr) == 0) {
        stop("No GC windows were produced; check the chromosome names.")
    }

    if (bin_scope == "chromosome") {
        per_chr <- lapply(per_chr, function(gr) {
            gr$GC_bin <- assign_gc_bins(gr$GC_bias)
            gr
        })
        return(do.call(c, per_chr))
    }

    res <- do.call(c, per_chr)
    res$GC_bin <- assign_gc_bins(res$GC_bias)
    return(res)
}


#' @title assign_gc_bins
#' @description Assign GC values to five bins delimited by their own
#' quantiles.
#' @param gc Numeric vector of GC fractions.
#' @return An integer vector of bin indices.
#' @importFrom stats quantile
#' @keywords internal
assign_gc_bins <- function(gc) {
    breaks <- stats::quantile(gc, probs = seq(0, 1, 1 / 5), na.rm = TRUE)
    findInterval(gc, breaks, rightmost.closed = TRUE)
}


#' @title computeGCgenome_helper
#' @description Compute GC content in 30 nt sliding windows across one
#' chromosome, tile by tile.
#' @param genome A \code{BSgenome} object.
#' @param chr Name of the sequence to scan.
#' @param chr_len Named vector of sequence lengths.
#' @param tile_size Number of windows scanned per tile.
#' @param sites Optional \code{GRanges} restricting the output.
#' @return A \code{GRanges} object with a \code{GC_bias} metadata column.
#' @import GenomicRanges Biostrings
#' @keywords internal
computeGCgenome_helper <- function(genome, chr, chr_len, tile_size = 5e6,
    sites = NULL) {
    win <- 30L
    nwin <- max(0L, as.integer(chr_len[chr]) - win + 1L)
    if (nwin == 0L) {
        return(GenomicRanges::GRanges())
    }

    keep_starts <- NULL
    if (!is.null(sites)) {
        on_chr <- sites[as.character(GenomicRanges::seqnames(sites)) == chr]
        if (length(on_chr) == 0) {
            return(GenomicRanges::GRanges())
        }
        keep_starts <- sort(unique(GenomicRanges::start(on_chr)))
        keep_starts <- keep_starts[keep_starts >= 1L & keep_starts <= nwin]
        if (length(keep_starts) == 0) {
            return(GenomicRanges::GRanges())
        }
    }

    tile_size <- max(1L, as.integer(tile_size))
    tile_from <- seq.int(1L, nwin, by = tile_size)

    pieces <- lapply(tile_from, function(from) {
        to <- min(from + tile_size - 1L, nwin)
        # A window starting at `to` ends at `to + win - 1`, so the
        # sequence slice has to run that far.
        tile_seq <- Biostrings::getSeq(
            genome, chr,
            start = from, end = min(to + win - 1L, as.integer(chr_len[chr]))
        )
        if (length(tile_seq) < win) {
            return(NULL)
        }
        nucfreqs <- Biostrings::letterFrequencyInSlidingView(
            tile_seq,
            view.width = win, letters = c("A", "C", "G", "T")
        )
        totals <- rowSums(nucfreqs)
        starts <- seq.int(from, length.out = nrow(nucfreqs))
        valid <- totals > 0
        if (!is.null(keep_starts)) {
            valid <- valid & starts %in% keep_starts
        }
        if (!any(valid)) {
            return(NULL)
        }
        gc <- rowSums(nucfreqs[valid, 2:3, drop = FALSE]) / totals[valid]
        gr <- GenomicRanges::GRanges(
            seqnames = chr,
            ranges = IRanges::IRanges(start = starts[valid], width = win)
        )
        gr$GC_bias <- gc
        gr
    })

    pieces <- pieces[!vapply(pieces, is.null, logical(1))]
    if (length(pieces) == 0) {
        return(GenomicRanges::GRanges())
    }
    do.call(c, pieces)
}
