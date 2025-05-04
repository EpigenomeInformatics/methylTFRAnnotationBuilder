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
    chr_len = chr_len,
    genome = genome,
    mc.cores = threads
  )
  return(unlist(gc_dist))
}

#' @title get_gcdist
#' @description  get the GC distribution for a given chromosome
#' @param chr The chromosome name
#' @param chr_len The length of the chromosome
#' @param genome The genome object to use
#' @import GenomicRanges Biostrings
#' @return GC content values for the chromosome
#' @export
get_gcdist <- function(chr, chr_len, genome) {
  qgr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = seq(1, chr_len[chr], 30), width = 30)
  )
  # last interval might be out of chromosome length
  qgr <- qgr[-length(qgr)]
  dna_seq <- Biostrings::getSeq(genome, qgr)
  # calculate GC content
  nucfreqs <- Biostrings::letterFrequency(dna_seq, c("A", "C", "G", "T"))
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
  mat[cbind(bins[, 2], bins[, 1])] <- 1
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
  tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")
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

#' @title computeGCgenome
#' @description  compute GC content for the entire genome
#' @param genome The genome object to use
#' @param cores The number of cores to use for parallel processing
#' @return A GRanges object with GC content values for the entire genome
#' @import GenomicRanges Biostrings
#' @export
#' @importFrom BiocParallel bplapply
computeGCgenome <- function(genome, cores = 1) {
  chr_len <- seqlengths(genome)
  chr_names <- names(chr_len)[1:24]

  # Automatically select the parallel backend
  param <- bpparam()
  param$workers <- cores # Set the number of cores/workers

  # Use bplapply for parallel processing
  t_qgr <- do.call(c, bplapply(chr_names, compute_gc_genome, chr_len, BPPARAM = param))
  return(t_qgr)
}
