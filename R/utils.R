#' @title prepareMotifmatchr
#' @description Prepare objects for a \code{motifmatchr} analysis
#' @param genome character string specifying genome assembly
#' @param motifs either a character string ("jaspar", "encode", "cisbp", etc.) or an object containing PWMs
#' @return A list containing objects to be used as arguments for \code{motifmatchr}
#' @author Fabian Mueller
#' @import GenomeInfoDb TFBSTools
#' @export
prepareMotifmatchr <- function(genome, motifs) {
  res <- list()
  motifs <- tolower(motifs)
  # Get genome object
  genomeObj <- genome
  if (!inherits(genomeObj, "BSgenome")) {
    genomeObj <- getGenomeObject(genome)
  }
  spec <- GenomeInfoDb::organism(genomeObj)

  # Prepare motif PWMs
  motifL <- TFBSTools::PWMatrixList()

  if (is.character(motifs)) {
    if ("jaspar2020" %in% motifs) {
      opts <- list(species = 9606, collection = "CORE")
      mlCur <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
      if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))) {
        names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
      }
      motifL <- c(motifL, TFBSTools::toPWM(mlCur))
    }
    if (grepl("jaspar2018", motifs)) {
      opts <- list(species = spec, collection = "CORE")
      mlCur <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
      if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))) {
        names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
      }
      motifL <- c(motifL, TFBSTools::toPWM(mlCur))
    }
    if ("jaspar_vert" %in% motifs) {
      opts <- list(tax_group = "vertebrates", collection = "CORE")
      mlCur <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
      if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))) {
        names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
      }
      motifL <- c(motifL, TFBSTools::toPWM(mlCur))
    }
    if ("jaspar2016" %in% motifs) {
      motifL <- c(motifL, TFBSTools::toPWM(chromVAR::getJasparMotifs(species = spec)))
    }
    if ("homer" %in% motifs) {
      if (!requireNamespace("chromVARmotifs", quietly = TRUE)) stop("chromVARmotifs package is required but not available.")
      data("homer_pwms", package = "chromVARmotifs", envir = environment())
      motifL <- c(motifL, chromVARmotifs::homer_pwms)
    }
    if ("encode" %in% motifs) {
      if (!requireNamespace("chromVARmotifs", quietly = TRUE)) stop("chromVARmotifs package is required but not available.")
      data("encode_pwms", package = "chromVARmotifs", envir = environment())
      motifL <- c(motifL, chromVARmotifs::encode_pwms)
    }
    if ("cisbp" %in% motifs) {
      if (!requireNamespace("chromVARmotifs", quietly = TRUE)) stop("chromVARmotifs package is required but not available.")
      if (spec == "Mus musculus") {
        data("mouse_pwms_v1", package = "chromVARmotifs", envir = environment())
        motifL <- c(motifL, chromVARmotifs::mouse_pwms_v1)
      } else if (spec == "Homo sapiens") {
        data("human_pwms_v1", package = "chromVARmotifs", envir = environment())
        motifL <- c(motifL, chromVARmotifs::human_pwms_v1)
      } else {
        warning(sprintf("Could not find cisBP annotation for species: %s", spec))
      }
    }
    if ("cisbp_v2" %in% motifs) {
      if (!requireNamespace("chromVARmotifs", quietly = TRUE)) stop("chromVARmotifs package is required but not available.")
      if (spec == "Mus musculus") {
        data("mouse_pwms_v2", package = "chromVARmotifs", envir = environment())
        motifL <- c(motifL, chromVARmotifs::mouse_pwms_v2)
      } else if (spec == "Homo sapiens") {
        data("human_pwms_v2", package = "chromVARmotifs", envir = environment())
        motifL <- c(motifL, chromVARmotifs::human_pwms_v2)
      } else {
        warning(sprintf("Could not find cisBP v2 annotation for species: %s", spec))
      }
    }
    if (length(motifL) < 1) {
      stop(sprintf("No motifs were loaded. Unsupported motifs: %s", paste(motifs, collapse = ", ")))
    }
  } else if (inherits(motifs, "PWMatrixList") || inherits(motifs, "PFMatrixList")) {
    motifL <- motifs
  } else {
    stop(sprintf("Unsupported value for motifs: %s", class(motifs)))
  }

  res[["genome"]] <- genomeObj
  res[["motifs"]] <- motifL
  return(res)
}

#' @title getGenomeObject
#' @description Retrieve the appropriate \code{BSgenome} object for a given assembly string
#' @param assembly string specifying the assembly (e.g., "hg19", "hg38")
#' @param adjChrNames should chromosome names be adjusted (e.g., adding "chr" prefix)?
#' @return \code{BSgenome} object
#' @author Fabian Mueller
#' @import GenomicRanges
#' @export
getGenomeObject <- function(assembly, adjChrNames = TRUE) {
  mainREnum <- "^([1-9][0-9]?|[XYM]|MT)$"

  if (assembly %in% c("hg19")) {
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
    res <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else if (assembly %in% c("grch37", "grch37_chr")) {
    requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)
    res <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  } else if (assembly %in% c("hg38", "hg38_chr")) {
    requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
    res <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  } else if (assembly %in% c("grch38", "grch38_chr")) {
    requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)
    res <- BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens
  } else if (assembly %in% c("mm9")) {
    requireNamespace("BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE)
    res <- BSgenome.Mmusculus.UCSC.mm9::Mmusculus
  } else if (assembly %in% c("mm10")) {
    requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)
    res <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
  } else {
    stop(sprintf("Unknown assembly: %s", assembly))
  }

  if (adjChrNames) {
    prep <- grepl(mainREnum, GenomicRanges::seqnames(res))
    GenomicRanges::seqnames(res)[prep] <- paste0("chr", GenomicRanges::seqnames(res)[prep])
    if ("chrMT" %in% GenomicRanges::seqnames(res) && !"chrM" %in% GenomicRanges::seqnames(res)) {
      GenomicRanges::seqnames(res)[GenomicRanges::seqnames(res) == "chrMT"] <- "chrM"
    }
  }

  return(res)
}
#' @title findTFBindSites
#' @description Find TF binding sites from a genome and motif list (parallelized across chromosomes)
#' @param genome BSgenome object
#' @param motifs PWMatrixList or PFMatrixList of motifs
#' @param BPPARAM Parallel backend (default: automatically selected by bpparam())
#' @return A GRangesList of binding sites for each motif
#' @import GenomicRanges motifmatchr Biostrings BiocParallel
#' @export
findTFBindSites <- function(genome, motifs, BPPARAM = BiocParallel::bpparam()) {
  if (!inherits(genome, "BSgenome")) {
    stop("genome must be a BSgenome object.")
  }
  if (!inherits(motifs, "PWMatrixList") && !inherits(motifs, "PFMatrixList")) {
    stop("motifs must be a PWMatrixList or PFMatrixList.")
  }

  # Filter for main chromosomes
  all_seqnames <- GenomeInfoDb::seqnames(genome)
  seqNames <- grep("^chr[0-9XY]+$", all_seqnames, value = TRUE)

  tf_binding_list <- BiocParallel::bplapply(names(motifs), function(mo) {
    binding_sites_motif <- list()

    for (chr in seqNames) {
      chr_seq <- Biostrings::getSeq(genome, chr)
      motif_ix <- motifmatchr::matchMotifs(motifs[mo], chr_seq, out = "positions")
      hits <- unlist(motif_ix[[mo]])
      if (length(hits) == 0) next

      strand <- mcols(hits)$strand
      start_pos <- start(hits)
      motif_width <- unique(width(hits))
      if (length(motif_width) != 1) {
        warning(sprintf("Multiple motif widths found for motif %s; using first width.", mo))
        motif_width <- motif_width[1]
      }

      curr_motif <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = start_pos, width = motif_width),
        strand = strand
      )

      mcols(curr_motif)$score <- mcols(hits)$score
      curr_motif <- resize(curr_motif, width = motif_width + 400, fix = "center")

      # Store results for this chromosome
      if (is.null(binding_sites_motif[[chr]])) {
        binding_sites_motif[[chr]] <- curr_motif
      } else {
        binding_sites_motif[[chr]] <- c(binding_sites_motif[[chr]], curr_motif)
      }
    }

    # Combine results across chromosomes for this motif
    Reduce(function(x, y) c(x, y), binding_sites_motif)
  }, BPPARAM = BPPARAM)

  # Convert list to named list by motif
  names(tf_binding_list) <- names(motifs)

  tf_binding_list <- GenomicRanges::GRangesList(tf_binding_list)
  return(tf_binding_list)
}

#' @title computeGCgenome_helper
#' @description  compute GC content for a given chromosome
#' @param genome The genome object to use
#' @param chr The chromosome name
#' @param chr_len The length of the chromosome
#' @return A GRanges object with GC content values
#' @keywords internal
#' @import GenomicRanges Biostrings
computeGCgenome_helper <- function(genome, chr, chr_len) {
  qgr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = seq(1, chr_len[chr], 30), width = 30)
  )
  qgr <- qgr[-length(qgr)]
  dna_seq <- getSeq(genome, qgr)
  nucfreqs <- letterFrequency(dna_seq, c("A", "C", "G", "T"))
  gc_tmp <- na.omit(rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs))
  gc_bin <- seq(0, 1, length.out = 6)
  gcbin <- findInterval(gc_tmp, gc_bin, rightmost.closed = TRUE)
  values(qgr) <- DataFrame(GC_bias = gc_tmp, GC_bin = gcbin)
  return(qgr)
}
