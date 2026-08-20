#' @title taxonomyId
#' @description Map an organism name to its NCBI taxonomy identifier.
#' @details JASPAR filters matrices by taxonomy identifier, so a genome
#' has to be translated before it can be used to select motifs. An
#' unrecognised organism returns \code{NA}, which callers treat as
#' "cannot select by species" rather than as an error, since the
#' vertebrate collection is usually still usable.
#' @param organism Organism name, as returned by
#' \code{GenomeInfoDb::organism}.
#' @return A single integer taxonomy identifier, or \code{NA_integer_}.
#' @keywords internal
taxonomyId <- function(organism) {
  ids <- c(
    "Homo sapiens" = 9606L,
    "Mus musculus" = 10090L,
    "Rattus norvegicus" = 10116L,
    "Danio rerio" = 7955L,
    "Gallus gallus" = 9031L,
    "Drosophila melanogaster" = 7227L,
    "Caenorhabditis elegans" = 6239L,
    "Saccharomyces cerevisiae" = 4932L,
    "Arabidopsis thaliana" = 3702L
  )
  if (is.null(organism) || !organism %in% names(ids)) {
    return(NA_integer_)
  }
  unname(ids[organism])
}


#' @title standardChrs
#' @description The primary assembled chromosomes of a genome.
#' @details Earlier versions took the first 24 sequence names, which is
#' the human chromosome count. On mm10, which has 19 autosomes, that
#' reaches past chrY into unplaced scaffolds; on any genome with a
#' different sequence ordering it selects an arbitrary set. The
#' chromosome list is now derived from the genome itself.
#'
#' The mitochondrion is excluded by default. Its GC content and
#' methylation are both atypical, and at 16 kb it contributes nothing to
#' genome-wide quantiles while skewing nothing but itself.
#' @param genome A \code{BSgenome} object.
#' @param drop_mito Exclude the mitochondrial sequence.
#' @return A character vector of sequence names.
#' @import GenomeInfoDb
#' @export
standardChrs <- function(genome, drop_mito = TRUE) {
  chrs <- tryCatch(
    GenomeInfoDb::standardChromosomes(genome),
    error = function(e) as.character(GenomeInfoDb::seqnames(genome))
  )
  if (length(chrs) == 0) {
    chrs <- as.character(GenomeInfoDb::seqnames(genome))
  }
  if (drop_mito) {
    chrs <- grep("^(chr)?(M|MT)$", chrs, value = TRUE, invert = TRUE)
  }
  chrs
}


#' @title prepareMotifmatchr
#' @description Prepare objects for a \code{motifmatchr} analysis
#' @details The JASPAR collections are filtered by species. Earlier
#' versions hardcoded \code{species = 9606} for jaspar2020, so asking
#' for a mouse genome returned human matrices and the build completed
#' without any error. The species is now taken from the genome unless
#' one is given explicitly.
#' @param genome character string specifying genome assembly, or a
#' \code{BSgenome} object
#' @param motifs either a character string ("jaspar", "encode", "cisbp", etc.) or an object containing PWMs
#' @param species NCBI taxonomy identifier used to filter the JASPAR
#' collections. \code{NULL} derives it from \code{genome}.
#' @param jaspar_opts Optional named list passed to
#' \code{TFBSTools::getMatrixSet} in place of the constructed options.
#' Use it when species filtering is not what you want, for example
#' \code{list(tax_group = "vertebrates", collection = "CORE")} to take
#' the whole vertebrate collection rather than the matrices annotated to
#' one species. This matters for mouse, where the species-filtered set
#' is far smaller than the human one.
#' @return A list containing objects to be used as arguments for \code{motifmatchr}
#' @author Fabian Mueller
#' @import GenomeInfoDb TFBSTools
#' @export
prepareMotifmatchr <- function(genome, motifs, species = NULL,
                               jaspar_opts = NULL) {
  res <- list()
  motifs <- tolower(motifs)
  # Get genome object
  genomeObj <- genome
  if (!inherits(genomeObj, "BSgenome")) {
    genomeObj <- getGenomeObject(genome)
  }
  spec <- GenomeInfoDb::organism(genomeObj)

  if (is.null(species)) {
    species <- taxonomyId(spec)
    if (is.na(species)) {
      warning(sprintf(
        paste(
          "No taxonomy identifier is known for %s. Pass species, or",
          "jaspar_opts, to select JASPAR matrices explicitly."
        ), spec
      ))
    }
  }

  # Species filtering unless the caller overrides it wholesale.
  jasparOptions <- function() {
    if (!is.null(jaspar_opts)) {
      return(jaspar_opts)
    }
    if (is.na(species)) {
      return(list(tax_group = "vertebrates", collection = "CORE"))
    }
    list(species = species, collection = "CORE")
  }

  # Prepare motif PWMs
  motifL <- TFBSTools::PWMatrixList()

  if (is.character(motifs)) {
    if ("jaspar2020" %in% motifs) {
      opts <- jasparOptions()
      mlCur <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
      if (length(mlCur) == 0) {
        stop(sprintf(
          "JASPAR2020 returned no matrices for: %s",
          paste(names(opts), unlist(opts), sep = " = ", collapse = ", ")
        ))
      }
      if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))) {
        names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
      }
      motifL <- c(motifL, TFBSTools::toPWM(mlCur))
    }
    if ("jaspar2018" %in% motifs) {
      opts <- jasparOptions()
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
    # cisbpv2 is accepted as well as cisbp_v2. The hg38 annotation files
    # are named cisbpv2_*, so a set requested under that name has to
    # resolve here or a from-scratch build writes files the accessors
    # cannot find.
    if (any(c("cisbp_v2", "cisbpv2") %in% motifs)) {
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
#' @description Find TF binding sites for a set of motifs across a genome.
#' @details
#' Work is parallelised over chromosomes, and every motif is matched in a
#' single \code{matchMotifs} call per chromosome. The previous
#' implementation parallelised over motifs with the chromosome loop
#' inside, so each worker held a full chromosome sequence at the same
#' time and the sequence of every chromosome was re-read once per motif.
#' For a 600 motif set that was roughly 14,400 sequence loads where 24
#' suffice, with peak memory proportional to the number of workers.
#'
#' Binding sites are returned resized to \code{width + flank}. This width
#' is load-bearing: \code{methylTFR::computeDeviation} widens the stored
#' ranges by a further 130 bases and reads methylation across the result,
#' so storing bare motif positions would silently narrow every footprint.
#' @param genome BSgenome object
#' @param motifs PWMatrixList or PFMatrixList of motifs
#' @param BPPARAM Parallel backend (default: automatically selected by
#' \code{bpparam()})
#' @param keep_score Retain the motif match score as a metadata column.
#' \code{methylTFR} never reads it, and dropping it removes eight bytes
#' per binding site from the stored object.
#' @param flank Number of bases added to each motif width.
#' @param chromosomes Character vector of sequences to scan. Defaults to
#' the primary assembled chromosomes.
#' @return A GRangesList of binding sites for each motif
#' @import GenomicRanges motifmatchr Biostrings BiocParallel
#' @export
findTFBindSites <- function(genome, motifs, BPPARAM = BiocParallel::bpparam(),
    keep_score = TRUE, flank = 400, chromosomes = NULL) {
    if (!inherits(genome, "BSgenome")) {
        stop("genome must be a BSgenome object.")
    }
    if (!inherits(motifs, "PWMatrixList") && !inherits(motifs, "PFMatrixList")) {
        stop("motifs must be a PWMatrixList or PFMatrixList.")
    }
    if (!is.logical(keep_score) || length(keep_score) != 1) {
        stop("keep_score must be a single logical value.")
    }

    motif_names <- names(motifs)
    if (is.null(motif_names) || anyDuplicated(motif_names) > 0) {
        stop("motifs must have unique names.")
    }

    if (is.null(chromosomes)) {
        chromosomes <- standardChrs(genome)
    }
    if (length(chromosomes) == 0) {
        stop("No chromosomes selected; check the sequence names.")
    }

    # One worker per chromosome. Every motif is matched in the same pass,
    # so the sequence is read once rather than once per motif.
    per_chr <- BiocParallel::bplapply(chromosomes, function(chr) {
        chr_seq <- Biostrings::getSeq(genome, chr)
        motif_ix <- motifmatchr::matchMotifs(motifs, chr_seq, out = "positions")

        lapply(motif_names, function(mo) {
            hits <- unlist(motif_ix[[mo]])
            if (length(hits) == 0) {
                return(NULL)
            }
            motif_width <- unique(width(hits))
            if (length(motif_width) != 1) {
                warning(sprintf(
                    "Multiple motif widths found for motif %s; using the first.",
                    mo
                ))
                motif_width <- motif_width[1]
            }
            gr <- GenomicRanges::GRanges(
                seqnames = chr,
                ranges = IRanges::IRanges(
                    start = start(hits), width = motif_width
                ),
                strand = mcols(hits)$strand
            )
            if (keep_score) {
                mcols(gr)$score <- mcols(hits)$score
            }
            GenomicRanges::resize(
                gr,
                width = motif_width + flank, fix = "center"
            )
        })
    }, BPPARAM = BPPARAM)

    # Collect each motif across chromosomes. A motif with no hits anywhere
    # yields an empty GRanges rather than NULL, which previously produced
    # a NULL element in the GRangesList.
    tf_binding_list <- lapply(seq_along(motif_names), function(i) {
        parts <- lapply(per_chr, function(x) x[[i]])
        parts <- parts[!vapply(parts, is.null, logical(1))]
        if (length(parts) == 0) {
            empty <- GenomicRanges::GRanges()
            if (keep_score) {
                mcols(empty)$score <- numeric(0)
            }
            return(empty)
        }
        do.call(c, parts)
    })
    names(tf_binding_list) <- motif_names

    GenomicRanges::GRangesList(tf_binding_list)
}
