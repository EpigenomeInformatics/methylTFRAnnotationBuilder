#' @title  Build Annotations
#' @description Build annotation RDS files from either motif set names or provided GRanges.
#' @param annotations Character vector of motif set names (e.g., "cisbp", "encode") OR a named list of GRanges objects
#' representing TF binding sites.
#' @param annotations_name Name for the annotation set (default: NULL), if annotations
#' is a GRangesList, this parameter is mandatory.
#' @param pkg.base.dir Package directory; files are written to its
#' \code{inst/extdata} subdirectory.
#' @param chunk_size Number of motifs processed per checkpoint file. An
#' interrupted run resumes from the last completed chunk.
#' @param cores Number of cores to use for parallel processing (default: 10)
#' @param enhancer Optional \code{GRanges} of distal regulatory regions.
#' When given, the GC frequency tables use only binding sites overlapping
#' them and are written as \code{<set>_distal_motif_gcfreq.rds}.
#' @param keep_score Retain the motif match score on each binding site.
#' \code{methylTFR} never reads it; dropping it removes eight bytes per
#' site from the stored object.
#' @param tile_size Number of windows scanned per tile when computing
#' genome-wide GC content.
#' @param step Distance between consecutive GC window starts. The
#' default of 30 matches the window width, giving non-overlapping
#' tiles, which is what the published annotations use.
#' @param bin_scope Whether GC bin boundaries are quantiles taken across
#' the genome (default) or per chromosome. Genome scope matches the
#' quantiles used for the motif GC frequency tables below.
#' @param chromosomes Character vector of sequences to use, for both
#' motif matching and the genome-wide GC scan. \code{NULL} takes the
#' primary assembled chromosomes of \code{genome}.
#' @param species NCBI taxonomy identifier used to filter the JASPAR
#' collections. \code{NULL} derives it from \code{genome}. Only consulted
#' when motif matching actually runs, that is, when no binding-site file
#' is already present.
#' @param jaspar_opts Optional named list passed to
#' \code{TFBSTools::getMatrixSet} in place of species filtering, for
#' example \code{list(tax_group = "vertebrates", collection = "CORE")}.
#' @param genome A \code{BSgenome} object.
#' @author Irem B. Gunduz
#' @import GenomicRanges Biostrings parallel motifmatchr logger BiocParallel
#' @importFrom S4Vectors metadata
#' @export
#' @return Invisibly \code{NULL}. Called for its side effect of writing
#' \code{<set>_tf_bindsites.rds}, \code{<set>_motif_gcfreq.rds} and
#' \code{genomewide_GC_<assembly>.rds} to \code{inst/extdata}.
build_annotations <- function(
    annotations, annotations_name = NULL,
    pkg.base.dir, chunk_size = 10, genome, cores = 10, enhancer = NULL,
    keep_score = TRUE, tile_size = 5e6, step = 30L,
    bin_scope = c("genome", "chromosome"), chromosomes = NULL,
    species = NULL, jaspar_opts = NULL) {
  bin_scope <- match.arg(bin_scope)
  outdir <- file.path(pkg.base.dir, "inst/extdata")
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (is.character(annotations)) {
    mode <- "motifsets"
  } else if (inherits(annotations, "GRangesList")) {
    mode <- "GRangesList"
  } else {
    stop("annotations must be either a character vector of motif set names or a GRangesList.")
  }
  if (!inherits(genome, "BSgenome")) {
    stop("genome must be a BSgenome object.")
  }
  if (!is.numeric(cores)) {
    # Set the default number of cores to 1 if cores is not numeric
    cores <- 1
  }
  if (mode == "GRangesList") {
    if (is.null(annotations_name)) {
      stop("annotations_name must be provided when annotations is a GRangesList.")
    }
  }
  if (mode == "GRangesList") {
    tf_bindsites_list <- list()
    annotations_name <- tolower(annotations_name)
    tf_bindsites_list[[annotations_name]] <- annotations
  }
  if (mode == "motifsets") {
    tf_bindsites_list <- list()
    for (set_name in tolower(annotations)) {
      tf_file <- file.path(outdir, paste0(set_name, "_tf_bindsites.rds"))
      if (!file.exists(tf_file)) {
        log_info("Building annotation for motifset: {set_name}")
        prep <- prepareMotifmatchr(genome, set_name,
          species = species, jaspar_opts = jaspar_opts
        )
        log_info("{set_name}: {length(prep$motifs)} motifs for ",
          "{GenomeInfoDb::organism(genome)}")

        tf_bindsites <- findTFBindSites(genome, prep$motifs,
          BPPARAM = BiocParallel::MulticoreParam(workers = cores),
          chromosomes = chromosomes,
          keep_score = keep_score
        )
        tf_bindsites_list[[set_name]] <- tf_bindsites


        log_info("Saving GRanges for: {set_name}")
        saveRDS(tf_bindsites, file = tf_file)
        log_info("Saved annotation for {set_name}.")
      } else {
        log_info("Loading existing annotation for motifset: {set_name}")
        tf_bindsites <- readRDS(tf_file)
        tf_bindsites_list[[set_name]] <- tf_bindsites
        log_info("Loaded annotation for {set_name}.")
      }
    }
  }
  assembly <- unique(genome(genome))
  genome_gc_path <- file.path(outdir, paste0("genomewide_GC_", assembly, ".rds"))
  if (!file.exists(genome_gc_path)) {
    # Compute the GC dist
    log_info("Computing GC dist for the genome ...")
    gc_genome <- computeGCgenome(
      genome = genome, cores = cores, tile_size = tile_size,
      step = step, bin_scope = bin_scope,
      chromosomes = chromosomes
    )

    # Save the GC dist
    log_info("Saving GC dist for the genome ...")
    saveRDS(gc_genome, file = genome_gc_path)
  } else {
    log_info("GC dist for the genome already exists. Skipping computation.")
    gc_genome <- readRDS(genome_gc_path)
  }

  # Bin boundaries for the motif GC frequency tables.
  #
  # These MUST be the boundaries the genome table's own GC_bin column
  # was assigned with. The observed side of a deviation score reaches
  # the table through addGCBintoMethylome(), which reads GC_bin; the
  # expected side is built from the motif frequency tables binned here.
  # computeGCgenome() records the boundaries it used, so use those
  # rather than recomputing them here and hoping the two agree.
  gc_bin <- S4Vectors::metadata(gc_genome)$gc_breaks
  if (is.list(gc_bin)) {
    stop(
      "This genome GC table was built with bin_scope = \"chromosome\", ",
      "which gives each chromosome its own boundaries. The motif GC ",
      "frequency tables are genome-wide, so the two sides of a ",
      "deviation score would be binned differently. Rebuild the table ",
      "with bin_scope = \"genome\"."
    )
  }
  if (is.null(gc_bin)) {
    log_warn("Genome GC table records no bin boundaries; deriving them ",
      "from the stored GC content. Correct for a genome-wide table, ",
      "which is what older files are. Rebuild it to remove the ",
      "assumption.")
    gc_bin <- gcBreaks(gc_genome$GC_bias)
  }


  for (set_name in names(tf_bindsites_list)) {
    tf_bindsites <- tf_bindsites_list[[set_name]]
    if (!is.null(enhancer)) {
      merged_file <- file.path(
        outdir, paste0(tolower(set_name), "_distal_motif_gcfreq.rds")
      )
    } else {
      merged_file <- file.path(
        outdir, paste0(tolower(set_name), "_motif_gcfreq.rds")
      )
    }

    if (!file.exists(merged_file)) {
      # The distal and unrestricted passes over the same motif set must
      # not share a chunk cache: chunks are reused if present, so an
      # interrupted pass would leak into the other. Separate directories
      # mean neither has to be cleared, so an interrupted run resumes
      # from its completed chunks instead of starting over.
      temp_dir <- file.path(
        pkg.base.dir, "temp",
        paste0(set_name, if (is.null(enhancer)) "" else "_distal")
      )
      if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
      }

      num_chunks <- ceiling(length(tf_bindsites) / chunk_size)
      motif_list <- names(tf_bindsites)

      for (chunk_idx in 1:num_chunks) {
        chunk_filename <- file.path(temp_dir, paste0("motif_gcfreq_chunk_", chunk_idx, ".rds"))
        if (!file.exists(chunk_filename)) {
          start_idx <- (chunk_idx - 1) * chunk_size + 1
          end_idx <- min(chunk_idx * chunk_size, length(motif_list))
          chunk_motifs <- motif_list[start_idx:end_idx]

          motif_gcfreq_chunk <- parallel::mclapply(chunk_motifs,
            processMotifs2Matrix,
            tf_bindsites = tf_bindsites,
            gc_bin = gc_bin, genome = genome,
            enhancer = enhancer, mc.cores = cores
          )
          names(motif_gcfreq_chunk) <- chunk_motifs
          saveRDS(motif_gcfreq_chunk, file = chunk_filename)
        }
      }

      chunk_files <- list.files(path = temp_dir, pattern = "motif_gcfreq_chunk_\\d+\\.rds", full.names = TRUE)
      merged_data <- do.call(c, lapply(chunk_files, readRDS))
      merged_data <- merged_data[names(tf_bindsites)]
      log_info("Saving merged data for {set_name} ...")
      saveRDS(merged_data, file = merged_file)
      log_success("Saved merged data for {set_name}.")

      if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE)
        log_info("Temporary directory '{temp_dir}' has been removed.")
      }
    }
  }
  log_info("All annotations built successfully.")
  return(invisible(NULL))
}
