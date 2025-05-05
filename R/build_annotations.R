#' @title  Build Annotations
#' @description Build annotation RDS files from either motif set names or provided GRanges.
#' @param annotations Character vector of motif set names (e.g., "cisbp", "encode") OR a named list of GRanges objects
#' @param pkg.base.dir Output directory to save RDS files (default: inst/extdata)
#' @param chunk_size Number of annotations to save per RDS (default: 10)
#' @param cores Number of cores to use for parallel processing (default: 10)
#' @param enhancer GRanges object for distal regions (default: NULL)
#' @author Irem Gunduz
#' @param genome Character genome assembly (e.g., "hg38") or BSgenome object
#' @import GenomicRanges Biostrings parallel motifmatchr logger BiocParallel
#' @export
#' @return NULL
build_annotations <- function(annotations, pkg.base.dir, chunk_size = 10, genome, cores = 10, enhancer = NULL) {
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
  if (!inherits(genome, "BSgenome") && !is.character(genome)) {
    stop("genome must be either a BSgenome object or a character assembly string.")
  }
  if (!is.numeric(cores)) {
    # Set the default number of cores to 1 if cores is not numeric
    cores <- 1
  }
  tf_bindsites_list <- annotations
  if (mode == "motifsets") {
    for (set_name in annotations) {
      log_info("Building annotation for motifset: {set_name}")
      if (is.character(genome)) {
        genome <- tolower(genome)
      }

      prep <- prepareMotifmatchr(genome, set_name)
      assembly <- unique(GenomeInfoDb::genome(prep$genome))
      genome <- prep$genome

      tf_bindsites <- findTFBindSites(genome, prep$motifs, BPPARAM = BiocParallel::MulticoreParam(workers = cores))
      tf_bindsites_list[[set_name]] <- tf_bindsites

      tf_file <- file.path(outdir, paste0(set_name, "_tf_bindsites.rds"))
      if (!file.exists(tf_file)) {
        log_info("Saving GRanges for: {set_name}")
        saveRDS(tf_bindsites, file = tf_file)
        log_info("Saved annotation for {set_name}.")
      }
    }
  } else if (mode == "GRangesList") {
    for (set_name in names(annotations)) {
      tf_bindsites <- annotations[[set_name]]
      log_info("Saving provided GRanges for: {set_name}")

      if (any(width(tf_bindsites) < 400)) {
        stop(sprintf(
          "Some ranges in '%s' are too short (<400bp) to compute GC frequencies. Please ensure all ranges are at least 400bp long or adjust the width and try again.",
          set_name
        ))
      }

      tf_file <- file.path(outdir, paste0(set_name, "_tf_bindsites.rds"))
      if (!file.exists(tf_file)) {
        saveRDS(tf_bindsites, file = tf_file)
        log_info("Saved annotation for {set_name}.")
      }
    }
  }
  genome_gc_path <- file.path(outdir, paste0("genome_wide_GC_", assembly, ".rds"))
  if (!file.exists(genome_gc_path)) {
    # Compute the GC dist
    log_info("Computing GC dist for the genome ...")
    gc_genome <- computeGCgenome(genome = genome, cores = cores)

    # Save the GC dist
    log_info("Saving GC dist for the genome ...")
    saveRDS(gc_genome, file = genome_gc_path)
  } else {
    log_info("GC dist for the genome already exists. Skipping computation.")
  }

  # Compute the GC dist for the genome for TFBS usage
  gc_dist <- calculate_gcdist(genome = genome, threads = cores)
  gc_bin <- quantile(gc_dist, probs = seq(0, 1, 1 / 5))


  for (set_name in names(tf_bindsites_list)) {
    tf_bindsites <- tf_bindsites_list[[set_name]]
    if (!is.null(enhancer)) {
      merged_file <- file.path(outdir, paste0(set_name, "_distal_motif_gcfreq.rds"))
    } else {
      merged_file <- file.path(outdir, paste0(set_name, "_motif_gcfreq.rds"))
    }

    if (!file.exists(merged_file)) {
      temp_dir <- file.path(pkg.base.dir, "temp", set_name)
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
