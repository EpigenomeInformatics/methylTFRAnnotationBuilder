#' @title  Build Annotations
#' @description Build annotation RDS files from either motif set names or provided GRanges.
#' @param annotations Character vector of motif set names (e.g., "cisbp", "encode") OR a named list of GRanges objects
#' @param pkg.base.dir Output directory to save RDS files (default: inst/extdata)
#' @param chunk_size Number of annotations to save per RDS (default: 10)
#' @param cores Number of cores to use for parallel processing (default: 10)
#' @param enhancer GRanges object for distal regions (default: NULL)
#' @author Irem Gunduz
#' @param genome Character genome assembly (e.g., "hg38") or BSgenome object
#' @import GenomicRanges Biostrings parallel motifmatchr logger
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
  if (mode == "motifsets") {
    for (set_name in annotations) {
      log_info("Building annotation for motifset: {set_name}")

      # Get the BSgenome object and PWM object
      prep <- prepareMotifmatchr(genome, set_name)

      # Get tf_bindsites GRanges object
      tf_bindsites <- findTFBindSites(prep$genome, prep$motifs, cores)
      genome <- prep$genome

      # Save the GRanges object
      log_info("Saving GRanges for: {set_name}")
      saveRDS(tf_bindsites, file = file.path(outdir, paste0(set_name, "_tf_bindsites.rds")))
      log_info("Saved annotation for {set_name}.")
    }
  } else if (mode == "GRangesList") { # need to fix it for GRangesList
    for (set_name in names(annotations)) {
      tf_bindsites <- annotations[[set_name]]
      log_info("Saving provided GRanges for: {set_name}")

      if (any(width(tf_bindsites) < 400)) {
        stop(sprintf(
          "Some ranges in '%s' are too short (<400bp) to compute GC frequencies. Please ensure all ranges are at least 400bp long or adjust the width and try again.",
          set_name
        ))
      }
      saveRDS(tf_bindsites, file = file.path(outdir, paste0(set_name, "_tf_bindsites.rds")))
      log_info("Saved annotation for {set_name}.")
    }

    # Compute the GC dist for the genome
    gc_dist <- calculate_gcdist(genome = genome, threads = cores)
    gc_bin <- quantile(gc_dist, probs = seq(0, 1, 1 / 5))

    for (set_name in annotations) {
      # Compute the GC frequency for the provided GRangesList
      # Create a temp directory for GC frequency calculation with set_name
      temp_dir <- file.path(pkg.base.dir, "temp", set_name)
      if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
      }

      # Define the number of chunks to process
      num_chunks <- ceiling(length(tf_bindsites) / chunk_size)
      motif_list <- names(tf_bindsites)

      # Process motifs in chunks
      for (chunk_idx in 1:num_chunks) {
        chunk_filename <- file.path(temp_dir, paste0("motif_gcfreq_chunk_", chunk_idx, ".rds"))
        if (!file.exists(chunk_filename)) {
          # Define the motif list for the current chunk
          start_idx <- (chunk_idx - 1) * chunk_size + 1
          end_idx <- min(chunk_idx * chunk_size, length(motif_list))
          chunk_motifs <- motif_list[start_idx:end_idx]

          # Process the motifs in the chunk
          motif_gcfreq_chunk <- parallel::mclapply(chunk_motifs,
            processMotifs2Matrix,
            tf_bindsites = tf_bindsites,
            gc_bin = gc_bin, genome = genome,
            enhancer = enhancer, mc.cores = cores
          )
        }
        names(motif_gcfreq_chunk) <- chunk_motifs
        saveRDS(motif_gcfreq_chunk, file = chunk_filename)
      }

      # List all chunk files in the temp directory
      chunk_files <- list.files(path = temp_dir, pattern = "motif_gcfreq_chunk_\\d+\\.rds", full.names = TRUE)

      # Combine the results from all chunks
      merged_data <- do.call(c, lapply(chunk_files, readRDS))
      merged_data <- merged_data[names(tf_bindsites)]

      if (!is.null(enhancer)) {
        merged_file <- file.path(outdir, paste0(set_name, "_distal_motif_gcfreq.rds"))
      } else {
        merged_file <- file.path(outdir, paste0(set_name, "_motif_gcfreq.rds"))
      }

      # Remove the temp directory temp_dir
      if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE)
        log_info("Temporary directory '{temp_dir}' has been removed.")
      } else {
        log_info("Temporary directory '{temp_dir}' does not exist.")
      }
    }
    log_info("All annotations built successfully.")
    return(invisible(NULL))
  }
}
