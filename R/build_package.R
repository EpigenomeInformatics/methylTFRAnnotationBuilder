#' @title createMethylTFRPackageScaffold
#' @description Creates a scaffold folder structure for a methylTFR annotation package.
#' @param assembly The genome assembly, e.g., "hg38".
#' @param dest     Destination directory where the package should be created.
#' @param motifSets a character vector of motif sets to be included in the package.
#' @return Invisibly, \code{TRUE} if the package directory and its \code{DESCRIPTION} file were successfully created;
#'         \code{FALSE} otherwise.
#' @author methylTFR Authors
#' Inspired by \code{createPackageScaffold} from the \code{RnBeadsAnnotationCreator} package.
#' @examples
#' createMethylTFRPackageScaffold("hg38")
#' @export
createMethylTFRPackageScaffold <- function(assembly, dest = getwd(), motifSets = c("JASPAR2020")) {
  pkg.name <- paste0("methylTFRAnnotation", assembly)
  desc <- c(
    Package = pkg.name,
    Title = paste("methylTFR Annotations for", assembly),
    Description = "methylTFR annotation package",
    Author = "methylTFRAnnotationBuilder",
    Date = format(Sys.Date(), format = "%Y-%m-%d"),
    License = "Artistic-2.0",
    Version = "0.1"
  )

  pkg.base.dir <- file.path(dest, pkg.name)
  if (dir.exists(pkg.base.dir)) {
    stop("Package directory already exists")
  } else {
    dir.create(pkg.base.dir)
  }

  ## Create the folder structure
  for (dname in c("R", "inst", "man", "temp")) {
    if (!dir.create(file.path(pkg.base.dir, dname), showWarnings = FALSE, recursive = TRUE)) {
      return(invisible(FALSE))
    }
  }

  # Create the extdata folder within the inst directory
  if (!dir.create(file.path(pkg.base.dir, "inst", "extdata"), showWarnings = FALSE, recursive = TRUE)) {
    return(invisible(FALSE))
  }

  ## Create DESCRIPTION file
  desc.lines <- paste(names(desc), desc, sep = ": ")
  writeLines(desc.lines, file.path(pkg.base.dir, "DESCRIPTION"))

  ## Create NAMESPACE file
  writeLines(c(""), file.path(pkg.base.dir, "NAMESPACE"))


  ## Create R script files for your functions
  functions <- list(
    getTFbindsites = "Retrieve Transcription Factor Binding sites (TFBS) annotation stored in this package",
    getGCfreq = "Load the GC frequency table stored in this package",
    getGenomeGC = "Load the genome-wide GC calculation annotation stored in this package"
  )


  writeLines(desc.lines, file.path(pkg.base.dir, "DESCRIPTION"))

  ## Create NAMESPACE file
  writeLines(c(""), file.path(pkg.base.dir, "NAMESPACE"))

  ## Create R script files for your functions
  functions <- list(
    getTFbindsites = "Retrieve Transcription Factor Binding sites (TFBS) annotation stored in this package",
    getGCfreq = "Load the GC frequency table stored in this package",
    getGenomeGC = "Load the genome-wide GC calculation annotation stored in this package"
  )

  for (func_name in names(functions)) {
    func_description <- functions[[func_name]]

    # Define function-specific logic
    if (func_name == "getTFbindsites" || func_name == "getGCfreq") {
      r_script <- paste0(
        "#' @title ", func_name, "\n",
        "#' @description ", func_description, "\n",
        "#' @param motifSet The motif set to use, default is \"", motifSets[1], "\", other options are ", paste0("\"", motifSets[-1], "\"", collapse = ", "), "\n",
        "#' @return \\code{GRangesList} object with score\n",
        "#' @export \n",
        func_name, " <- function(motifSet = \"", motifSets[1], "\"){\n",
        "  if(!tolower(motifSet) %in% tolower(motifSets)){\n",
        "    stop(\"Invalid motif set, please use one of the following: \", paste(motifSets, collapse = \", \"))\n",
        "  }\n",
        "  motifSet <- tolower(motifSet)\n",
        "  res <- readRDS(system.file(\"extdata\", paste0(motifSet, \"_", ifelse(func_name == "getTFbindsites", "tf_bindsites", "motif_gcfreq"), "\", \".rds\"), package=.PKG_NAME))\n",
        "  return(res)\n",
        "}\n"
      )
    } else if (func_name == "getGenomeGC") {
      r_script <- paste0(
        "#' @title ", func_name, "\n",
        "#' @description ", func_description, "\n",
        "#' @param assembly The genome assembly to use, no default value is set\n",
        "#' @return \\code{GRangesList} object with score\n",
        "#' @export \n",
        func_name, " <- function(assembly){\n",
        "  assembly <- tolower(assembly)\n",
        "  res <- readRDS(system.file(\"extdata\", paste0(\"genomewide_GC_\", assembly, \".rds\"), package=.PKG_NAME))\n",
        "  return(res)\n",
        "}\n"
      )
    } else {
      # Default case for other functions
      r_script <- paste0(
        "#' @title ", func_name, "\n",
        "#' @description ", func_description, "\n",
        "#' @return \\code{GRangesList} object with score\n",
        "#' @export \n",
        func_name, " <- function(){\n",
        "  res <- readRDS(system.file(\"extdata\", paste0(\"", tolower(func_name), "\", \".rds\"), package=.PKG_NAME))\n",
        "  return(res)\n",
        "}\n"
      )
    }

    # Write the script to the appropriate file
    r_file <- paste0(pkg.base.dir, "/R/", func_name, ".R")
    writeLines(r_script, r_file)
  }

  # Generate Roxygen documentation using roxygen2::roxygenise
  roxygen2::roxygenise(pkg.base.dir)

  ## Create NEWS file
  fname <- system.file(paste0("extdata/NEWS.", assembly), package = "methylTFRAnnotationCreator")
  if (file.exists(fname)) {
    txt <- scan(fname, "", sep = "\n", na.strings = character(), quiet = TRUE, blank.lines.skip = FALSE)
  } else {
    txt <- paste0("methylTFRAnnotations", assembly, " ", desc["Version"])
    txt <- c(txt, paste(rep("=", nchar(txt)), collapse = ""))
    txt <- c(txt, "", paste0("* Initial release of methylTFRAnnotations", assembly, "."))
  }
  cat(txt, file = paste0(pkg.base.dir, "/NEWS.md"), sep = "\n")

  invisible(TRUE)
}
