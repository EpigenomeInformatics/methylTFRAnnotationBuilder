#' @title createMethylTFRPackageScaffold
#' @description Creates a scaffold folder structure for a methylTFR annotation package.
#' @details
#' The generated accessors validate their argument against the motif sets
#' the package was actually built with. Earlier versions wrote a check
#' against \code{motifSets}, a variable that exists only while the
#' scaffold is being generated and is undefined inside the finished
#' package, so the generated accessor failed at run time. The valid
#' names are now written into the generated file as a literal.
#'
#' File names are normalised to lower case on both sides, here and in
#' \code{build_annotations}. Previously the accessor upper-cased the
#' motif set before building the file name while \code{build_annotations}
#' wrote the name as supplied, so any set whose name was not already
#' upper case resolved to a missing file and failed inside
#' \code{readRDS("")}.
#' @param assembly The genome assembly, e.g., "hg38".
#' @param dest     Destination directory where the package should be created.
#' @param motifSets a character vector of motif sets to be included in the package.
#' @param version Version string for the generated package.
#' @param assembly_id Genome assembly identifier used to name the
#' genome-wide GC file, defaulting to \code{tolower(assembly)}. These
#' are usually the same, but they must be given separately when the
#' package is named after something other than the bare assembly: the
#' GC file is named by \code{build_annotations} from the genome object
#' itself, so a package called methylTFRAnnotationHg38Test built on hg38
#' still needs to look for genomewide_GC_hg38.rds.
#' @return Invisibly, \code{TRUE} if the package directory and its \code{DESCRIPTION} file were successfully created;
#'         \code{FALSE} otherwise.
#' @author methylTFR Authors
#' Inspired by \code{createPackageScaffold} from the \code{RnBeadsAnnotationCreator} package.
#' @examples
#' \donttest{
#' createMethylTFRPackageScaffold("hg38", dest = tempdir())
#' }
#' @export
createMethylTFRPackageScaffold <- function(assembly, dest = getwd(),
    motifSets = c("JASPAR2020"), version = "0.99.0",
    assembly_id = tolower(assembly)) {
    pkg.name <- paste0("methylTFRAnnotation", assembly)
    motifSets <- tolower(motifSets)
    if (length(motifSets) == 0 || anyDuplicated(motifSets) > 0) {
        stop("motifSets must be a non-empty vector of unique names.")
    }

    desc <- c(
        Package = pkg.name,
        Title = paste("methylTFR Annotations for", assembly),
        Description = paste(
            "Precomputed transcription factor binding sites, motif GC",
            "frequency tables and a genome-wide GC distribution for use",
            "with the methylTFR package."
        ),
        Author = "methylTFRAnnotationBuilder",
        Maintainer = "Irem B. Gunduz <irembgunduz@gmail.com>",
        Date = format(Sys.Date(), format = "%Y-%m-%d"),
        License = "Artistic-2.0",
        Encoding = "UTF-8",
        Version = version,
        Depends = "R (>= 4.4.0)",
        Imports = "GenomicRanges",
        biocViews = "AnnotationData, Genome, Homo_sapiens",
        RoxygenNote = "7.3.2"
    )

    pkg.base.dir <- file.path(dest, pkg.name)
    if (dir.exists(pkg.base.dir)) {
        stop("Package directory already exists")
    }
    dir.create(pkg.base.dir)

    ## Create the folder structure
    for (dname in c("R", "inst", "man", "temp")) {
        if (!dir.create(file.path(pkg.base.dir, dname),
            showWarnings = FALSE, recursive = TRUE
        )) {
            return(invisible(FALSE))
        }
    }
    if (!dir.create(file.path(pkg.base.dir, "inst", "extdata"),
        showWarnings = FALSE, recursive = TRUE
    )) {
        return(invisible(FALSE))
    }

    ## DESCRIPTION
    desc.lines <- paste(names(desc), desc, sep = ": ")
    writeLines(desc.lines, file.path(pkg.base.dir, "DESCRIPTION"))

    ## Shared header: package name and the motif sets that were built,
    ## written as literals so the finished package does not depend on
    ## anything from the generating session.
    header <- paste0(
        '.PKG_NAME <- "', pkg.name, '"\n',
        ".MOTIF_SETS <- c(",
        paste0('"', motifSets, '"', collapse = ", "), ")\n",
        '.ASSEMBLY <- "', tolower(assembly_id), '"\n\n',
        "#' @keywords internal\n",
        ".resolve_extdata <- function(file) {\n",
        '    path <- system.file("extdata", file, package = .PKG_NAME)\n',
        '    if (!nzchar(path)) {\n',
        '        stop("Annotation file not found in ", .PKG_NAME, ": ", file)\n',
        "    }\n",
        "    path\n",
        "}\n\n",
        "#' @keywords internal\n",
        ".check_motif_set <- function(motifSet) {\n",
        "    motifSet <- tolower(motifSet)\n",
        "    if (length(motifSet) != 1 || !motifSet %in% .MOTIF_SETS) {\n",
        '        stop("Invalid motif set. Available: ",\n',
        '            paste(.MOTIF_SETS, collapse = ", "))\n',
        "    }\n",
        "    motifSet\n",
        "}\n"
    )

    accessor <- function(fun, suffix, description) {
        paste0(
            "#' @title ", fun, "\n",
            "#' @description ", description, "\n",
            "#' @param motifSet Motif set to load. One of: ",
            paste(motifSets, collapse = ", "), ".\n",
            "#' @return The stored annotation object.\n",
            "#' @export\n",
            fun, ' <- function(motifSet = "', motifSets[1], '") {\n',
            "    motifSet <- .check_motif_set(motifSet)\n",
            '    readRDS(.resolve_extdata(paste0(motifSet, "_', suffix,
            '.rds")))\n',
            "}\n"
        )
    }

    scripts <- list(
        getTFbindsites = accessor(
            "getTFbindsites", "tf_bindsites",
            paste(
                "Retrieve transcription factor binding sites stored in",
                "this package."
            )
        ),
        getGCfreq = accessor(
            "getGCfreq", "motif_gcfreq",
            paste(
                "Load the motif GC frequency table stored in this",
                "package."
            )
        ),
        getGenomeGC = paste0(
            "#' @title getGenomeGC\n",
            "#' @description Load the genome-wide GC distribution stored",
            " in this package.\n",
            "#' @param assembly Genome assembly. Defaults to the",
            " assembly this package was built for.\n",
            "#' @return A \\code{GRanges} object with GC_bias and GC_bin.\n",
            "#' @export\n",
            "getGenomeGC <- function(assembly = .ASSEMBLY) {\n",
            "    assembly <- tolower(assembly)\n",
            '    readRDS(.resolve_extdata(paste0("genomewide_GC_",\n',
            '        assembly, ".rds")))\n',
            "}\n"
        )
    )

    # Shared constants and helpers live in one file rather than being
    # repeated in each accessor.
    writeLines(header, file.path(pkg.base.dir, "R", "aaa-utils.R"))
    for (fun in names(scripts)) {
        writeLines(
            scripts[[fun]],
            file.path(pkg.base.dir, "R", paste0(fun, ".R"))
        )
    }

    ## Generate documentation and NAMESPACE
    roxygen2::roxygenise(pkg.base.dir)

    ## NEWS
    txt <- paste0("methylTFRAnnotation", assembly, " ", desc[["Version"]])
    txt <- c(txt, paste(rep("=", nchar(txt)), collapse = ""))
    txt <- c(txt, "", paste0(
        "* Initial release of methylTFRAnnotation", assembly, "."
    ))
    cat(txt, file = file.path(pkg.base.dir, "NEWS.md"), sep = "\n")

    invisible(TRUE)
}
