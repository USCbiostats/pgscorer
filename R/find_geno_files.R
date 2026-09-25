#' Find the Genotype Files compute_prs() Would Use
#'
#' Looks in a directory the same way \code{\link{compute_prs}} does and reports
#' which genetic data files it finds, which chromosome(s) each one holds, and
#' anything likely to cause trouble — without running any scoring. Use it to
#' check a directory before a long run.
#'
#' Discovery is identical to \code{compute_prs()}: a VCF is a \code{*.vcf.gz}
#' with a \code{.tbi} index, a BinaryDosage file is any \code{*.bdose} (Format 5
#' with its \code{.bdose.bdi}, or single-file Format 4, optionally with a saved
#' \code{.bdinfo}). When \code{format} is \code{NULL} and both kinds are
#' present, BinaryDosage is used. Chromosomes are read from each file's own
#' metadata, so this reads every BinaryDosage file's information (a full parse
#' of the data file for Format 4 files without a \code{.bdinfo}, which can be
#' slow for large files — the report flags these).
#'
#' Unlike \code{compute_prs()}, a chromosome found in more than one file is
#' reported (\code{duplicated = TRUE}) rather than stopping with an error.
#'
#' @param geno_dir Directory to search. Defaults to the current working
#'   directory.
#' @param format   \code{"vcf"} or \code{"bdose"} to look only at that type;
#'   \code{NULL} (default) autodetects exactly as \code{compute_prs()} does.
#' @param verbose  Logical; if \code{TRUE} (default) print a report of what was
#'   found, including a list of potential problems.
#'
#' @return Invisibly, a data.frame with one row per file and chromosome:
#'   \describe{
#'     \item{chrom}{Chromosome label without any \code{"chr"} prefix — the form
#'       matched against the PGS files.}
#'     \item{contig}{Chromosome label as stored in the file.}
#'     \item{file}{File name.}
#'     \item{format}{\code{"vcf"} or \code{"bdose"}.}
#'     \item{version}{BinaryDosage format number (4 or 5); \code{NA} for VCF.}
#'     \item{info}{Where the file's information came from: \code{"tbi"} (VCF),
#'       \code{"bdi"} (Format 5), \code{"bdinfo"} (saved \code{.bdinfo}), or
#'       \code{"parsed"} (read from the whole data file).}
#'     \item{n_samples}{Number of samples in the file.}
#'     \item{duplicated}{\code{TRUE} if the chromosome appears in more than one
#'       file (\code{compute_prs()} would stop).}
#'     \item{path}{Full path.}
#'   }
#'   The attribute \code{"issues"} holds the problems that were found, as a
#'   character vector (empty if none).
#'
#' @examples
#' \dontrun{
#' find_geno_files("path/to/genotypes")
#' files <- find_geno_files("path/to/genotypes", format = "vcf", verbose = FALSE)
#' }
#' @export
find_geno_files <- function(geno_dir = ".", format = NULL, verbose = TRUE) {
  stopifnot(
    "tabixr is not installed" = requireNamespace("tabixr", quietly = TRUE),
    "geno_dir does not exist" = dir.exists(geno_dir),
    "format must be NULL, 'vcf', or 'bdose'" =
      is.null(format) || (is.character(format) && length(format) == 1L && format %in% c("vcf", "bdose"))
  )

  geno    <- .scan_geno_files(geno_dir, format, verbose = FALSE)
  samples <- attr(geno, "samples")
  unused  <- attr(geno, "unused")
  chosen  <- geno$format[1L]

  geno$duplicated <- geno$chrom %in% geno$chrom[duplicated(geno$chrom)]
  geno <- geno[order(match(geno$chrom, .chrom_sort(geno$chrom)), geno$file),
               c("chrom", "contig", "file", "format", "version", "info",
                 "n_samples", "duplicated", "path")]
  rownames(geno) <- NULL

  # ---- Things compute_prs() would stop on or be slow at -----------------------
  issues <- character(0)
  for (ch in unique(geno$chrom[geno$duplicated]))
    issues <- c(issues, sprintf("Chromosome %s is in more than one file (%s); compute_prs() would stop.",
                                ch, paste(geno$file[geno$chrom == ch], collapse = ", ")))
  slow <- unique(geno$file[geno$info %in% "parsed"])
  if (length(slow) > 0L)
    issues <- c(issues, sprintf(
      "%d Format 4 file(s) have no .bdinfo and are parsed in full each time (slow for large files); create one with saveRDS(getbdinfo(file), \"<name>.bdinfo\"): %s",
      length(slow), paste(slow, collapse = ", ")))
  sets <- unique(lapply(samples, function(s) sort(s)))
  if (length(sets) > 1L)
    issues <- c(issues, "Sample IDs differ between files; scores need every file to contain the same samples.")
  attr(geno, "issues") <- issues

  if (verbose) {
    label  <- c(vcf = "VCF", bdose = "BinaryDosage")
    n_file <- length(unique(geno$file))
    cat(sprintf("Genotype files in %s\n", geno_dir))
    cat(sprintf("  Format used : %s%s\n", label[[chosen]],
                if (is.null(format)) " (autodetected)" else " (requested)"))
    cat(sprintf("  Files       : %d, covering %d chromosome(s)\n", n_file, length(unique(geno$chrom))))
    if (length(unused$paths) > 0L)
      cat(sprintf("  Not used    : %d %s file(s) also present%s\n", length(unused$paths),
                  label[[unused$format]],
                  if (is.null(format)) "" else sprintf(" (format = '%s' was requested)", format)))
    cat("\n")
    shown <- geno[, c("chrom", "contig", "file", "version", "info", "n_samples", "duplicated")]
    shown$version <- ifelse(is.na(shown$version), "-", as.character(shown$version))
    print(shown, row.names = FALSE)
    cat("\n")
    if (length(issues) == 0L) {
      cat("No problems found.\n")
    } else {
      cat("Potential problems:\n")
      cat(paste0("  - ", issues, "\n"), sep = "")
    }
  }

  invisible(geno)
}
