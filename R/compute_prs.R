#' Compute Polygenic Risk Scores from PGS Catalog Models
#'
#' Reads one or more PGS Catalog scoring files, determines which chromosomes
#' they need, and scores each model against the genotype files found in
#' \code{geno_dir} — either BGZF/Tabix VCF or BinaryDosage (Format 4 or 5).
#' Positions are unioned across models on each chromosome so a shared position
#' is only queried once.
#'
#' Genotype files are found by extension, so file names are otherwise
#' arbitrary: every \code{*.vcf.gz} with a \code{.tbi} index is a VCF, and every
#' \code{*.bdose} is a BinaryDosage file (Format 5, which also needs its
#' \code{.bdose.bdi}, or the single-file Format 4; see below for a faster
#' way to use Format 4 files). The chromosome(s) in
#' each file are read from the file's own metadata (VCF: the contig names in
#' the Tabix index; BinaryDosage: the chromosome column of the header or
#' \code{.bdi}). A leading
#' \code{"chr"} is ignored when matching against the PGS files' chromosome
#' labels. If two files (of the format in use) contain the same chromosome,
#' an error is raised. If a model references a chromosome for which no
#' genotype file exists, every SNP on that chromosome is reported as unmatched
#' (with a warning) rather than causing an error.
#'
#' The genotype files need not contain exactly the same samples (for example
#' subjects dropped from some chromosomes for poor imputation quality). Only
#' samples present in every file that is actually queried are scored, since a
#' chromosome's contribution can only be added for samples every file has; the
#' others are listed in \code{excluded_samples} and, if \code{verbose}, counted
#' in the log. If no sample is common to all the files, an error is raised.
#'
#' Format 4 (and earlier) BinaryDosage files have no companion information
#' file, so \code{BinaryDosage::getbdinfo()} must parse the whole data file to
#' get their chromosomes, samples, and SNP offsets. To avoid that, save the
#' result next to the data file with \code{.bdose} replaced by \code{.bdinfo},
#' e.g. \code{saveRDS(getbdinfo("chr1.bdose"), "chr1.bdinfo")}; if that file
#' exists it is used instead. Format 5 files are unaffected (they use their
#' \code{.bdi}, and any \code{.bdinfo} beside them is ignored). The path stored
#' inside a \code{.bdinfo} is always replaced with the actual location of the
#' data file, so the files may be moved or copied after it was created. A
#' \code{.bdinfo} that cannot be read, or whose SNP offsets are inconsistent
#' with the data file, is ignored with a warning and the data file is parsed
#' instead. Regenerate it whenever the \code{.bdose} file is rewritten; a
#' \code{.bdinfo} left over from an older version of the data file cannot be
#' detected beyond those consistency checks.
#'
#' On Windows, VCF batches are queried in-process; on other platforms each VCF
#' batch runs in its own \code{Rscript} subprocess to bound peak memory.
#' BinaryDosage files are always queried in-process (random-access reads, no
#' large sequential scan to isolate).
#'
#' @param geno_dir   Directory containing the genotype files: BGZF VCFs with a
#'   \code{.vcf.gz.tbi} index and/or BinaryDosage \code{.bdose} files (Format 4,
#'   or Format 5 with its \code{.bdose.bdi}). Defaults to the current working
#'   directory.
#' @param format     \code{"vcf"} or \code{"bdose"} to force the genotype file
#'   type; \code{NULL} (default) autodetects from what is present in
#'   \code{geno_dir}. If both types are present and \code{format} is
#'   \code{NULL}, BinaryDosage files are used and a message is printed.
#' @param pgs_files  Character vector of paths to PGS Catalog scoring files
#'   (\code{.txt.gz}). If \code{NULL} (default), all files matching
#'   \code{^PGS.*\\.txt\\.gz$} in \code{pgs_dir} are used.
#' @param pgs_dir    Directory searched for PGS files when \code{pgs_files} is
#'   \code{NULL}. Defaults to the current working directory.
#' @param batch_size Integer; positions per VCF query batch (ignored for
#'   BinaryDosage input). Default 10000.
#' @param output_dir Directory in which to write one \code{.rds} file per
#'   model, named \code{pgs_<model>_prs.rds}. Set to \code{NULL} to skip
#'   saving. Defaults to the current working directory.
#' @param verbose    Logical; if \code{TRUE} (default) print progress and
#'   summary messages.
#'
#' @return A named list (invisibly), one element per model, each itself a list
#'   with:
#'   \describe{
#'     \item{prs}{Named numeric vector of PRS values, one per scored sample
#'       (those present in every genotype file that was queried).}
#'     \item{excluded_samples}{Character vector of sample IDs left out because
#'       they were missing from at least one queried file (the same for every
#'       model); empty when all files have the same samples.}
#'     \item{unmatched_by_chr}{Named integer vector, one entry per chromosome
#'       referenced by any model or found in \code{geno_dir}: the number of
#'       this model's SNPs on that chromosome with no matching genotype
#'       record (including every SNP on a chromosome with no genotype file).}
#'     \item{unmatched_rsIDs}{Named list, one character vector of unmatched
#'       rsIDs per chromosome (same chromosome set as \code{unmatched_by_chr}).}
#'     \item{excluded_chr_counts}{Named integer vector, one entry per
#'       chromosome: SNPs found in the genotype data but excluded because the
#'       effect allele matched neither REF nor ALT.}
#'   }
#'
#' @export
compute_prs <- function(geno_dir   = ".",
                        format     = NULL,
                        pgs_files  = NULL,
                        pgs_dir    = ".",
                        batch_size = 10000L,
                        output_dir = ".",
                        verbose    = TRUE) {

  stopifnot(
    "tabixr is not installed"          = requireNamespace("tabixr",     quietly = TRUE),
    "data.table is not installed"      = requireNamespace("data.table", quietly = TRUE),
    "geno_dir does not exist"          = dir.exists(geno_dir),
    "batch_size must be a positive integer" = is.numeric(batch_size) && batch_size >= 1L,
    "format must be NULL, 'vcf', or 'bdose'" = is.null(format) || (is.character(format) && format %in% c("vcf", "bdose"))
  )
  batch_size <- as.integer(batch_size)

  # ---- Discover PGS files -----------------------------------------------------
  if (is.null(pgs_files)) {
    pgs_files <- sort(list.files(pgs_dir, pattern = "^PGS.*\\.txt\\.gz$", full.names = TRUE))
    if (length(pgs_files) == 0L)
      stop(sprintf("No PGS model files (^PGS.*\\.txt\\.gz$) found in: %s", pgs_dir))
  } else {
    missing_pgs <- pgs_files[!file.exists(pgs_files)]
    if (length(missing_pgs) > 0L)
      stop(sprintf("PGS file(s) not found:\n  %s", paste(missing_pgs, collapse = "\n  ")))
  }
  if (verbose) {
    cat(sprintf("Found %d PGS model file(s):\n", length(pgs_files)))
    cat(paste0("  ", basename(pgs_files), "\n"), sep = "")
  }

  # ---- Discover genotype files and resolve format -----------------------------
  geno <- .discover_geno_files(geno_dir, format, verbose)
  if (verbose) {
    cat(sprintf("\nUsing %s genotype files from %s (%d chromosome(s) available):\n",
                geno$format[1L], geno_dir, nrow(geno)))
    cat(sprintf("  chr%-3s <- %s\n", geno$chrom, basename(geno$path)), sep = "")
  }

  # ---- Read PGS models ----------------------------------------------------------
  model_ids <- sub("_.*", "", basename(pgs_files))
  pgs_list  <- setNames(vector("list", length(pgs_files)), model_ids)
  for (i in seq_along(pgs_files)) {
    pgs <- readPGSmodel(pgs_files[i], verbose = FALSE)
    pgs$chr_name <- .norm_chrom(pgs$chr_name)
    pgs_list[[i]] <- pgs[, c("rsID", "chr_name", "chr_position", "effect_allele", "effect_weight")]
    if (verbose) cat(sprintf("  %-12s : %d SNPs\n", model_ids[i], nrow(pgs_list[[i]])))
  }

  # ---- Determine every chromosome to account for -------------------------------
  model_chroms <- unique(unlist(lapply(pgs_list, `[[`, "chr_name")))
  all_chroms   <- .chrom_sort(unique(c(geno$chrom, model_chroms)))
  needed_chroms <- .chrom_sort(model_chroms)   # only these are actually queried

  missing_chroms <- setdiff(needed_chroms, geno$chrom)
  for (chrom in missing_chroms)
    warning(sprintf("No genotype file found for chromosome %s - all SNPs on chr%s will be treated as unmatched.",
                    chrom, chrom))

  # ---- Samples to score: IDs common to every file that will be queried --------
  file_samples <- attr(geno, "samples")
  used_paths   <- unique(geno$path[geno$chrom %in% setdiff(needed_chroms, missing_chroms)])
  if (length(used_paths) > 0L) {
    common           <- .common_samples(file_samples[used_paths], verbose)
    scoring_samples  <- common$common
    excluded_samples <- common$excluded
  } else {   # nothing will be queried (every needed chromosome lacks a file)
    scoring_samples  <- file_samples[[1L]]
    excluded_samples <- character(0)
  }

  # ---- Per-model accumulators ----------------------------------------------------
  zero_chrom_vec <- setNames(integer(length(all_chroms)), all_chroms)
  empty_rsid_list <- setNames(vector("list", length(all_chroms)), all_chroms)
  for (nm in names(empty_rsid_list)) empty_rsid_list[[nm]] <- character(0)

  prs_total           <- setNames(vector("list", length(pgs_list)), model_ids)
  unmatched_by_chr     <- setNames(rep(list(zero_chrom_vec),  length(pgs_list)), model_ids)
  unmatched_rsIDs      <- setNames(rep(list(empty_rsid_list), length(pgs_list)), model_ids)
  excluded_chr_counts  <- setNames(rep(list(zero_chrom_vec),  length(pgs_list)), model_ids)

  for (i in seq_along(pgs_list)) prs_total[[i]] <- setNames(numeric(length(scoring_samples)), scoring_samples)

  # ---- Missing-chromosome SNPs are entirely unmatched --------------------------
  for (chrom in missing_chroms) {
    for (m in model_ids) {
      snps_m <- pgs_list[[m]][pgs_list[[m]]$chr_name == chrom, ]
      if (nrow(snps_m) == 0L) next
      unmatched_by_chr[[m]][chrom]  <- nrow(snps_m)
      unmatched_rsIDs[[m]][[chrom]] <- snps_m$rsID
    }
  }

  # ---- Score chromosome by chromosome -------------------------------------------
  t_start <- proc.time()[["elapsed"]]

  for (chrom in setdiff(needed_chroms, missing_chroms)) {
    geno_row     <- geno[geno$chrom == chrom, ]
    models_here  <- model_ids[vapply(pgs_list, function(p) any(p$chr_name == chrom), logical(1L))]

    positions <- sort(unique(unlist(lapply(pgs_list[models_here], function(p) p$chr_position[p$chr_name == chrom]))))
    if (verbose)
      cat(sprintf("\n--- chr%s : %d unique position(s) across %d model(s) ---\n",
                  chrom, length(positions), length(models_here)))

    dosage_chrom <- if (geno_row$format[1L] == "vcf") {
      .query_vcf_dosage(geno_row$path[1L], geno_row$contig[1L], positions, batch_size, verbose)
    } else {
      .query_bd_dosage(geno_row$path[1L], geno_row$contig[1L], positions, verbose)
    }

    for (m in models_here) {
      snps_m <- pgs_list[[m]][pgs_list[[m]]$chr_name == chrom, ]

      if (is.null(dosage_chrom) || nrow(dosage_chrom) == 0L) {
        unmatched_by_chr[[m]][chrom]  <- nrow(snps_m)
        unmatched_rsIDs[[m]][[chrom]] <- snps_m$rsID
        next
      }

      found_pos      <- snps_m$chr_position %in% dosage_chrom$POS
      unmatched_by_chr[[m]][chrom]  <- sum(!found_pos)
      unmatched_rsIDs[[m]][[chrom]] <- snps_m$rsID[!found_pos]
      snps_found <- snps_m[found_pos, ]
      if (nrow(snps_found) == 0L) next

      merged <- merge(dosage_chrom, snps_found,
                      by.x = "POS", by.y = "chr_position", sort = FALSE)

      is_ref    <- merged$effect_allele == merged$REF
      is_alt    <- merged$effect_allele == merged$ALT
      n_neither <- sum(!is_ref & !is_alt)
      excluded_chr_counts[[m]][chrom] <- n_neither
      merged <- merged[is_ref | is_alt, ]
      if (nrow(merged) == 0L) next
      is_ref <- merged$effect_allele == merged$REF

      # REF-effect formula: weight * (2 - dosage) = -weight * dosage + 2 * weight
      adj_weight     <- ifelse(is_ref, -merged$effect_weight, merged$effect_weight)
      ref_correction <- 2 * sum(merged$effect_weight[is_ref])

      # Only the samples common to every file; columns are in scoring_samples order,
      # so this chromosome's scores line up with the running totals.
      dos_mat     <- as.matrix(merged[, scoring_samples, drop = FALSE])
      prs_chrom   <- drop(t(dos_mat) %*% adj_weight) + ref_correction
      prs_total[[m]] <- prs_total[[m]] + prs_chrom

      if (verbose) {
        cat(sprintf("  %-12s : %d matched, %d excluded (allele mismatch), %d unmatched\n",
                    m, nrow(merged), n_neither, sum(!found_pos)))
      }
    }
  }

  t_elapsed <- proc.time()[["elapsed"]] - t_start
  if (verbose) cat(sprintf("\nScoring complete : %.1f s\n", t_elapsed))

  # ---- Assemble return value and save -------------------------------------------
  results <- setNames(vector("list", length(pgs_list)), model_ids)
  for (m in model_ids) {
    results[[m]] <- list(
      prs                 = prs_total[[m]],
      unmatched_by_chr    = unmatched_by_chr[[m]],
      unmatched_rsIDs     = unmatched_rsIDs[[m]],
      excluded_chr_counts = excluded_chr_counts[[m]],
      excluded_samples    = excluded_samples
    )
  }

  if (!is.null(output_dir)) {
    for (m in model_ids) {
      out_file <- file.path(output_dir, sprintf("pgs_%s_prs.rds", m))
      saveRDS(results[[m]], out_file)
      if (verbose) cat(sprintf("Saved : %s\n", out_file))
    }
  }

  invisible(results)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Strip a leading "chr" (case-insensitive) and coerce to character.
.norm_chrom <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}

# Sort chromosome labels: numeric ones in numeric order, then the rest
# alphabetically (X, Y, MT, ...).
.chrom_sort <- function(chroms) {
  chroms   <- unique(chroms)
  num      <- suppressWarnings(as.integer(chroms))
  is_num   <- !is.na(num)
  c(chroms[is_num][order(num[is_num])], sort(chroms[!is_num]))
}

# Data files whose index file (.tbi / .bdi) exists in dir. The index is what
# identifies a genotype file, so the data file names themselves are arbitrary.
# An index with no data file next to it is skipped with a warning.
.find_indexed_files <- function(dir, index_ext, data_ext) {
  all_files <- list.files(dir, full.names = TRUE)
  idx       <- all_files[endsWith(all_files, paste0(data_ext, index_ext))]
  data_files <- substr(idx, 1L, nchar(idx) - nchar(index_ext))
  orphan <- !file.exists(data_files)
  if (any(orphan))
    warning(sprintf("Index file(s) without a matching data file were ignored:\n  %s",
                    paste(basename(idx[orphan]), collapse = "\n  ")))
  sort(data_files[!orphan])
}

# BinaryDosage files in dir. Format 5 files have a companion .bdi; Format 4
# files are a single self-contained file, so the .bdose file itself is what
# identifies them. A .bdi with no .bdose next to it is skipped with a warning.
.find_bdose_files <- function(dir) {
  all_files <- list.files(dir, full.names = TRUE)
  bdose     <- all_files[endsWith(all_files, ".bdose")]
  bdi       <- all_files[endsWith(all_files, ".bdose.bdi")]
  orphan    <- !file.exists(substr(bdi, 1L, nchar(bdi) - nchar(".bdi")))
  if (any(orphan))
    warning(sprintf("Index file(s) without a matching data file were ignored:\n  %s",
                    paste(basename(bdi[orphan]), collapse = "\n  ")))
  sort(bdose)
}

# A saved getbdinfo() result is only trusted if it has the pieces getsnp()
# needs and its SNP offsets are consistent with the data file it is paired with.
.valid_bdinfo <- function(x, bdose_path) {
  is.list(x) && inherits(x, "genetic-info") &&
    all(c("filename", "samples", "snps", "indices") %in% names(x)) &&
    is.data.frame(x$snps) && length(x$indices) == nrow(x$snps) &&
    (length(x$indices) == 0L || max(x$indices) < file.size(bdose_path))
}

# BinaryDosage file information (the getbdinfo() result) for one .bdose file.
#
# Format 5 keeps its information in the companion .bdi, which getbdinfo() reads
# cheaply. For earlier formats getbdinfo() has to parse the whole data file, so
# if <name>.bdinfo (a saveRDS() of a previous getbdinfo() result) sits next to
# <name>.bdose it is used instead. The stored 'filename' is always replaced
# with the actual path: getsnp() reads from it, and a stale path (files moved
# or copied since the .bdinfo was made) silently returns zeros, not an error.
# A .bdinfo that cannot be read or does not match the data file is ignored
# with a warning and the data file is parsed instead.
.get_bdinfo <- function(path) {
  stopifnot("BinaryDosage is not installed" = requireNamespace("BinaryDosage", quietly = TRUE))

  bdinfo_file <- sub("\\.bdose$", ".bdinfo", path)
  if (!file.exists(paste0(path, ".bdi")) && file.exists(bdinfo_file)) {
    saved <- tryCatch(readRDS(bdinfo_file), error = function(e) NULL)
    if (.valid_bdinfo(saved, path)) {
      saved$filename <- path
      attr(saved, "info_source") <- "bdinfo"
      return(saved)
    }
    warning(sprintf("Ignoring %s: it could not be read as a getbdinfo() result matching %s. Reading %s instead.",
                    basename(bdinfo_file), basename(path), basename(path)), call. = FALSE)
  }

  bdinfo <- tryCatch(
    BinaryDosage::getbdinfo(bdfiles = path),
    error = function(e)
      stop(sprintf(paste0("Could not read BinaryDosage file %s (Format 4 and 5 files are supported; ",
                          "Format 5 needs its %s.bdi alongside): %s"),
                   path, basename(path), conditionMessage(e)), call. = FALSE)
  )
  attr(bdinfo, "info_source") <- if (file.exists(paste0(path, ".bdi"))) "bdi" else "parsed"
  bdinfo
}

# What one genotype file contains: the chromosome label(s) exactly as stored in
# the file (e.g. "chr1" or "1"), its sample IDs, where the information came
# from (VCF: "tbi"; BinaryDosage: "bdi", "bdinfo", or "parsed" = read from the
# whole data file), and the BinaryDosage format number (NA for VCF).
.file_info <- function(path, format) {
  if (format == "vcf")
    return(list(contigs = tabixr::vcf_seqnames(path), samples = tabixr::vcf_samples(path),
                info = "tbi", version = NA_integer_))

  bdinfo  <- .get_bdinfo(path)
  contigs <- unique(as.character(bdinfo$snps$chromosome))
  if (any(is.na(contigs) | !nzchar(contigs)))
    stop(sprintf("BinaryDosage file %s has SNPs with no chromosome recorded", path), call. = FALSE)
  version <- bdinfo$additionalinfo$format
  list(contigs = contigs, samples = as.character(bdinfo$samples$sid),
       info = attr(bdinfo, "info_source"),
       version = if (is.null(version)) NA_integer_ else as.integer(version))
}

# Discover genotype files in geno_dir, resolve which format to use, and work
# out which chromosome(s) each file holds. Returns a data.frame with one row
# per (file, chromosome): chrom (label without "chr"), contig (label as stored
# in the file), file, format, version, info, n_samples, path. Duplicate
# chromosomes are NOT rejected here (see .discover_geno_files). Attributes:
# "samples" (list of sample IDs per file path) and "unused" (list with the
# format and paths of any files of the other format that were not used).
.scan_geno_files <- function(geno_dir, format, verbose) {
  vcf_files   <- .find_indexed_files(geno_dir, ".tbi", ".vcf.gz")
  bdose_files <- .find_bdose_files(geno_dir)
  has_vcf   <- length(vcf_files)   > 0L
  has_bdose <- length(bdose_files) > 0L

  if (!is.null(format)) {
    chosen <- format
    if (chosen == "vcf"   && !has_vcf)
      stop(sprintf("format='vcf' requested but no indexed VCF files (*.vcf.gz + *.vcf.gz.tbi) found in %s", geno_dir))
    if (chosen == "bdose" && !has_bdose)
      stop(sprintf("format='bdose' requested but no BinaryDosage files (*.bdose) found in %s", geno_dir))
  } else if (has_vcf && has_bdose) {
    chosen <- "bdose"
    if (verbose)
      cat(sprintf("Both VCF and BinaryDosage genotype files found in %s - using BinaryDosage.\n", geno_dir))
  } else if (has_bdose) {
    chosen <- "bdose"
  } else if (has_vcf) {
    chosen <- "vcf"
  } else {
    stop(sprintf("No genotype files (*.vcf.gz + .tbi, or *.bdose) found in %s", geno_dir))
  }

  files  <- if (chosen == "vcf") vcf_files else bdose_files
  infos  <- lapply(files, .file_info, format = chosen)
  geno   <- do.call(rbind, lapply(seq_along(files), function(i) {
    if (length(infos[[i]]$contigs) == 0L)
      stop(sprintf("No chromosomes found in %s", files[i]))
    data.frame(chrom = .norm_chrom(infos[[i]]$contigs), contig = infos[[i]]$contigs,
               file = basename(files[i]), format = chosen, version = infos[[i]]$version,
               info = infos[[i]]$info, n_samples = length(infos[[i]]$samples),
               path = files[i], stringsAsFactors = FALSE)
  }))
  attr(geno, "samples") <- setNames(lapply(infos, `[[`, "samples"), files)
  attr(geno, "unused")  <- list(format = if (chosen == "vcf") "bdose" else "vcf",
                                paths  = if (chosen == "vcf") bdose_files else vcf_files)
  geno
}

# .scan_geno_files() plus the check compute_prs() needs: stops if a chromosome
# appears in more than one file. Returns the rows ordered by chromosome.
.discover_geno_files <- function(geno_dir, format, verbose) {
  geno   <- .scan_geno_files(geno_dir, format, verbose)
  chosen <- geno$format[1L]

  dups <- unique(geno$chrom[duplicated(geno$chrom)])
  if (length(dups) > 0L) {
    detail <- vapply(dups, function(ch)
      sprintf("  chromosome %s: %s", ch, paste(basename(geno$path[geno$chrom == ch]), collapse = ", ")),
      character(1L))
    stop(sprintf("More than one %s file contains the same chromosome in %s:\n%s",
                 chosen, geno_dir, paste(detail, collapse = "\n")))
  }

  samples <- attr(geno, "samples")
  geno    <- geno[order(match(geno$chrom, .chrom_sort(geno$chrom))), ]
  attr(geno, "samples") <- samples
  geno
}

# Sample IDs to score: those present in EVERY file that will be queried.
# Files often differ slightly (subjects dropped from some chromosomes for poor
# imputation quality), and scores can only be summed across chromosomes over
# subjects every file has. Order follows the first file. Stops if the
# intersection is empty, which means the files are not from the same cohort or
# their sample IDs are formatted differently.
# 'samples_by_file' is a named list (name = path) of sample ID vectors.
.common_samples <- function(samples_by_file, verbose) {
  common   <- Reduce(intersect, samples_by_file)
  excluded <- setdiff(unique(unlist(samples_by_file, use.names = FALSE)), common)

  if (length(common) == 0L)
    stop(paste0("No sample IDs are common to all the genotype files used, so nothing can be scored. ",
                "Check that the files come from the same cohort and that sample IDs are formatted ",
                "the same way in each:\n",
                paste(sprintf("  %s: %d sample(s), e.g. %s", basename(names(samples_by_file)),
                              lengths(samples_by_file),
                              vapply(samples_by_file, function(s) paste(head(s, 2L), collapse = ", "), character(1L))),
                      collapse = "\n")),
         call. = FALSE)

  if (verbose && length(excluded) > 0L)
    cat(sprintf("\nSample IDs: scoring the %d sample(s) present in all %d genotype file(s) used; %d sample(s) missing from at least one file are excluded.\n",
                length(common), length(samples_by_file), length(excluded)))

  list(common = common, excluded = excluded)
}

# ---- VCF dosage retrieval (batched; subprocess-isolated off Windows) --------

.query_vcf_dosage <- function(vcf_path, chrom_internal, positions, batch_size, verbose) {
  if (length(positions) == 0L) return(NULL)
  batches <- split(positions, ceiling(seq_along(positions) / batch_size))

  dosage_all <- if (.Platform$OS.type != "windows")
    .query_batches_subprocess(vcf_path, chrom_internal, batches, verbose)
  else
    .query_batches_inprocess(vcf_path, chrom_internal, batches, verbose)

  dosage_all
}

.query_batches_inprocess <- function(vcf_path, chrom, batches, verbose) {
  n_batches  <- length(batches)
  dosage_all <- NULL
  # tabixr's data.frame can rename sample columns ("1001" -> "X1001"); the header has the true IDs.
  sample_ids <- tabixr::vcf_samples(vcf_path)

  for (b in seq_along(batches)) {
    batch_pos <- batches[[b]]
    if (verbose)
      cat(sprintf("  Batch %d/%d : %d positions", b, n_batches, length(batch_pos)))

    vcf_hits <- tabixr::query_vcf_positions(vcf_path, chrom, batch_pos)
    matched  <- vcf_hits[vcf_hits$POS %in% batch_pos, ]
    rm(vcf_hits)

    if (nrow(matched) == 0L) {
      if (verbose) cat(" - 0 hits\n")
      next
    }
    if (verbose) cat(sprintf(" - %d VCF row(s)\n", nrow(matched)))

    dos        <- extract_dosage(matched, sample_ids)
    rm(matched)
    dosage_all <- if (is.null(dosage_all)) dos else rbind(dosage_all, dos)
    rm(dos)
  }

  dosage_all
}

.query_batches_subprocess <- function(vcf_path, chrom, batches, verbose) {
  n_batches       <- length(batches)
  batch_dos_files <- character(n_batches)
  tmp             <- tempdir()
  worker_script   <- system.file("extdata", "batch_worker.R", package = "pgscorer")
  rscript         <- file.path(R.home("bin"), "Rscript")

  for (b in seq_along(batches)) {
    pos_rds <- file.path(tmp, sprintf("_pgscorer_pos_%02d.rds", b))
    dos_rds <- file.path(tmp, sprintf("_pgscorer_dos_%02d.rds", b))
    batch_dos_files[b] <- dos_rds

    saveRDS(batches[[b]], pos_rds)
    if (verbose)
      cat(sprintf("  Batch %d/%d : %d positions\n", b, n_batches, length(batches[[b]])))

    status <- system2(
      rscript,
      args = c(shQuote(worker_script),
               shQuote(pos_rds), shQuote(vcf_path), chrom, shQuote(dos_rds))
    )
    file.remove(pos_rds)
    if (status != 0L)
      stop(sprintf("Batch %d subprocess failed (exit status %d)", b, status))
  }

  dosage_parts <- Filter(Negate(is.null), lapply(batch_dos_files, readRDS))
  for (f in batch_dos_files) if (file.exists(f)) file.remove(f)
  dosage_all <- if (length(dosage_parts) > 0L) do.call(rbind, dosage_parts) else NULL
  rm(dosage_parts)
  invisible(gc())
  dosage_all
}

# ---- BinaryDosage dosage retrieval (in-process, random-access by SNP) -------

.query_bd_dosage <- function(bdose_path, contig, positions, verbose) {
  if (length(positions) == 0L) return(NULL)
  stopifnot("BinaryDosage is not installed" = requireNamespace("BinaryDosage", quietly = TRUE))

  bdinfo <- .get_bdinfo(bdose_path)
  snps   <- bdinfo$snps
  idx    <- which(snps$chromosome == contig & snps$location %in% positions)

  if (verbose)
    cat(sprintf("  %d/%d SNP(s) in %s match requested positions\n",
                length(idx), nrow(snps), basename(bdose_path)))
  if (length(idx) == 0L) return(NULL)

  samp_ids <- bdinfo$samples$sid
  dos_mat  <- t(vapply(idx, function(i) {
    BinaryDosage::getsnp(bdinfo, i, dosageonly = TRUE)$dosage
  }, FUN.VALUE = numeric(length(samp_ids))))
  colnames(dos_mat) <- samp_ids

  cbind(
    data.frame(
      POS = snps$location[idx],
      ID  = snps$snpid[idx],
      REF = snps$reference[idx],
      ALT = snps$alternate[idx],
      stringsAsFactors = FALSE
    ),
    as.data.frame(dos_mat, stringsAsFactors = FALSE)
  )
}
