#' Compute Polygenic Risk Scores from PGS Catalog Models
#'
#' Unions SNP positions from one or more PGS Catalog scoring files, queries a
#' BGZF-compressed VCF file once per position (in batches), then computes a
#' PRS for every model from the shared dosage table.
#'
#' On Windows the VCF is queried in-process in sequential batches.  On other
#' platforms each batch is run in its own \code{Rscript} subprocess, isolating
#' batch memory so peak footprint is bounded to one batch at a time.
#'
#' @param vcf_path   Path to the BGZF-compressed VCF (\code{.vcf.gz}).  A
#'   matching \code{.tbi} Tabix index must exist alongside it.
#' @param pgs_files  Character vector of paths to PGS Catalog scoring files
#'   (\code{.txt.gz}).  If \code{NULL} (default), all files matching
#'   \code{^PGS.*\\.txt\\.gz$} in \code{pgs_dir} are used.
#' @param pgs_dir    Directory searched for PGS files when \code{pgs_files} is
#'   \code{NULL}.  Defaults to the current working directory.
#' @param chrom      Chromosome label to restrict scoring (e.g. \code{"21"}).
#'   \code{NULL} (default) auto-detects from the VCF index; an error is raised
#'   when the VCF contains more than one sequence and \code{chrom} is
#'   \code{NULL}.
#' @param batch_size Integer; positions per VCF query batch.  Default 10000.
#' @param output_dir Directory in which to write one \code{.rds} file per model
#'   named \code{pgs_chr<chrom>_<model>_prs_combined.rds}.  Set to \code{NULL}
#'   to skip saving.  Defaults to the current working directory.
#' @param verbose    Logical; if \code{TRUE} (default) print progress and
#'   summary messages.
#'
#' @return A named list of numeric vectors (invisibly), one element per model.
#'   Each vector is named by sample identifier and contains the PRS value.
#'
#' @export
compute_prs <- function(vcf_path,
                        pgs_files  = NULL,
                        pgs_dir    = ".",
                        chrom      = NULL,
                        batch_size = 10000L,
                        output_dir = ".",
                        verbose    = TRUE) {

  stopifnot(
    "tabixr is not installed"          = requireNamespace("tabixr",     quietly = TRUE),
    "data.table is not installed"      = requireNamespace("data.table", quietly = TRUE),
    "vcf_path must be a single string" = is.character(vcf_path) && length(vcf_path) == 1L,
    "VCF file not found"               = file.exists(vcf_path),
    ".tbi index not found"             = file.exists(paste0(vcf_path, ".tbi")),
    "batch_size must be a positive integer" = is.numeric(batch_size) && batch_size >= 1L
  )
  batch_size <- as.integer(batch_size)

  # ---- Discover PGS files ----------------------------------------------------
  if (is.null(pgs_files)) {
    pgs_files <- sort(list.files(pgs_dir, pattern = "^PGS.*\\.txt\\.gz$",
                                 full.names = TRUE))
    if (length(pgs_files) == 0L)
      stop(sprintf("No PGS model files (^PGS.*\\.txt\\.gz$) found in: %s", pgs_dir))
  } else {
    missing_pgs <- pgs_files[!file.exists(pgs_files)]
    if (length(missing_pgs) > 0L)
      stop(sprintf("PGS file(s) not found:\n  %s",
                   paste(missing_pgs, collapse = "\n  ")))
  }
  if (verbose) {
    cat(sprintf("Found %d PGS model file(s):\n", length(pgs_files)))
    cat(paste0("  ", basename(pgs_files), "\n"), sep = "")
  }

  # ---- Detect chromosome -----------------------------------------------------
  seqs <- tabixr::vcf_seqnames(vcf_path)
  if (is.null(chrom)) {
    if (length(seqs) == 1L) {
      chrom <- seqs
      if (verbose) cat(sprintf("Detected chromosome : %s\n", chrom))
    } else {
      stop(sprintf(
        "VCF contains %d chromosomes (%s).\nSet chrom= to select one.",
        length(seqs), paste(seqs, collapse = ", ")
      ))
    }
  } else if (!chrom %in% seqs) {
    stop(sprintf("chrom='%s' not found in VCF index.\nAvailable: %s",
                 chrom, paste(seqs, collapse = ", ")))
  }

  # ---- Read PGS models -------------------------------------------------------
  if (verbose) cat(sprintf("\nReading PGS models (chr%s only)...\n", chrom))
  model_ids <- sub("_.*", "", basename(pgs_files))
  pgs_list  <- setNames(vector("list", length(pgs_files)), model_ids)

  for (i in seq_along(pgs_files)) {
    pgs           <- readPGSmodel(pgs_files[i], verbose = FALSE)
    pgs           <- pgs[pgs$chr_name == chrom,
                         c("chr_position", "effect_allele", "effect_weight")]
    pgs_list[[i]] <- pgs
    if (verbose) cat(sprintf("  %-12s : %d SNPs\n", model_ids[i], nrow(pgs)))
  }

  # ---- Union positions -------------------------------------------------------
  all_positions <- sort(unique(unlist(lapply(pgs_list, `[[`, "chr_position"))))
  if (verbose)
    cat(sprintf("\nUnique positions across all models : %d\n", length(all_positions)))

  # ---- Query VCF in batches --------------------------------------------------
  batches   <- split(all_positions,
                     ceiling(seq_along(all_positions) / batch_size))
  n_batches <- length(batches)

  use_subprocess <- .Platform$OS.type != "windows"
  if (verbose) {
    mode_label <- if (use_subprocess) "subprocess per batch" else "in-process"
    cat(sprintf("Querying VCF in %d batch(es) of up to %d positions (%s)...\n",
                n_batches, batch_size, mode_label))
  }

  t_vcf_start <- proc.time()[["elapsed"]]
  dosage_all  <- if (use_subprocess)
    .query_batches_subprocess(vcf_path, chrom, batches, verbose)
  else
    .query_batches_inprocess(vcf_path, chrom, batches, verbose)
  t_vcf <- proc.time()[["elapsed"]] - t_vcf_start

  heap_mb <- function() { g <- gc(verbose = FALSE); sum(g[, "used"]) * 8 / 1e6 }
  if (verbose)
    cat(sprintf("VCF query complete : %.1f s  |  %d dosage row(s)  |  %.1f MB heap\n",
                t_vcf, nrow(dosage_all), heap_mb()))

  # ---- Score each model ------------------------------------------------------
  fixed_cols    <- c("POS", "ID", "REF", "ALT")
  samp_cols     <- setdiff(names(dosage_all), fixed_cols)
  summary_rows  <- vector("list", length(pgs_list))
  prs_list      <- setNames(vector("list", length(pgs_list)), model_ids)
  t_score_start <- proc.time()[["elapsed"]]

  for (i in seq_along(pgs_list)) {
    model_id <- model_ids[i]
    pgs      <- pgs_list[[i]]
    if (verbose) cat(sprintf("\n--- %s ---\n", model_id))

    merged <- merge(dosage_all, pgs,
                    by.x = "POS", by.y = "chr_position", sort = FALSE)

    is_ref    <- merged$effect_allele == merged$REF
    is_alt    <- merged$effect_allele == merged$ALT
    n_neither <- sum(!is_ref & !is_alt)
    n_matched <- nrow(merged)
    if (n_neither > 0L)
      message(sprintf("  WARNING: %d row(s) where effect_allele matches neither REF nor ALT — excluded.",
                      n_neither))
    merged <- merged[is_ref | is_alt, ]
    is_ref <- merged$effect_allele == merged$REF

    # REF-effect formula: weight * (2 - dosage) = -weight * dosage + 2 * weight
    adj_weight     <- ifelse(is_ref, -merged$effect_weight, merged$effect_weight)
    ref_correction <- 2 * sum(merged$effect_weight[is_ref])

    dos_mat    <- as.matrix(merged[, samp_cols])
    prs        <- drop(t(dos_mat) %*% adj_weight) + ref_correction
    names(prs) <- samp_cols

    if (verbose) {
      cat(sprintf("  Model SNPs    : %d\n", nrow(pgs)))
      cat(sprintf("  VCF rows      : %d\n", n_matched))
      cat(sprintf("  Excluded      : %d\n", n_neither))
      cat(sprintf("  Scored        : %d\n", nrow(merged)))
      cat(sprintf("  REF-effect    : %d\n", sum(is_ref)))
      cat(sprintf("  ALT-effect    : %d\n", sum(!is_ref)))
    }

    prs_list[[i]]     <- prs
    summary_rows[[i]] <- data.frame(
      model      = model_id,
      model_snps = nrow(pgs),
      matched    = nrow(merged),
      ref_effect = sum(is_ref),
      alt_effect = sum(!is_ref),
      subjects   = length(prs),
      stringsAsFactors = FALSE
    )
  }

  t_score <- proc.time()[["elapsed"]] - t_score_start

  # ---- Summary ---------------------------------------------------------------
  if (verbose) {
    cat(sprintf("\n=== Timing ===\n"))
    cat(sprintf("VCF query : %.1f s\n", t_vcf))
    cat(sprintf("Scoring   : %.1f s\n", t_score))
    cat(sprintf("Total     : %.1f s\n", t_vcf + t_score))
    cat("\n")
    print(do.call(rbind, summary_rows), row.names = FALSE)
  }

  # ---- Save results ----------------------------------------------------------
  if (!is.null(output_dir)) {
    if (verbose) cat("\n")
    for (i in seq_along(prs_list)) {
      out_file <- file.path(
        output_dir,
        sprintf("pgs_chr%s_%s_prs_combined.rds", chrom, model_ids[i])
      )
      saveRDS(prs_list[[i]], out_file)
      if (verbose) cat(sprintf("Saved : %s\n", out_file))
    }
  }

  invisible(prs_list)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.query_batches_inprocess <- function(vcf_path, chrom, batches, verbose) {
  n_batches  <- length(batches)
  dosage_all <- NULL

  for (b in seq_along(batches)) {
    batch_pos <- batches[[b]]
    if (verbose)
      cat(sprintf("  Batch %d/%d : %d positions", b, n_batches, length(batch_pos)))

    vcf_hits <- tabixr::query_vcf_positions(vcf_path, chrom, batch_pos)
    matched  <- vcf_hits[vcf_hits$POS %in% batch_pos, ]
    rm(vcf_hits)

    if (nrow(matched) == 0L) {
      if (verbose) cat(" — 0 hits\n")
      next
    }
    if (verbose) cat(sprintf(" — %d VCF row(s)\n", nrow(matched)))

    dos        <- extract_dosage(matched)
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
      cat(sprintf("  Batch %d/%d : %d positions\n",
                  b, n_batches, length(batches[[b]])))

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
  dosage_all <- do.call(rbind, dosage_parts)
  rm(dosage_parts)
  invisible(gc())
  dosage_all
}
