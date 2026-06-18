#' Read a PGS Catalog Scoring File
#'
#' Reads a PGS Catalog scoring file (.txt or .txt.gz), skips \code{##} comment
#' lines, normalises harmonised columns (\code{hm_chr}/\code{hm_pos}) to
#' \code{chr_name}/\code{chr_position}, fills missing \code{other_allele} and
#' \code{rsID} fields, and returns a six-column data.frame.
#'
#' @param file_path Path to the PGS Catalog file (.txt or .txt.gz).
#' @param verbose Logical; if \code{TRUE} (default) print progress messages.
#'
#' @return A data.frame with columns \code{rsID}, \code{chr_name},
#'   \code{chr_position}, \code{effect_allele}, \code{other_allele},
#'   \code{effect_weight}.
#' @export
readPGSmodel <- function(file_path, verbose = TRUE) {
  if (!is.character(file_path) || length(file_path) != 1L)
    stop("file_path must be a single character string")
  if (!file.exists(file_path))
    stop(sprintf("File does not exist: %s", file_path))
  if (verbose) message(sprintf("Reading PGS model from: %s", file_path))

  pgs_data <- tryCatch({
    con    <- file(file_path, "r")
    n_skip <- 0L
    while (TRUE) {
      line <- readLines(con, n = 1L)
      if (length(line) == 0L || !grepl("^##", line)) break
      n_skip <- n_skip + 1L
    }
    close(con)

    dt <- data.table::fread(file_path,
                            header      = TRUE,
                            sep         = "\t",
                            skip        = n_skip,
                            quote       = "",
                            check.names = FALSE)
    as.data.frame(dt)
  }, error = function(e) {
    stop(sprintf("Error reading %s: %s", file_path, e$message))
  })

  if (verbose)
    message(sprintf("Read %d rows, %d columns", nrow(pgs_data), ncol(pgs_data)))

  # Prefer harmonised chromosome / position columns when present.
  if ("hm_chr" %in% names(pgs_data)) {
    pgs_data$chr_name <- pgs_data$hm_chr
    if (verbose) message("Using hm_chr for chr_name")
  } else if (!"chr_name" %in% names(pgs_data)) {
    if (verbose) message("Warning: neither hm_chr nor chr_name found")
  }

  if ("hm_pos" %in% names(pgs_data)) {
    pgs_data$chr_position <- pgs_data$hm_pos
    if (verbose) message("Using hm_pos for chr_position")
  } else if (!"chr_position" %in% names(pgs_data)) {
    if (verbose) message("Warning: neither hm_pos nor chr_position found")
  }

  # Drop rows with unusable positions.
  if ("chr_position" %in% names(pgs_data)) {
    bad <- is.na(pgs_data$chr_position) | pgs_data$chr_position == 0
    n_bad <- sum(bad, na.rm = TRUE)
    if (n_bad > 0L) {
      pgs_data <- pgs_data[!bad, ]
      if (verbose) message(sprintf("Removed %d rows with NA/0 chr_position", n_bad))
    }
  }

  # Fill other_allele.
  if (!"other_allele" %in% names(pgs_data)) {
    pgs_data$other_allele <-
      if ("hm_inferOtherAllele" %in% names(pgs_data))
        pgs_data$hm_inferOtherAllele
      else
        NA_character_
  }

  # Fill rsID.
  make_composite <- function(rows) {
    paste(rows$chr_name, rows$chr_position, rows$effect_allele,
          ifelse(is.na(rows$other_allele), "X", rows$other_allele),
          sep = ":")
  }
  if ("rsID" %in% names(pgs_data)) {
    empty <- pgs_data$rsID == "" | is.na(pgs_data$rsID)
    if (any(empty)) pgs_data$rsID[empty] <- make_composite(pgs_data[empty, ])
  } else {
    pgs_data$rsID <- make_composite(pgs_data)
  }

  required <- c("rsID", "chr_name", "chr_position",
                 "effect_allele", "other_allele", "effect_weight")
  missing  <- setdiff(required, names(pgs_data))
  if (length(missing) > 0L)
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))

  pgs_data[, required]
}
