#' Extract DS (Dosage) Values from a VCF Data Frame
#'
#' Takes the data.frame returned by \code{\link[tabixr]{query_vcf_positions}}
#' (or \code{\link[tabixr]{query_vcf}}) and returns a slimmed data.frame with
#' columns \code{POS}, \code{ID}, \code{REF}, \code{ALT}, and one numeric
#' column per sample containing the \code{DS} dosage value.
#'
#' The \code{DS} field position is read from the \code{FORMAT} column of the
#' first row, so the function works regardless of FORMAT field order
#' (\code{GT:DS:GP}, \code{GT:GP:DS}, etc.). Rows where \code{DS} is missing
#' (\code{"."}) are returned as \code{NA}.
#'
#' Sample columns of the data.frame from \code{tabixr} can have altered names:
#' R's \code{data.frame()} turns an ID such as \code{"1001"} into
#' \code{"X1001"}. Pass the true IDs (in header order, as returned by
#' \code{tabixr::vcf_samples()}) as \code{sample_ids} to name the output columns
#' correctly.
#'
#' @param df A data.frame as returned by \code{tabixr::query_vcf_positions()}.
#' @param sample_ids Optional character vector of sample IDs to use as the names
#'   of the sample columns, one per sample column of \code{df} and in the same
#'   order. Default \code{NULL} uses the column names of \code{df}.
#'
#' @return A data.frame with columns \code{POS} (integer), \code{ID},
#'   \code{REF}, \code{ALT} (all character), and one numeric column per sample.
#' @export
extract_dosage <- function(df, sample_ids = NULL) {
  fmt_col   <- match("FORMAT", names(df))
  samp_cols <- names(df)[(fmt_col + 1L):ncol(df)]
  if (is.null(sample_ids)) {
    out_names <- samp_cols
  } else if (length(sample_ids) == length(samp_cols)) {
    out_names <- as.character(sample_ids)
  } else {
    stop(sprintf("sample_ids has %d element(s) but the data.frame has %d sample column(s)",
                 length(sample_ids), length(samp_cols)))
  }

  if (nrow(df) == 0L) {
    out <- df[, c("POS", "ID", "REF", "ALT")]
    for (col in out_names) out[[col]] <- numeric(0)
    return(out)
  }

  fmt      <- df$FORMAT[1L]
  fields   <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
  ds_idx   <- match("DS", fields)
  if (is.na(ds_idx))
    stop("DS field not found in FORMAT: ", fmt)

  # Extract the ds_idx-th colon-delimited token from each sample column.
  pattern <- paste0("^(?:[^:]*:){", ds_idx - 1L, "}([^:]*).*$")
  ds_vals <- lapply(samp_cols, function(col) {
    as.numeric(sub(pattern, "\\1", df[[col]], perl = TRUE))
  })

  cbind(
    df[, c("POS", "ID", "REF", "ALT")],
    # check.names = FALSE so as.data.frame() does not rename the sample columns
    # itself (e.g. "1001" -> "X1001").
    as.data.frame(setNames(ds_vals, out_names), stringsAsFactors = FALSE, check.names = FALSE),
    stringsAsFactors = FALSE
  )
}
