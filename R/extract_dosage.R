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
#' @param df A data.frame as returned by \code{tabixr::query_vcf_positions()}.
#'
#' @return A data.frame with columns \code{POS} (integer), \code{ID},
#'   \code{REF}, \code{ALT} (all character), and one numeric column per sample.
#' @export
extract_dosage <- function(df) {
  if (nrow(df) == 0L) {
    fmt_col   <- match("FORMAT", names(df))
    samp_cols <- names(df)[(fmt_col + 1L):ncol(df)]
    out <- df[, c("POS", "ID", "REF", "ALT")]
    for (col in samp_cols) out[[col]] <- numeric(0)
    return(out)
  }

  fmt      <- df$FORMAT[1L]
  fields   <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
  ds_idx   <- match("DS", fields)
  if (is.na(ds_idx))
    stop("DS field not found in FORMAT: ", fmt)

  fmt_col   <- match("FORMAT", names(df))
  samp_cols <- names(df)[(fmt_col + 1L):ncol(df)]

  # Extract the ds_idx-th colon-delimited token from each sample column.
  pattern <- paste0("^(?:[^:]*:){", ds_idx - 1L, "}([^:]*).*$")
  ds_vals <- lapply(samp_cols, function(col) {
    as.numeric(sub(pattern, "\\1", df[[col]], perl = TRUE))
  })

  cbind(
    df[, c("POS", "ID", "REF", "ALT")],
    as.data.frame(setNames(ds_vals, samp_cols), stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )
}
