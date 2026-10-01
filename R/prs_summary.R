#' Summarize PRS Results as an HTML Report
#'
#' Writes an HTML file summarizing the results returned by
#' \code{\link{compute_prs}}. The report opens with a statement of what was
#' done, followed by the list of PGS Catalog models that were fitted (with the
#' number of SNPs in each) and a table of SNPs by chromosome for each model.
#'
#' The results must contain \code{n_snps} and \code{n_snps_by_chr}, as returned
#' by the current \code{compute_prs()}.
#'
#' @param results   A named list as returned by \code{compute_prs()} (or read
#'   back from an \code{.rds} file saved from it): one element per model. The
#'   list names are the model IDs.
#' @param trait     Character; the trait the models were fitted for.
#' @param data_name Character; a name for the data set the models were fitted
#'   on.
#' @param file      Path of the HTML file to write. Default
#'   \code{"prs_summary.html"}.
#'
#' @return The path of the HTML file, invisibly.
#'
#' @examples
#' \dontrun{
#' results <- readRDS("results.rds")
#' prs_summary(results, trait = "Trait 1", data_name = "Sample",
#'             file = "prs_summary.html")
#' }
#' @export
prs_summary <- function(results, trait, data_name, file = "prs_summary.html") {
  if (!is.list(results) || length(results) == 0L)
    stop("results must be a non-empty list as returned by compute_prs()", call. = FALSE)
  if (is.null(names(results)) || any(is.na(names(results)) | !nzchar(names(results))))
    stop("every element of results must be named (the model ID)", call. = FALSE)
  for (arg in c("trait", "data_name", "file")) {
    v <- get(arg)
    if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v))
      stop(sprintf("%s must be a single non-empty string", arg), call. = FALSE)
  }

  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }

  models <- names(results)
  if (!all(vapply(results, function(r) !is.null(r$n_snps) && !is.null(r$n_snps_by_chr), logical(1L))))
    stop("results lack n_snps / n_snps_by_chr (saved by an older compute_prs()); ",
         "add them with add_snp_counts.R or re-run compute_prs()", call. = FALSE)

  n_snps <- vapply(results, function(r) as.numeric(r$n_snps), numeric(1L))
  fmt    <- function(x) formatC(x, format = "d", big.mark = ",")

  # SNPs by chromosome: one row per chromosome, one column per model. The
  # chromosome set is the union across models, in the order first seen.
  chroms <- unique(unlist(lapply(results, function(r) names(r$n_snps_by_chr)), use.names = FALSE))
  by_chr <- vapply(results, function(r) {
    v <- r$n_snps_by_chr[chroms]
    v[is.na(v)] <- 0L
    as.numeric(v)
  }, numeric(length(chroms)))
  by_chr <- matrix(by_chr, nrow = length(chroms), dimnames = list(chroms, models))

  # SNPs not used in the score, same layout: those with no genotype record
  # (unmatched) plus those whose effect allele matched neither REF nor ALT.
  skipped <- vapply(results, function(r) {
    v <- as.numeric(r$unmatched_by_chr[chroms]) + as.numeric(r$excluded_chr_counts[chroms])
    v[is.na(v)] <- 0
    v
  }, numeric(length(chroms)))
  skipped <- matrix(skipped, nrow = length(chroms), dimnames = list(chroms, models))
  # Percent of SNPs not included; "<0.1%" rather than "0.0%" when some, but very few, are.
  pct <- function(n, s) {
    p <- ifelse(n > 0, 100 * s / n, 0)
    ifelse(s > 0 & p < 0.1, "&lt;0.1%", sprintf("%.1f%%", p))
  }
  cell <- function(n, s)
    paste0("<td>", fmt(n), "<br><span class=\"skip\">", fmt(s), "</span> (", pct(n, s), ")</td>",
           collapse = "")

  chr_table <- c(
    "<table>",
    paste0("<thead><tr><th>Chromosome</th>",
           paste0("<th>", esc(models), "</th>", collapse = ""), "</tr></thead>"),
    "<tbody>",
    vapply(seq_along(chroms), function(i)
      paste0("<tr><th scope=\"row\">", esc(chroms[i]), "</th>",
             cell(by_chr[i, ], skipped[i, ]), "</tr>"),
      character(1L)),
    "</tbody>",
    paste0("<tfoot><tr><th scope=\"row\">Total</th>",
           cell(colSums(by_chr), colSums(skipped)), "</tr></tfoot>"),
    "</table>",
    paste0("<p class=\"note\">Black: SNPs in the model. <span class=\"skip\">Dark red</span>, with the percent of the total in parentheses: ",
           "SNPs not included in the PRS calculation (no matching genotype record, or an ",
           "effect allele matching neither REF nor ALT).</p>")
  )

  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    sprintf("<title>PRS Summary: %s</title>", esc(trait)),
    "<style>",
    "body { font-family: system-ui, sans-serif; max-width: 48rem; margin: 2rem auto; padding: 0 1rem; line-height: 1.5; color: #222; }",
    "h1 { font-size: 1.6rem; } h2 { font-size: 1.2rem; margin-top: 2rem; }",
    "li { font-family: ui-monospace, Consolas, monospace; }",
    "table { border-collapse: collapse; font-variant-numeric: tabular-nums; }",
    "th, td { padding: 0.25rem 0.9rem; border-bottom: 1px solid #ddd; text-align: right; }",
    "thead th { border-bottom: 2px solid #888; } tbody th, tfoot th { text-align: left; }",
    ".skip { color: #8b0000; font-size: 0.9em; } .note { font-size: 0.9rem; color: #555; }",
    "tfoot th, tfoot td { border-top: 2px solid #888; border-bottom: none; font-weight: 600; }",
    "</style>",
    "</head>",
    "<body>",
    "<h1>Polygenic Score Summary</h1>",
    sprintf(paste0("<p>The following is a summary of the results of fitting the PGS ",
                   "catalog models for <strong>%s</strong> on the <strong>%s</strong> data.</p>"),
            esc(trait), esc(data_name)),
    sprintf("<h2>Models fitted (%d)</h2>", length(models)),
    "<ul>",
    sprintf("<li>%s &mdash; %s SNPs</li>", esc(models), fmt(n_snps)),
    "</ul>",
    "<h2>SNPs by chromosome</h2>",
    chr_table,
    "</body>",
    "</html>"
  )
  writeLines(html, file, useBytes = TRUE)
  invisible(file)
}
