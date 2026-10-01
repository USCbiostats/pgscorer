#' Summarize PRS Results as an HTML Report
#'
#' Writes an HTML file summarizing the results returned by
#' \code{\link{compute_prs}}. The report opens with a statement of what was
#' done, followed by the list of PGS Catalog models that were fitted.
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
    "</style>",
    "</head>",
    "<body>",
    "<h1>Polygenic Score Summary</h1>",
    sprintf(paste0("<p>The following is a summary of the results of fitting the PGS ",
                   "catalog models for <strong>%s</strong> on the <strong>%s</strong> data.</p>"),
            esc(trait), esc(data_name)),
    sprintf("<h2>Models fitted (%d)</h2>", length(models)),
    "<ul>",
    sprintf("<li>%s</li>", esc(models)),
    "</ul>",
    "</body>",
    "</html>"
  )
  writeLines(html, file, useBytes = TRUE)
  invisible(file)
}
