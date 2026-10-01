#' Tabulate PRS Values: One Row per Subject, One Column per Model
#'
#' Turns the list returned by \code{\link{compute_prs}} into a data.frame with a
#' column of subject IDs and one column of PRS values per model, named after the
#' model.
#'
#' Scores are matched to subjects by ID, not by position. All models from one
#' \code{compute_prs()} call cover the same subjects, but if \code{results} was
#' assembled from several runs (e.g. \code{c(run1, run2)}) whose subject sets
#' differ, every subject appearing in any model is kept, with \code{NA} for the
#' models that did not score it, and a warning is given. Subjects excluded from
#' a run because they were missing from a genotype file are not in it; they are
#' listed in \code{results[[1]]$excluded_samples}.
#'
#' Column names are the model names exactly as given (they are not altered to
#' be syntactically valid R names).
#'
#' @param results   A named list as returned by \code{compute_prs()} (or read
#'   back from an \code{.rds} file saved from it): one element per model, each
#'   a list with a named numeric vector \code{prs}. The list names become the
#'   column names.
#' @param id_column Name of the column holding the subject IDs. Default
#'   \code{"subject_id"}.
#'
#' @return A data.frame with one row per subject: the ID column first, then one
#'   numeric column per model, in the order of \code{results}.
#'
#' @examples
#' \dontrun{
#' results <- compute_prs("path/to/genotypes", pgs_dir = "path/to/models")
#' prs_df  <- prs_table(results)
#' head(prs_df)
#'
#' # or from a saved run
#' prs_df <- prs_table(readRDS("results.rds"))
#' }
#' @export
prs_table <- function(results, id_column = "subject_id") {
  if (!is.list(results) || length(results) == 0L)
    stop("results must be a non-empty list as returned by compute_prs()", call. = FALSE)
  if (is.null(names(results)) || any(is.na(names(results)) | !nzchar(names(results))))
    stop("every element of results must be named (the model ID); the names become the column names",
         call. = FALSE)
  dup <- unique(names(results)[duplicated(names(results))])
  if (length(dup) > 0L)
    stop(sprintf("duplicate model name(s) in results: %s", paste(dup, collapse = ", ")), call. = FALSE)
  if (!is.character(id_column) || length(id_column) != 1L || is.na(id_column) || !nzchar(id_column))
    stop("id_column must be a single non-empty string", call. = FALSE)
  if (id_column %in% names(results))
    stop(sprintf("a model is named '%s', the same as id_column; choose a different id_column", id_column),
         call. = FALSE)

  scores <- lapply(names(results), function(nm) {
    prs <- results[[nm]]$prs
    if (!is.numeric(prs) || is.null(names(prs)))
      stop(sprintf("results[['%s']]$prs must be a named numeric vector of scores (as returned by compute_prs())", nm),
           call. = FALSE)
    prs
  })
  names(scores) <- names(results)

  ids <- Reduce(union, lapply(scores, names))
  mat <- do.call(cbind, lapply(scores, function(p) unname(p[ids])))
  colnames(mat) <- names(scores)

  n_missing <- sum(is.na(mat)) - sum(vapply(scores, function(p) sum(is.na(p)), numeric(1L)))
  if (n_missing > 0L)
    warning(sprintf("models cover different subjects: %d score(s) are NA because a model did not score that subject",
                    n_missing), call. = FALSE)

  out <- data.frame(ids, mat, check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  names(out)[1L] <- id_column
  out
}
