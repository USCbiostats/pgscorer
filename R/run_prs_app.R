#' Explore PRS Results in a Shiny App
#'
#' Opens an interactive app for results saved from \code{\link{compute_prs}}
#' (an \code{.rds} file holding the list it returns). It has three tabs:
#' \describe{
#'   \item{Summary}{A statement of what was done, with a trait name you can edit
#'     (it defaults to the name of the results file), the models fitted with
#'     their SNP counts, and a table of SNPs by chromosome for each model with
#'     the number not included in the PRS in dark red.}
#'   \item{Histogram}{Pick one model and see the distribution of its PRS.}
#'   \item{Scatter}{Pick two models and see their PRS plotted against each
#'     other, with the R-squared between them.}
#' }
#' This suits results with many models, where a single static report becomes
#' unwieldy. The SNP counts need results from the current \code{compute_prs()}
#' (see \code{n_snps}); without them the Summary tab omits those tables but the
#' other two tabs work as usual.
#'
#' @param file Optional path to a results \code{.rds} file to open at start-up.
#'   When \code{NULL} (default), choose a file in the app.
#' @param launch.browser Passed to \code{shiny::runApp()}; whether to open the
#'   app in a browser. Default \code{TRUE}.
#'
#' @return Called for its side effect of running the app; returns when the app
#'   is stopped.
#'
#' @examples
#' \dontrun{
#' run_prs_app()
#' run_prs_app("results.rds")
#' }
#' @export
run_prs_app <- function(file = NULL, launch.browser = TRUE) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("the shiny package is required: install.packages(\"shiny\")", call. = FALSE)
  if (!is.null(file) && !file.exists(file))
    stop(sprintf("file not found: %s", file), call. = FALSE)
  # Results files can be large (one PRS vector per model); lift Shiny's 5 MB upload cap.
  old <- options(shiny.maxRequestSize = 4 * 1024^3)
  on.exit(options(old), add = TRUE)
  shiny::runApp(.prs_app(file), launch.browser = launch.browser)
}

# Check that an object read from a file looks like compute_prs() results; return
# it, or stop with a message the app can show.
.check_results <- function(x) {
  if (!is.list(x) || length(x) == 0L || is.null(names(x)) ||
      any(is.na(names(x)) | !nzchar(names(x))))
    stop("not a results list: expected a named list with one element per model", call. = FALSE)
  ok <- vapply(x, function(r) is.list(r) && is.numeric(r$prs), logical(1L))
  if (!all(ok))
    stop("not a results list: every model needs a numeric `prs` element", call. = FALSE)
  x
}

# TRUE when every model carries the SNP counts the summary tables need.
.has_snp_counts <- function(results) {
  all(vapply(results, function(r)
    !is.null(r$n_snps) && !is.null(r$n_snps_by_chr) &&
      !is.null(r$unmatched_by_chr) && !is.null(r$excluded_chr_counts), logical(1L)))
}

# The Summary tab's body (HTML) for `results` and the chosen trait name.
.summary_html <- function(results, trait) {
  esc    <- .html_esc
  fmt    <- .fmt_int
  models <- names(results)
  n_subj <- length(results[[1L]]$prs)
  excl   <- length(results[[1L]]$excluded_samples)

  out <- c(
    sprintf(paste0("<p>The following is a summary of the results of fitting the PGS catalog ",
                   "models for <strong>%s</strong>.</p>"), esc(trait)),
    sprintf("<p>%d model%s fitted; PRS computed for %s subjects%s.</p>",
            length(models), if (length(models) == 1L) "" else "s", fmt(n_subj),
            if (excl > 0L) sprintf(" (%s left out for missing genotype data)", fmt(excl)) else ""))

  if (!.has_snp_counts(results))
    return(paste(c(out,
      paste0("<p class=\"note\">SNP counts are not available: these results were saved by an older ",
             "version of <code>compute_prs()</code>. Add them with <code>add_snp_counts.R</code> ",
             "or re-run the models.</p>"),
      "<h4>Models fitted</h4>",
      sprintf("<p class=\"models\">%s</p>", paste(esc(models), collapse = ", "))), collapse = "\n"))

  n_snps  <- vapply(results, function(r) as.numeric(r$n_snps), numeric(1L))
  skipped <- vapply(results, function(r)
    sum(r$unmatched_by_chr) + sum(r$excluded_chr_counts), numeric(1L))
  pct <- ifelse(n_snps > 0, 100 * skipped / n_snps, 0)
  pct <- ifelse(skipped > 0 & pct < 0.1, "&lt;0.1%", sprintf("%.1f%%", pct))

  models_table <- c(
    "<table class=\"prs\">",
    "<thead><tr><th>Model</th><th>SNPs</th><th>Not included</th></tr></thead>",
    "<tbody>",
    sprintf("<tr><th scope=\"row\">%s</th><td>%s</td><td><span class=\"skip\">%s</span> (%s)</td></tr>",
            esc(models), fmt(n_snps), fmt(skipped), pct),
    "</tbody>",
    "</table>")

  paste(c(out,
    "<h4>Models fitted</h4>", models_table,
    "<h4>SNPs by chromosome</h4>",
    "<div class=\"scroll\">", .chr_snp_table_html(results), "</div>"), collapse = "\n")
}

.prs_app <- function(file = NULL) {
  css <- "
    table.prs, .scroll table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
    table.prs th, table.prs td, .scroll th, .scroll td { padding: 0.25rem 0.9rem; border-bottom: 1px solid #ddd; text-align: right; white-space: nowrap; }
    table.prs thead th, .scroll thead th { border-bottom: 2px solid #888; }
    table.prs tbody th, .scroll tbody th, .scroll tfoot th { text-align: left; }
    .scroll tfoot th, .scroll tfoot td { border-top: 2px solid #888; border-bottom: none; font-weight: 600; }
    .scroll { overflow-x: auto; max-width: 100%; }
    .scroll tbody th, .scroll tfoot th { position: sticky; left: 0; background: #fff; }
    .skip { color: #8b0000; font-size: 0.9em; }
    .note { font-size: 0.9rem; color: #555; }
    .models { font-family: ui-monospace, Consolas, monospace; }
    .stat { font-size: 1rem; font-variant-numeric: tabular-nums; margin: 0.5rem 0 1rem; }"

  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(css))),
    shiny::titlePanel("Polygenic risk score results"),
    shiny::fileInput("file", "Results file (.rds saved from compute_prs)", accept = ".rds"),
    shiny::tabsetPanel(
      shiny::tabPanel("Summary",
        shiny::br(),
        shiny::textInput("trait", "Trait name", value = ""),
        shiny::uiOutput("summary")),
      shiny::tabPanel("Histogram",
        shiny::br(),
        shiny::selectInput("hist_model", "Model", choices = NULL),
        shiny::plotOutput("hist", height = "420px"),
        shiny::textOutput("hist_stats", container = function(...) shiny::div(class = "stat", ...))),
      shiny::tabPanel("Scatter",
        shiny::br(),
        shiny::fluidRow(
          shiny::column(3, shiny::selectInput("x_model", "Model on x-axis", choices = NULL)),
          shiny::column(3, shiny::selectInput("y_model", "Model on y-axis", choices = NULL))),
        shiny::textOutput("r2", container = function(...) shiny::div(class = "stat", ...)),
        shiny::plotOutput("scatter", height = "560px")))
  )

  server <- function(input, output, session) {
    # The source of the results: an uploaded file, else the file given at start-up.
    src <- shiny::reactive({
      if (!is.null(input$file))
        list(path = input$file$datapath, name = input$file$name)
      else if (!is.null(file))
        list(path = file, name = basename(file))
      else NULL
    })

    results <- shiny::reactive({
      s <- shiny::req(src())
      x <- tryCatch(.check_results(readRDS(s$path)),
                    error = function(e) conditionMessage(e))
      shiny::validate(shiny::need(is.list(x), paste0("Could not use ", s$name, ": ", x)))
      x
    })

    # New file: the trait name defaults to the file's name, and the model pickers refill.
    shiny::observeEvent(src(), {
      shiny::updateTextInput(session, "trait", value = tools::file_path_sans_ext(src()$name))
    })
    shiny::observeEvent(results(), {
      m <- names(results())
      shiny::updateSelectInput(session, "hist_model", choices = m, selected = m[1L])
      shiny::updateSelectInput(session, "x_model", choices = m, selected = m[1L])
      shiny::updateSelectInput(session, "y_model", choices = m, selected = m[min(2L, length(m))])
    })

    prompt <- "Choose a results file above to begin."

    output$summary <- shiny::renderUI({
      shiny::validate(shiny::need(!is.null(src()), prompt))
      # Blank (or not yet reported by the browser) falls back to the file's name.
      trait <- if (is.null(input$trait)) "" else trimws(input$trait)
      if (!nzchar(trait)) trait <- tools::file_path_sans_ext(src()$name)
      shiny::HTML(.summary_html(results(), trait))
    })

    output$hist <- shiny::renderPlot({
      shiny::validate(shiny::need(!is.null(src()), prompt))
      m <- shiny::req(input$hist_model)
      shiny::req(m %in% names(results()))
      .prs_histogram_plot(results()[[m]]$prs, m)
    }, res = 110)

    output$hist_stats <- shiny::renderText({
      m <- shiny::req(input$hist_model)
      shiny::req(m %in% names(results()))
      p <- results()[[m]]$prs
      sprintf("n = %s   mean = %s   SD = %s", .fmt_int(length(p)),
              signif(mean(p), 3), signif(stats::sd(p), 3))
    })

    # The two selected models' PRS over the subjects scored by both.
    pair <- shiny::reactive({
      shiny::validate(shiny::need(!is.null(src()), prompt))
      a <- shiny::req(input$x_model); b <- shiny::req(input$y_model)
      shiny::req(a %in% names(results()), b %in% names(results()))
      shiny::validate(shiny::need(a != b, "Choose two different models."))
      list(a = a, b = b, d = .prs_pair(results(), a, b))
    })

    output$r2 <- shiny::renderText({
      p <- pair()
      sprintf("R\u00b2 = %.3f   (%s subjects)", stats::cor(p$d$x, p$d$y)^2, .fmt_int(length(p$d$x)))
    })

    output$scatter <- shiny::renderPlot({
      p <- pair()
      .prs_scatter_plot(p$d$x, p$d$y, p$a, p$b)
    }, res = 110)
  }

  shiny::shinyApp(ui, server)
}
