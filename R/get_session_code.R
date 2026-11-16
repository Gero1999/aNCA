#' Generate a session script code in R that can replicate the App outputs
#'
#' @param template_path Path to the R script template (e.g., script_template.R)
#' @param setts_obj The read settings object with all analysis specifications
#' or the full session object from the App.
#' @param output_path Path to write the resulting script file (e.g., "output_script.R")
#' @returns The output_path (invisibly)
#' @keywords internal
#' @noRd
get_code <- function(
  setts_obj,
  output_path,
  template_path = system.file("www/templates/script_template.R", package = "aNCA")
) {

  # Helper to get value from settings_list by path (e.g., 'settings$method')
  get_session_value <- function(path) {
    parts <- strsplit(path, "\\$")[[1]]
    obj <- setts_obj
    for (p in parts) {
      if (inherits(obj[[p]], "reactive")) {
        obj <- obj[[p]]()
      } else {
        obj <- obj[[p]]
      }
      if (is.null(obj)) {
        return(NULL)
      }
    }
    obj
  }
  # Read template
  script <- readLines(template_path, warn = FALSE) %>%
    paste(collapse = "\n")

  # Find all settings_list$...
  pattern <- "settings_list(\\$[a-zA-Z0-9_]+(\\(\\))?(\\$[a-zA-Z0-9_]+)*)"
  matches <- gregexpr(pattern, script, perl = TRUE)[[1]]
  if (matches[1] == -1) {
    stop(
      "Template has no placeholders (settings_list...) to substitute.",
      "This may be due to an incorrect file path, a missing template, ",
      "or a modified template without placeholders."
    )
  }

  # Replace each match with deparsed value
  for (i in rev(seq_along(matches))) {
    start <- matches[i]
    len <- attr(matches, "match.length")[i]
    matched <- substr(script, start, start + len - 1)
    # Extract the path after settings_list$
    path <- sub("^settings_list\\$", "", matched)
    value <- get_session_value(path)

    deparsed <- clean_deparse(value)
    script <- paste0(
      substr(script, 1, start - 1),
      deparsed,
      substr(script, start + len, nchar(script))
    )
  }

  # Split back into lines
  script_lines <- strsplit(script, "\n")[[1]]
  writeLines(script_lines, output_path)
  invisible(output_path)
}

#' Convert settings into faithful, readable R code (internal)
#'
#' Base deparsing preserves literal escaping, names, classes, typed missing
#' values, and floating-point precision. Tibbles are represented as data frames
#' because generated analysis scripts do not require their display classes.
#'
#' @param obj An R object to serialize into a reconstructing expression.
#' @returns One character string containing runnable R code.
#' @keywords internal
#' @noRd
clean_deparse <- function(obj) {
  if (inherits(obj, "tbl_df")) obj <- as.data.frame(obj)
  # Avoid niceNames: it can emit unescaped backslashes inside argument names.
  control <- c("keepNA", "keepInteger", "showAttributes", "quoteExpressions", "digits17")
  paste(deparse(obj, width.cutoff = 100, backtick = TRUE, control = control), collapse = "\n")
}

#' Generate a session script from a settings YAML file
#'
#' Reads a settings YAML file and generates an R script that can reproduce the
#' session. Mapping, filters, ratio table, and units are all read from the YAML.
#'
#' @param settings_file_path Path to the YAML file containing the settings.
#' @param data_path Path to the data file to be referenced in the script.
#' @param template_path Path to the R script template file. By default, uses
#'   the one installed from your aNCA package version.
#' @param output_path Path to write the resulting script file.
#'
#' @returns Invisibly returns the output_path.
#' @export
get_settings_code <- function(
  settings_file_path,
  data_path,
  output_path = "settings_code.R",
  template_path = system.file("www/templates/script_template.R", package = "aNCA")
) {
  settings <- read_settings(settings_file_path)

  session <- list(
    settings = settings[["settings"]],
    slope_rules = settings[["slope_rules"]],
    data_path = data_path,
    mapping = settings[["mapping"]],
    applied_filters = settings[["filters"]],
    ratio_table = settings[["settings"]][["ratio_table"]],
    units_table = settings[["settings"]][["units"]],
    extra_vars_to_keep = c(
      settings[["mapping"]][["Grouping_Variables"]], "DOSEA", "ATPTREF", "ROUTE"
    ),
    time_duplicate_rows = settings[["time_duplicate_keys"]]
  )

  get_code(
    template_path = template_path,
    setts_obj = session,
    output_path = output_path
  )
  invisible(output_path)
}

#' Generate a session script from a Shiny session object
#'
#' This function generates an R script that can reproduce the outputs of a Shiny app session.
#' It extracts the session's user data, substitutes it into a script template, and writes
#' the result to a file.
#'
#' @param session The Shiny session object containing user data (typically from a running app).
#' @param output_path Path to write the resulting script file (e.g., "output_script.R").
#' @param template_path Path to the R script template file. By default, uses the one
#' installed from your aNCA package version.
#' @returns Invisibly returns the output_path.
#' @keywords Internal
#' @noRd
get_session_code <- function(
  session,
  output_path,
  template_path = system.file("www/templates/script_template.R", package = "aNCA")
) {
  get_code(
    template_path = template_path,
    setts_obj = session$userData,
    output_path = output_path
  )
  invisible(output_path)
}
