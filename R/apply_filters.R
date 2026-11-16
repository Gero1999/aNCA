#' Apply Filters to a Dataset
#'
#' This function applies a set of filters to a dataset. Each filter specifies a column,
#' condition, and value to filter the dataset.
#'
#' @param data A data frame containing the raw data to be filtered.
#' @param filters  A list of filters, where each filter is a list containing
#'                 the column, condition, and value.
#'
#' @returns A data frame containing the filtered data.
#'
#' @details
#' The function iterates over the list of filters and applies each filter to the dataset.
#' The supported conditions are `==`, `>`, `<`, `>=`, `<=`, `!=` and `min-max`.
#' Equality/exclusion accept sets of values. A range accepts two numeric bounds
#' or a signed range string such as `-2-2` or `1e-3-2e-3`. Invalid filters error.
#'
#' @examples
#' # Example usage:
#' data <- mtcars
#' filters <- list(
#'   list(column = "mpg", condition = ">", value = "20"),
#'   list(column = "cyl", condition = "==", value = "6")
#' )
#' filtered_data <- apply_filters(data, filters)
#'
#' @importFrom dplyr filter between
#' @importFrom rlang sym
#' @export
apply_filters <- function(data, filters) {
  for (filter in filters) {
    # check if filter is valid #
    if (is.null(filter)) next
    .validate_filter(filter, data)
    column <- filter$column
    value <- filter$value

    # Strip "(NA)" sentinel before numeric coercion
    na_result <- strip_na_sentinel(value)
    value <- na_result$value
    include_na <- na_result$include_na

    if (is.numeric(data[[column]]) && filter$condition != "min-max") {
      converted <- suppressWarnings(as.numeric(value))
      if (any(is.na(converted))) stop("Invalid numeric filter value for ", column, ".")
      value <- converted
    }

    data <- apply_single_filter(data, column, filter$condition, value, include_na)
  }

  data
}

# Validate before evaluation so malformed filters cannot silently change selection.
.validate_filter <- function(filter, data) {
  missing_fields <- setdiff(c("column", "condition", "value"), names(filter))
  if (length(missing_fields) > 0) {
    stop("Missing required filter fields: ", paste0(missing_fields, collapse = ", "))
  }
  column <- filter$column
  if (length(column) != 1L || is.na(column)) stop("Filter must specify one column.")
  if (!column %in% names(data)) stop("Data is missing filtered column: ", column)
  if (length(filter$condition) != 1L ||
        !filter$condition %in% c("==", "!=", ">", "<", ">=", "<=", "min-max")) {
    stop("Unsupported filter condition.")
  }
  if (filter$condition %in% c(">", "<", ">=", "<=") && length(filter$value) != 1L) {
    stop("Comparison filters require one value.")
  }
}

#' Strip the "(NA)" sentinel from a filter value vector
#'
#' @param value Character vector of filter values.
#'
#' @returns A list with `value` (cleaned vector) and `include_na` (logical).
#'
#' @keywords internal
strip_na_sentinel <- function(value) {
  include_na <- "(NA)" %in% value
  if (include_na) value <- setdiff(value, "(NA)")
  list(value = value, include_na = include_na)
}

#' Apply a single filter condition to a data frame
#'
#' @param data A data frame.
#' @param column Column name to filter on.
#' @param condition Filter condition string.
#' @param value Filter value(s).
#' @param include_na Logical; whether to include NA rows.
#'
#' @returns A filtered data frame.
#'
#' @keywords internal
apply_single_filter <- function(data, column, condition, value, include_na) {
  switch(
    condition,
    "==" = filter_eq(data, column, value, include_na),
    ">" = filter(data, !!sym(column) > value),
    "<" = filter(data, !!sym(column) < value),
    ">=" = filter(data, !!sym(column) >= value),
    "<=" = filter(data, !!sym(column) <= value),
    "!=" = filter_neq(data, column, value, include_na),
    "min-max" = filter_min_max(data, column, value)
  )
}

#' @keywords internal
filter_eq <- function(data, column, value, include_na) {
  if (include_na) {
    filter(data, !!sym(column) %in% value | is.na(!!sym(column)))
  } else {
    filter(data, !!sym(column) %in% value)
  }
}

#' @keywords internal
filter_neq <- function(data, column, value, include_na) {
  if (include_na) {
    filter(data, !(!!sym(column) %in% value) | is.na(!!sym(column)))
  } else {
    filter(data, !is.na(!!sym(column)) & !(!!sym(column) %in% value))
  }
}

#' @keywords internal
filter_min_max <- function(data, column, value) {
  number <- "[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?"
  pattern <- paste0("^\\s*(", number, ")\\s*-\\s*(", number, ")\\s*$")
  if (length(value) == 1L && !is.na(value)) {
    matched <- regmatches(value, regexec(pattern, value, perl = TRUE))[[1]]
    value <- if (length(matched) == 3L) matched[2:3] else NA_real_
  }
  min_max <- suppressWarnings(as.numeric(value))
  if (length(min_max) != 2L || any(!is.finite(min_max)) || min_max[1] > min_max[2]) {
    stop("Invalid min-max value for ", column, ": supply two finite, ascending bounds.")
  }

  filter(data, between(!!sym(column), min_max[1], min_max[2]))
}
