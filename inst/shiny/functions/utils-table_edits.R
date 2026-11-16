#' Check a browser-supplied row/column index before indexing a data frame.
.is_table_edit_index <- function(index, size) {
  if (!is.numeric(index) || length(index) != 1L) return(FALSE)
  isTRUE(is.finite(index) & index == trunc(index) & index >= 1 & index <= size)
}

#' Resolve named and legacy positional column references from a table editor.
.table_edit_column <- function(table, column) {
  if (.is_table_edit_index(column, ncol(table))) {
    names(table)[column]
  } else {
    column
  }
}

#' Apply a valid edit from one declared table column.
#' Invalid or stale browser events leave the table unchanged.
.apply_table_edit <- function(table, edit, column) {
  if (!is.data.frame(table) || !is.list(edit)) return(table)
  row <- edit$row
  if (!.is_table_edit_index(row, nrow(table))) return(table)
  if (!identical(.table_edit_column(table, edit$column), column)) return(table)
  if (!column %in% names(table)) return(table)
  if (!is.atomic(edit$value) || length(edit$value) != 1L) return(table)
  table[row, column] <- edit$value
  table
}

#' Bind one edit observer per current column and retire removed-column observers.
#' Cell changes do not register new observers. Retained browser input values are
#' not replayed when a column is recreated after a schema change.
.bind_table_edits <- function(input, table, columns, prefix = "",
                              session = shiny::getDefaultReactiveDomain()) {
  registry <- new.env(parent = emptyenv())
  registry$observers <- list()
  add_column <- function(column) {
    input_id <- paste0(prefix, column)
    initial_edit <- isolate(input[[input_id]])
    first_event <- TRUE
    observeEvent(input[[input_id]], {
      edit <- input[[input_id]]
      if (first_event) {
        first_event <<- FALSE
        if (identical(edit, initial_edit)) return()
      }
      updated <- .apply_table_edit(table(), edit, column)
      if (!identical(updated, table())) table(updated)
    })
  }
  sync_columns <- function(current) {
    for (column in setdiff(names(registry$observers), current)) {
      registry$observers[[column]]$destroy()
      registry$observers[[column]] <- NULL
    }
    for (column in setdiff(current, names(registry$observers))) {
      registry$observers[[column]] <- add_column(column)
    }
  }
  sync_columns(isolate(columns()))
  observeEvent(columns(), sync_columns(columns()), ignoreNULL = FALSE)
  session$onSessionEnded(function() {
    lapply(registry$observers, function(observer) observer$destroy())
    registry$observers <- list()
  })
  invisible(registry)
}
