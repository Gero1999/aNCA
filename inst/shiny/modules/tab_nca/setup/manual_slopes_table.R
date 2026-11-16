#' Manual Slopes Table UI for Slope Selection
#'
#' UI module for displaying and editing the manual slopes table (inclusion/exclusion rules).
#' Provides buttons to add/remove rules and a reactable table for editing.
#'
#' @param id Shiny module id
#' @param help_widget Optional Shiny UI element to display in the header row (e.g., help dropdown)
#' @return Shiny UI element (tagList of fluidRows)
manual_slopes_table_ui <- function(id, help_widget = NULL) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(
        width = 10,
        div(
          style = "display: flex; gap: 8px; align-items: center;",
          actionButton(ns("add_rule"), "+ Exclusion/Selection", class = "btn-success"),
          actionButton(ns("remove_rule"), "- Remove selected rows", class = "btn-warning")
        )
      ),
      if (!is.null(help_widget)) column(width = 2, help_widget)
    ),
    fluidRow(
      reactableOutput(ns("manual_slopes"))
    )
  )
}

# Validate imported profile keys without letting malformed settings disconnect the session.
.manual_slopes_match_data <- function(rules, data) {
  match_cols <- setdiff(names(rules), c("TYPE", "RANGE", "REASON"))
  valid <- is.data.frame(rules) && length(match_cols) > 0 &&
    all(match_cols %in% names(data$conc$data))
  if (!valid) return(FALSE)
  tryCatch(
    nrow(semi_join(rules, data$conc$data, by = match_cols)) == nrow(rules),
    error = function(e) FALSE
  )
}

#' Manual Slopes Table Server for Slope Selection
#'
#' Server module for managing the manual slopes table (inclusion/exclusion rules).
#' Handles adding/removing/editing rules, table's reactivity, and optional override logic.
#'
#' @param id Shiny module id
#' @param pknca_data Reactive providing the current PKNCA data object
#' @param manual_slopes_override Optional reactive providing a table to override manual slopes
#' @return List with:
#'   - manual_slopes: reactiveVal containing the current manual slopes table
#'   - refresh_reactable: reactiveVal for triggering table re-render
manual_slopes_table_server <- function(
  id, pknca_data, manual_slopes_override = NULL
) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    # Get group columns for the current PKNCA data (for table structure)
    slopes_pknca_groups <- reactive({
      req(pknca_data())
      pknca_data()$intervals %>%
        select(any_of(c(group_vars(pknca_data()))))
    })

    # manual_slopes: stores the current table of user rules (inclusion/exclusion).
    # Stays NULL until the user adds a rule or a settings override is applied.
    manual_slopes <- reactiveVal(NULL)

    refresh_reactable <- reactiveVal(0)
    pending_override <- reactiveVal(NULL)
    warned_override <- FALSE
    if (is.null(manual_slopes_override)) manual_slopes_override <- reactive(NULL)

    # An imported value is an explicit restore, independent of ordinary data refreshes.
    # Keep it pending until matching data are available, and clear the previous rules now.
    observeEvent(manual_slopes_override(), {
      pending_override(list(rules = manual_slopes_override()))
      warned_override <<- FALSE
      manual_slopes(NULL)
      refresh_reactable(refresh_reactable() + 1)
    }, ignoreNULL = FALSE, priority = 100)

    observe({
      request <- pending_override()
      req(!is.null(request))
      rules <- request$rules
      if (is.null(rules) || (is.data.frame(rules) && nrow(rules) == 0)) {
        pending_override(NULL)
        return()
      }
      req(pknca_data(), slopes_pknca_groups())
      valid <- .manual_slopes_match_data(rules, pknca_data())
      if (!valid) {
        if (!warned_override) {
          msg <- paste("Imported manual slopes do not match the current data;",
                       "previous rules were cleared.")
          log_warn(msg)
          showNotification(msg, type = "warning", duration = 5)
          warned_override <<- TRUE
        }
        return()
      }
      log_debug_list("Manual slopes override:", rules)
      pending_override(NULL)
      manual_slopes(rules)
      refresh_reactable(refresh_reactable() + 1)
    })

    # A user edit takes precedence over an unmatched pending import.
    observeEvent(manual_slopes(), {
      if (!is.null(pending_override())) pending_override(NULL)
    }, ignoreInit = TRUE, priority = -10)

    # Add a new row to the table when the user clicks the add button
    observeEvent(input$add_rule, {
      log_trace("{id}: adding manual slopes row")
      first_group <- slopes_pknca_groups()[1, ]
      time_col <- pknca_data()$conc$columns$time
      new_row <- cbind(
        first_group,
        data.frame(
          TYPE = "Exclusion",
          RANGE = paste0(
            inner_join(
              first_group,
              pknca_data()$conc$data,
              by = intersect(
                names(first_group),
                names(pknca_data()$conc$data)
              )
            )[[time_col]][2]
          ),
          REASON = ""
        )
      )

      old_rows <- if (!is.null(manual_slopes())) {
        filter(manual_slopes(), !is.na(TYPE))
      } else {
        new_row[0, ] # empty DF with matching columns
      }
      updated_data <- as.data.frame(
        rbind(old_rows, new_row),
        stringsAsFactors = FALSE
      )
      manual_slopes(updated_data)
      reset_reactable_memory()
      refresh_reactable(refresh_reactable() + 1)
    })

    # Remove selected rows from the table when the user clicks the remove button
    observeEvent(input$remove_rule, {
      log_trace("{id}: removing manual slopes row")
      req(manual_slopes())
      selected <- getReactableState("manual_slopes", "selected")
      req(selected)
      edited_slopes <- manual_slopes()[-selected, ]
      if (nrow(edited_slopes) == 0) edited_slopes <- NULL
      manual_slopes(edited_slopes)
      reset_reactable_memory()
      refresh_reactable(refresh_reactable() + 1)
    })

    # Render the manual slopes table (reactable)
    output$manual_slopes <- renderReactable({
      req(manual_slopes())
      log_trace("{id}: rendering slope edit data table")
      isolate({
        data <- manual_slopes()
      })
      # Define columns: group columns (dynamic), then TYPE/RANGE/REASON (fixed)
      fixed_columns <- list(
        TYPE = colDef(
          cell = dropdown_extra(
            id = ns("edit_TYPE"),
            choices = c("Selection", "Exclusion"),
            class = "dropdown-extra"
          ),
          width = 200
        ),
        RANGE = colDef(
          cell = text_extra(
            id = ns("edit_RANGE")
          )
        ),
        REASON = colDef(
          cell = text_extra(
            id = ns("edit_REASON")
          ),
          width = 400
        )
      )
      dynamic_columns <- lapply(colnames(slopes_pknca_groups()), function(col) {
        colDef(
          cell = dropdown_extra(
            id = ns(paste0("edit_", col)),
            choices = unique(slopes_pknca_groups()[[col]]),
            class = "dropdown-extra"
          ),
          width = 150
        )
      })
      names(dynamic_columns) <- colnames(slopes_pknca_groups())
      all_columns <- c(dynamic_columns, fixed_columns)
      reactable(
        data = data,
        defaultColDef = colDef(align = "center"),
        columns = all_columns,
        selection = "multiple",
        defaultExpanded = TRUE,
        borderless = TRUE,
        theme = reactableTheme(
          rowSelectedStyle = list(
            backgroundColor = "#eee",
            boxShadow = "inset 2px 0 0 0 #ffa62d"
          )
        )
      )
    }) %>%
      shiny::bindEvent(refresh_reactable())

    # Reconcile handlers only when the editable schema changes, not on cell edits.
    edit_observers <- .bind_table_edits(
      input, manual_slopes, reactive(names(manual_slopes())), prefix = "edit_",
      session = session
    )

    # Output: manual_slopes (reactiveVal) and refresh_reactable (for UI updates)
    list(
      manual_slopes = manual_slopes,
      refresh_reactable = refresh_reactable
    )
  })
}
