# Production modules share one bounded registry of edit observers.
local({
  library(shiny)
  library(reactable)
  library(reactable.extras)
  library(logger)
  shiny_dir <- system.file("shiny", package = "aNCA")
  source(file.path(shiny_dir, "functions", "utils.R"), local = TRUE)
  source(file.path(shiny_dir, "functions", "utils-table_edits.R"), local = TRUE)
  source(file.path(shiny_dir, "modules", "tab_nca", "setup", "manual_slopes_table.R"),
         local = TRUE)
}, envir = parent.env(environment()))

describe("table edit validation", {
  it("ignores stale rows and undeclared columns without adding cells", {
    original <- data.frame(value = "A", other = "B")
    invalid <- list(
      list(row = 0, column = "value", value = "X"),
      list(row = 2, column = "value", value = "X"),
      list(row = NA_real_, column = "value", value = "X"),
      list(row = 1.5, column = "value", value = "X"),
      list(row = c(1, 2), column = "value", value = "X"),
      list(row = 1, column = "other", value = "X"),
      list(row = 1, column = "value", value = c("X", "Y")),
      list(row = 1, column = "value", value = NULL)
    )
    for (edit in invalid) expect_identical(.apply_table_edit(original, edit, "value"), original)
    expect_identical(.apply_table_edit(original, NULL, "value"), original)
    updated <- .apply_table_edit(original, list(row = 1, column = "value", value = ""), "value")
    expect_identical(updated$value, "")
    expect_identical(updated$other, "B")
    numeric_edit <- .apply_table_edit(original, list(row = 1, column = 1, value = "N"), "value")
    expect_equal(numeric_edit$value, "N")
    expect_identical(.apply_table_edit(original, list(row = 1, column = 2, value = "X"),
                                       "value"), original)
  })
})

describe("editable column observer lifecycle", {
  edit_test_server <- function(id) {
    moduleServer(id, function(input, output, session) {
      current <- reactiveVal(data.frame(value = "A", other = "B"))
      writes <- reactiveVal(0L)
      table <- function(value) {
        if (missing(value)) return(current())
        writes(writes() + 1L)
        current(value)
      }
      handlers <- .bind_table_edits(input, table, reactive(names(current())))
      list(current = current, handlers = handlers)
    })
  }

  it("keeps handlers stable and performs one write per changed edit", {
    registry <- NULL
    testServer(edit_test_server, {
      session$flushReact()
      registry <<- handlers
      originals <- handlers$observers
      for (i in seq_len(50)) {
        session$setInputs(value = list(row = 1, column = "value", value = paste0("edit", i)))
      }
      expect_identical(handlers$observers, originals)
      expect_equal(writes(), 50L)
      expect_equal(current()$value, "edit50")
      expect_equal(current()$other, "B")
    })
    expect_length(registry$observers, 0)
  })

  it("retires removed columns and does not replay browser values when re-added", {
    testServer(edit_test_server, {
      session$flushReact()
      removed <- handlers$observers$value
      retained <- handlers$observers$other
      session$setInputs(value = list(row = 1, column = "value", value = "OLD"))
      current(data.frame(other = "B"))
      session$flushReact()
      expect_named(handlers$observers, "other")
      expect_identical(handlers$observers$other, retained)
      current(data.frame(other = "B", value = "NEW"))
      session$flushReact()
      expect_false(identical(handlers$observers$value, removed))
      expect_equal(current()$value, "NEW")
      session$setInputs(value = list(row = 1, column = "value", value = "CURRENT"))
      expect_equal(current()$value, "CURRENT")
      expect_equal(writes(), 2L)
    })
  })
})

describe("manual slope editing handlers", {
  it("keeps one handler per column while rules are edited and removed", {
    original <- PKNCA_create_data_object(
      subset(adnca_example, USUBJID == "S1-01" & PARAM == "DrugA")
    )
    original$intervals <- unique(original$conc$data[group_vars(original)])
    rule <- original$intervals[1, , drop = FALSE]
    rule$TYPE <- "Exclusion"
    rule$RANGE <- "4.48"
    rule$REASON <- "Initial"
    testServer(
      manual_slopes_table_server,
      args = list(pknca_data = reactive(original), manual_slopes_override = reactive(rule)),
      {
        session$flushReact()
        initial <- edit_observers$observers
        expect_length(initial, ncol(rule))
        for (i in seq_len(25)) {
          session$setInputs(edit_REASON = list(row = 1, column = "REASON", value = paste0("R", i)))
        }
        expect_equal(manual_slopes()$REASON, "R25")
        expect_identical(edit_observers$observers, initial)
        session$setInputs(edit_REASON = list(row = 2, column = "REASON", value = "stale"))
        expect_equal(nrow(manual_slopes()), 1L)
        expect_equal(manual_slopes()$REASON, "R25")
        manual_slopes(NULL)
        session$flushReact()
        expect_length(edit_observers$observers, 0)
      }
    )
  })
})
