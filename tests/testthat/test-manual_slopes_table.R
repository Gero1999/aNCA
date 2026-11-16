local({
  library(shiny)
  library(purrr)
  library(logger)
  library(reactable)
  library(reactable.extras)
  source(system.file("shiny/modules/tab_nca/setup/manual_slopes_table.R", package = "aNCA"),
         local = TRUE)
  reset_reactable_memory <- function(...) NULL
  log_debug_list <- function(...) NULL
}, envir = parent.env(environment()))

manual_slopes_fixture <- function() {
  pk <- PKNCA_create_data_object(subset(adnca_example, USUBJID == "S1-01" & PARAM == "DrugA"))
  pk$intervals <- unique(pk$conc$data[group_vars(pk)])
  rule <- pk$intervals[1, , drop = FALSE]
  rule$TYPE <- "Exclusion"
  rule$RANGE <- "4.48"
  rule$REASON <- "Session A"
  list(data = pk, rule = rule)
}

describe("manual slope rules during session restore", {
  it("replaces restored rules, preserves edits on data refresh, and clears empty imports", {
    fixture <- manual_slopes_fixture()
    data <- reactiveVal(fixture$data)
    override <- reactiveVal(fixture$rule)
    shiny::testServer(
      manual_slopes_table_server,
      args = list(pknca_data = data, manual_slopes_override = override),
      {
        session$flushReact()
        expect_identical(manual_slopes(), fixture$rule)

        edited <- manual_slopes()
        edited$REASON <- "Interactive edit"
        manual_slopes(edited)
        refreshed <- fixture$data
        refreshed$conc$data$AVAL <- refreshed$conc$data$AVAL * 2
        data(refreshed)
        session$flushReact()
        expect_identical(manual_slopes(), edited)

        second <- fixture$rule
        second$RANGE <- "6.011"
        second$REASON <- "Session B"
        override(second)
        session$flushReact()
        expect_identical(manual_slopes(), second)

        override(NULL)
        session$flushReact()
        expect_null(manual_slopes())
        override(fixture$rule)
        session$flushReact()
        override(fixture$rule[0, , drop = FALSE])
        session$flushReact()
        expect_null(manual_slopes())
        expect_false(session$isClosed())
      }
    )
  })

  it("waits for matching data while clearing the previous session's rules", {
    fixture <- manual_slopes_fixture()
    data <- reactiveVal(fixture$data)
    override <- reactiveVal(fixture$rule)
    shiny::testServer(
      manual_slopes_table_server,
      args = list(pknca_data = data, manual_slopes_override = override),
      {
        session$flushReact()
        second <- fixture$rule
        second$USUBJID <- "Session-B-subject"
        override(second)
        session$flushReact()
        expect_null(manual_slopes())

        replacement <- fixture$data
        replacement$conc$data$USUBJID <- second$USUBJID
        replacement$intervals$USUBJID <- second$USUBJID
        data(replacement)
        session$flushReact()
        expect_identical(manual_slopes(), second)
        expect_false(session$isClosed())
      }
    )
  })

  it("applies an import when the initial data arrive and repeats identical restores", {
    fixture <- manual_slopes_fixture()
    data <- reactiveVal(NULL)
    override <- reactiveVal(fixture$rule)
    shiny::testServer(
      manual_slopes_table_server,
      args = list(pknca_data = data, manual_slopes_override = override),
      {
        session$flushReact()
        expect_null(manual_slopes())
        data(fixture$data)
        session$flushReact()
        expect_identical(manual_slopes(), fixture$rule)

        edited <- fixture$rule
        edited$REASON <- "Edit after restore"
        manual_slopes(edited)
        session$flushReact()
        override(NULL)
        override(fixture$rule)
        session$flushReact()
        expect_identical(manual_slopes(), fixture$rule)
      }
    )
  })

  it("lets a user edit supersede an unmatched pending import", {
    fixture <- manual_slopes_fixture()
    data <- reactiveVal(fixture$data)
    override <- reactiveVal(fixture$rule)
    shiny::testServer(
      manual_slopes_table_server,
      args = list(pknca_data = data, manual_slopes_override = override),
      {
        session$flushReact()
        unmatched <- fixture$rule
        unmatched$USUBJID <- "Different-subject"
        override(unmatched)
        session$flushReact()
        expect_null(manual_slopes())
        edited <- fixture$rule
        edited$REASON <- "User choice"
        manual_slopes(edited)
        session$flushReact()
        replacement <- fixture$data
        replacement$conc$data$USUBJID <- unmatched$USUBJID
        replacement$intervals$USUBJID <- unmatched$USUBJID
        data(replacement)
        session$flushReact()
        expect_identical(manual_slopes(), edited)
      }
    )
  })

  it("rejects incompatible profile-key types without disconnecting", {
    fixture <- manual_slopes_fixture()
    bad <- fixture$rule
    bad$USUBJID <- 1
    shiny::testServer(
      manual_slopes_table_server,
      args = list(pknca_data = reactive(fixture$data), manual_slopes_override = reactive(bad)),
      {
        session$flushReact()
        expect_null(manual_slopes())
        expect_false(session$isClosed())
      }
    )
  })

})
