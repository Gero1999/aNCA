local({
  library(shiny)
  library(logger)
  library(reactable)
  library(reactable.extras)
  env <- new.env(parent = environment())
  shiny_dir <- system.file("shiny", package = "aNCA")
  source(file.path(shiny_dir, "modules/tab_nca/setup/manual_slopes_table.R"), local = env)
  env$reset_reactable_memory <- function(...) NULL

  describe("manual slope edit handlers", {
    it("keeps registration bounded across edits, clearing, and adding group columns", {
      state <- new.env()
      state$registrations <- 0L
      env$observeEvent <- function(...) {
        state$registrations <- state$registrations + 1L
        eval.parent(substitute(shiny::observeEvent(...)))
      }
      testServer(
        env$manual_slopes_table_server,
        args = list(pknca_data = reactive(NULL), manual_slopes_override = reactive(NULL)),
        {
          session$flushReact()
          manual_slopes(data.frame(TYPE = "Exclusion", RANGE = "1", REASON = ""))
          session$flushReact()
          initial_count <- state$registrations
          for (i in 1:20) {
            session$setInputs(edit_RANGE = list(row = 1, column = "RANGE", value = as.character(i)))
            expect_identical(manual_slopes()$RANGE, as.character(i))
          }
          expect_identical(state$registrations, initial_count)
          manual_slopes(NULL)
          session$flushReact()
          session$setInputs(edit_RANGE = list(row = 2, column = "RANGE", value = "stale"))
          expect_null(manual_slopes())
          manual_slopes(data.frame(TYPE = "Exclusion", RANGE = "2", REASON = "", USUBJID = "A"))
          session$flushReact()
          expect_identical(state$registrations, initial_count + 1L)
          expect_identical(manual_slopes()$RANGE, "2")
          session$setInputs(edit_USUBJID = list(row = 1, column = "USUBJID", value = "B"))
          expect_identical(manual_slopes()$USUBJID, "B")
          expect_identical(state$registrations, initial_count + 1L)
          expect_false(session$isClosed())
        }
      )
    })
  })
})
