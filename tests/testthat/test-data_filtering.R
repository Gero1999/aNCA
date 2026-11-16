local({
  library(shiny)
  library(logger)
  env <- new.env(parent = environment())
  source(system.file("shiny/modules/tab_data/data_filtering.R", package = "aNCA"), local = env)
  env$reactable_server <- function(...) NULL

  describe("restoring filter panels", {
    it("clears old modules before collecting values and restores fewer filters", {
      state <- new.env()
      state$destroyed <- 0L
      env$.insert_filter_panel <- function(session, filters, filter_counter,
                                           filters_metadata, initial_values = NULL) {
        filter_counter(filter_counter() + 1L)
        value <- reactiveVal(initial_values)
        attr(value, "destroy") <- function() state$destroyed <- state$destroyed + 1L
        filters[[paste0("filter_", filter_counter())]] <- value
      }
      override <- reactiveVal(list(list(column = "SEX", condition = "==", value = "M")))
      data <- reactive(data.frame(SEX = c("M", "F"), AVAL = c(1, 2)))
      args <- list(raw_adnca_data = data, imported_filters = override)
      testServer(env$data_filtering_server, args = args, {
        session$setInputs(submit_filters = 1)
        expect_equal(filtered_data()$SEX, "M")
        expect_length(applied_filters_list(), 1)
        override(list())
        session$flushReact()
        session$setInputs(submit_filters = 2)
        expect_equal(filtered_data()$SEX, c("M", "F"))
        expect_length(applied_filters_list(), 0)
        expect_length(session$userData$applied_filters, 0)
        expect_identical(state$destroyed, 1L)
        expect_false(session$isClosed())
        override(list(list(column = "SEX", condition = "==", value = "F")))
        session$flushReact()
        session$setInputs(submit_filters = 3)
        expect_equal(filtered_data()$SEX, "F")
        expect_length(applied_filters_list(), 1)
        expect_identical(filter_counter(), 2)
        expect_false(session$isClosed())
      })
    })

    it("reports malformed filters without returning unfiltered data or closing the session", {
      override <- reactiveVal(list(list(column = "AVAL", condition = "min-max", value = "bad")))
      data <- reactive(data.frame(AVAL = 1:30))
      args <- list(raw_adnca_data = data, imported_filters = override)
      testServer(env$data_filtering_server, args = args, {
        session$setInputs(submit_filters = 1)
        expect_error(filtered_data(), "Invalid min-max value")
        expect_false(session$isClosed())
        override(list(list(column = "AVAL", condition = "min-max", value = "1-2")))
        session$flushReact()
        session$setInputs(submit_filters = 2)
        expect_equal(filtered_data()$AVAL, 1:2)
        expect_false(session$isClosed())
      })
    })
  })
})

local({
  library(shiny)
  library(shinyWidgets)
  env <- new.env(parent = environment())
  source(system.file("shiny/modules/tab_data/input_filter.R", package = "aNCA"), local = env)

  describe("removed filter module lifecycle", {
    it("stops reading old metadata and inputs after explicit destruction", {
      metadata <- reactiveVal(list(x = list(type = "text", choices = c("a", "b"))))
      testServer(env$input_filter_server, args = list(filters_metadata = metadata), {
        session$setInputs(column = "x", condition = "==", value_select = "a")
        expect_equal(applied_filters()$value, "a")
        attr(applied_filters, "destroy")()
        metadata(list())
        session$setInputs(column = "old", value_select = "b")
        session$flushReact()
        expect_null(applied_filters())
        expect_false(session$isClosed())
      })
    })
  })
})
