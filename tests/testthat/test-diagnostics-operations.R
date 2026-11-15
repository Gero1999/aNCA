local({
  library(shiny)
  library(logger)
  app <- new.env(parent = environment())
  shiny_dir <- system.file("shiny", package = "aNCA")
  source(file.path(shiny_dir, "functions", "diagnostics.R"), local = app)
  source(file.path(shiny_dir, "modules", "tab_nca.R"), local = app)

  describe("NCA operation recovery", {
    it("keeps HTML-like package warnings literal", {
      message <- '<img src=x onerror="alert(1)">'
      warning_text <- app$.parse_pknca_warning(simpleWarning(message))
      expect_identical(warning_text, message)
      expect_false(inherits(warning_text, "html"))
      rendered <- as.character(tags$p(warning_text))
      expect_match(rendered, "&lt;img", fixed = TRUE)
      expect_false(grepl("<img", rendered, fixed = TRUE))
    })

    it("clears stale results, cleans up errors and req(), and permits a successful retry", {
      model <- new.env(parent = emptyenv())
      model$mode <- "success"
      model$data <- list(units = data.frame(PPTESTCD = "CMAX", PPSTRESU = "ng/mL",
                                            conversion_factor = 1))
      model$result <- list(result = data.frame(PPTESTCD = "CMAX", PPORRES = 12),
                           data = model$data)
      app$nca_setup_server <- function(...) {
        list(processed_pknca_data = reactive(model$data),
             settings = reactive(list(data_imputation = list(), flags = list())),
             general_settings = list(), ratio_table = reactive(NULL), slope_rules = reactive(NULL))
      }
      app$.validate_nca_run <- function(...) {
        if (model$mode == "silent") req(FALSE)
        model$data
      }
      app$PKNCA_calculate_nca <- function(data, ...) {
        model$received <- data
        if (model$mode == "error") stop("private path /srv/study/patient.csv")
        model$result
      }
      app$mutate <- function(.data, ...) .data
      app$PKNCA_hl_rules_exclusion <- function(data, ...) data
      app$calculate_table_ratios <- function(data, ...) data
      app$remove_pp_not_requested <- identity
      app$loading_popup <- function(...) showModal(modalDialog("Loading"))
      for (name in c("reactable_server", "nca_results_server", "descriptive_statistics_server",
                     "additional_analysis_server", "parameter_plots_server")) {
        app[[name]] <- function(...) NULL
      }
      app$parameter_exclusions_server <- function(...) reactive(NULL)
      app$parameter_datasets_server <- function(...) reactive(list(adpp = NULL))
      local_mocked_bindings(check_valid_pknca_data = identity, .package = "aNCA")

      testServer(app$tab_nca_server, args = list(
        pknca_data = reactiveVal(model$data), extra_group_vars = reactive(character()),
        settings_override = reactive(NULL), auto_replay_ready = reactive(FALSE)
      ), {
        messages <- list()
        modals <- character()
        root_session <- session$rootScope()
        root_session$sendCustomMessage <- function(type, message) {
          messages[[length(messages) + 1L]] <<- list(type = type, message = message)
        }
        root_session$sendModal <- function(type, message) modals <<- c(modals, type)
        session$userData$diagnostics <- app$diagnostics_store()
        session$flushReact()
        session$setInputs(run_nca = 1)
        expect_identical(res_nca(), model$result)
        expect_identical(model$received, model$data)
        expect_true(session$userData$nca_ran)

        session$userData$results$CDISC <- list(adpp = "old")
        model$mode <- "error"
        session$setInputs(run_nca = 2)
        expect_null(res_nca())
        expect_false(session$userData$nca_ran)
        expect_null(session$userData$final_units)
        expect_null(session$userData$results$CDISC)
        expect_length(session$userData$diagnostics$snapshot()$entries, 1L)
        expect_identical(tail(messages, 1)[[1]]$type, "shinyjs-enable")
        Sys.sleep(0.6)
        session$flushReact()
        expect_identical(tail(modals, 1), "remove")

        model$mode <- "silent"
        session$setInputs(run_nca = 3)
        expect_null(res_nca())
        expect_length(session$userData$diagnostics$snapshot()$entries, 1L)
        expect_identical(tail(messages, 1)[[1]]$type, "shinyjs-enable")
        Sys.sleep(0.6)
        session$flushReact()
        expect_identical(tail(modals, 1), "remove")

        model$mode <- "success"
        session$setInputs(run_nca = 4)
        expect_identical(res_nca(), model$result)
        expect_true(session$userData$nca_ran)
        expect_identical(tail(messages, 1)[[1]]$type, "shinyjs-enable")
        pknca_data(NULL)
        session$flushReact()
        expect_null(res_nca())
        expect_false(session$userData$nca_ran)
      })
    })
  })
})
