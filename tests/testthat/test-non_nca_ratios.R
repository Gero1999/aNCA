# Exercise calculation/cache behavior without constructing browser table widgets.
local({
  library(shiny)
  library(purrr)
  source(system.file("shiny/modules/tab_nca/non_nca_ratios.R", package = "aNCA"),
         local = TRUE)
  selector_label <- function(...) NULL
  reactable_server <- function(...) NULL
}, envir = parent.env(environment()))

describe("matrix ratio calculation recovery", {
  it("notifies, invalidates old export results and recovers after unit correction", {
    fixture <- subset(adnca_example, USUBJID == "S1-01" & PARAM == "DrugA" & PCSPEC == "SERUM")
    fixture$PCSPEC <- "BLOOD"
    fixture$AVAL <- 1
    fixture$AVALU <- "ug/mL"
    plasma <- fixture
    plasma$PCSPEC <- "PLASMA"
    plasma$AVAL <- 1000
    plasma$AVALU <- "ng/mL"
    original <- PKNCA_create_data_object(rbind(fixture, plasma))
    data <- shiny::reactiveVal(original)
    notifications <- character()
    testthat::local_mocked_bindings(showNotification = function(ui, ...) {
      notifications <<- c(notifications, as.character(ui))
    }, .package = "shiny")

    shiny::testServer(non_nca_ratio_server, args = list(
      data = data, grouping_vars = shiny::reactive("NFRLT")
    ), {
      session$userData$results <- shiny::reactiveValues(additional_analysis = list())
      session$setInputs(add_pair = 1)
      session$setInputs(spec1_1 = "BLOOD", spec2_1 = "PLASMA", summary_groups = character())
      session$flushReact()
      session$elapse(600)
      session$flushReact()
      expect_true(all(results()$Ratio == 1))
      expect_gt(nrow(session$userData$results$additional_analysis$matrix_ratios), 0)

      for (bad_unit in c("hour", NA_character_, "")) {
        prior_notifications <- length(notifications)
        invalid <- original
        invalid$conc$data$AVALU[invalid$conc$data$PCSPEC == "PLASMA"] <- bad_unit
        data(invalid)
        session$flushReact()
        session$elapse(600)
        session$flushReact()
        expect_null(results())
        expect_null(session$userData$results$additional_analysis$matrix_ratios)
        expect_gt(length(notifications), prior_notifications)
        expect_match(tail(notifications, 1), "concentration units must be")
        expect_false(session$isClosed())

        data(original)
        session$flushReact()
        session$elapse(600)
        session$flushReact()
        expect_true(all(results()$Ratio == 1))
        expect_gt(nrow(session$userData$results$additional_analysis$matrix_ratios), 0)
      }
    })
  })
})
