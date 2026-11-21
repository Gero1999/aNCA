describe("NCA Results", {
  it("displays the calculated values when two session unit names map to VSSO", {
    library(shiny)
    env <- new.env(parent = environment())
    source(system.file("shiny/modules/tab_nca/nca_results.R", package = "aNCA"), local = env)
    env$reactable_server <- function(...) NULL
    env$selector_label <- function(...) NULL

    calculated <- FIXTURE_PKNCA_RES
    calculated$result <- calculated$result %>%
      filter(PPTESTCD == "CMAX", type_interval == "main", USUBJID %in% c(1, 4)) %>%
      arrange(USUBJID) %>%
      mutate(PPTESTCD = "VSSO", PPORRES = c(1500, 2750), PPORRESU = "mL",
             PPSTRES = PPORRES / 1000, PPSTRESU = "L")
    testServer(env$nca_results_server, args = list(
      pknca_data = reactive(calculated$data), res_nca = reactive(calculated),
      settings = reactive(list()), ratio_table = reactive(NULL), grouping_vars = reactive(NULL)
    ), {
      session$userData$units_table <- reactiveVal(data.frame(
        PPTESTCD = c("vss.obs", "vss.inf.obs"), PPORRESU = "mL",
        PPSTRESU = "L", conversion_factor = 0.001
      ))
      expect_no_warning(session$flushReact())
      session$setInputs(params = "VSSO")
      shown <- arrange(output_results(), USUBJID)
      expect_equal(shown$USUBJID, c(1, 4), ignore_attr = TRUE)
      expect_equal(shown[["VSSO[L]"]], c(1.5, 2.75), ignore_attr = TRUE)
    })
  })
})
