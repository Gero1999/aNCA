describe("PKNCA parameter compatibility", {
  it("resolves both spellings and checks minimum versions numerically", {
    meta <- subset(metadata_nca_parameters, PKNCA %in% c("vss.obs", "vss.int.last"))
    meta$min_pknca_version[meta$PKNCA == "vss.int.last"] <- "0.12.10"
    old <- c("vss.obs", "vss.int.last", "custom")
    catalog <- .pknca_parameter_catalog(meta, old, numeric_version("0.12.9"))
    expect_equal(.resolve_pknca_parameters(c("vss.inf.obs", "custom"), catalog, old),
                 c("vss.obs", "custom"))
    expect_error(.resolve_pknca_parameters("vss.int.last", catalog, old), "vss.int.last")

    for (registered in list(c(old, "vss.inf.obs"), c("vss.inf.obs", "vss.int.last"))) {
      catalog <- .pknca_parameter_catalog(meta, registered, numeric_version("0.12.10"))
      expect_equal(.resolve_pknca_parameters(c("vss.obs", "vss.int.last"), catalog, registered),
                   c("vss.inf.obs", "vss.int.last"))
    }
  })

  it("merges requests in prebuilt intervals row by row", {
    inputs <- data.frame(start = 0:2, end = 1:3, vss.obs = c(TRUE, FALSE, FALSE),
                         vss.inf.obs = c(FALSE, TRUE, FALSE))
    got <- .normalize_pknca_intervals(inputs, .pknca_parameter_catalog())
    preferred <- .resolve_pknca_parameters("vss.obs")
    expect_identical(got[[preferred]], c(TRUE, TRUE, FALSE))
    alias <- setdiff(c("vss.obs", "vss.inf.obs"), preferred)
    expect_false(any(got[[alias]]))
  })

  it("preserves both input spellings in CDISC translations", {
    catalog <- metadata_nca_parameters
    aliases <- subset(catalog, nzchar(PKNCA_preferred))
    expect_equal(translate_terms(aliases$PKNCA_preferred), translate_terms(aliases$PKNCA))
    expect_identical(.pknca_settings_keys(aliases$PKNCA_preferred), aliases$PKNCA)
  })

  it("calculates VSSIL only where PKNCA provides the distinct interval estimator", {
    example <- subset(adnca_example, USUBJID == "S1-01" & PCSPEC == "SERUM" & PARAM == "DrugA")
    data <- PKNCA_create_data_object(example)
    selections <- setNames(list("vss.int.last"), unique(.derive_study_types(data)$type))
    prepare <- function() {
      PKNCA_update_data_object(
        data, method = "lin up/log down", selected_analytes = "DrugA",
        selected_profile = unique(example$ATPTREF), selected_pcspec = "SERUM",
        parameter_selections = selections
      )
    }
    if ("vss.int.last" %in% names(PKNCA::get.interval.cols())) {
      raw <- PKNCA_calculate_nca(prepare())
      raw$result$PPTESTCD <- translate_terms(raw$result$PPTESTCD)
      result <- remove_pp_not_requested(raw)$result
      expect_equal(unique(result$PPTESTCD), "VSSIL")
      expect_true(any(is.finite(result$PPORRES)))
    } else {
      expect_error(prepare(), "Parameters unavailable with PKNCA .*vss.int.last")
    }
  })
})

describe("restoring parameter selections", {
  it("restores either spelling and reports unavailable parameters for review", {
    library(shiny)
    env <- new.env(parent = environment())
    source(system.file("shiny/modules/tab_nca/setup/parameter_selection.R", package = "aNCA"),
           local = env)
    env$reactable_server <- function(...) NULL
    env$runjs <- function(...) NULL
    example <- subset(adnca_example, USUBJID == "S1-01" & PCSPEC == "SERUM" & PARAM == "DrugA")
    data <- PKNCA_create_data_object(example)
    type <- unique(.derive_study_types(data)$type)
    overrides <- reactiveVal(setNames(list(c("vss.inf.obs", "vss.int.last", "absent.param")), type))

    testServer(env$parameter_selection_server, args = list(
      processed_pknca_data = reactive(data), parameter_override = overrides
    ), {
      session$flushReact()
      selected <- parameter_lists_by_type()[[type]]
      expect_true("vss.obs" %in% selected)
      expect_false("vss.inf.obs" %in% selected)
      available <- "vss.int.last" %in% names(PKNCA::get.interval.cols())
      expect_identical("vss.int.last" %in% selected, available)
      expect_match(skipped(), "absent.param", fixed = TRUE)
      expect_identical(grepl("vss.int.last", skipped(), fixed = TRUE), !available)

      overrides(setNames(list("vss.obs"), type))
      session$flushReact()
      expect_equal(parameter_lists_by_type()[[type]], "vss.obs")
      expect_length(skipped(), 0)
    })
  })
})
