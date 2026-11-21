describe("NCA calculation warnings", {
  it("explains why PKNCA skipped an interval with no data", {
    env <- new.env(parent = environment())
    env$HTML <- htmltools::HTML
    source(system.file("shiny/modules/tab_nca.R", package = "aNCA"), local = env)
    warning <- simpleWarning("Error with interval start=0, end=23.881: No data for interval")

    expect_equal(
      as.character(env$.parse_pknca_warning(warning)),
      paste(
        "NCA skipped interval start=0, end=23.881 because it has no eligible",
        "concentration samples. Check the sample times and any NCA exclusion reasons."
      )
    )
  })
})
