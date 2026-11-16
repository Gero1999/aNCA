describe("application dependency boundary", {
  it("keeps optional UI requirements out of core Imports", {
    description <- read.dcf(system.file("DESCRIPTION", package = "aNCA"))
    packages <- function(field) {
      trimws(sub("[[:space:]]*[(].*", "", strsplit(description[1, field], ",")[[1]]))
    }
    expect_false(any(.app_dependencies %in% packages("Imports")))
    expect_true(all(.app_dependencies %in% packages("Suggests")))
    expect_true(all(c("PKNCA", "units", "dplyr") %in% packages("Imports")))
  })

  it("reports the exact missing app packages without trying to install them", {
    available <- function(package, quietly = TRUE) {
      !package %in% c("reactable.extras", "shinyjqui")
    }
    expect_error(check_app_dependencies(available), "'reactable.extras', 'shinyjqui'", fixed = TRUE)
    expect_invisible(check_app_dependencies(function(...) TRUE))
  })
})
