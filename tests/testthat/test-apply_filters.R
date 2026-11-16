describe("apply_filters works correctly", {
  test_that("apply_filters works with ==", {
    filter <- list(
      list(
        column = "mpg",
        condition = "==",
        value = "21"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 2)
  })

  test_that("apply_fitlers works with >", {
    filter <- list(
      list(
        column = "mpg",
        condition = ">",
        value = "21"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 12)
  })

  test_that("apply_fitlers works with <", {
    filter <- list(
      list(
        column = "mpg",
        condition = "<",
        value = "21"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 18)
  })

  test_that("apply_fitlers works with >=", {
    filter <- list(
      list(
        column = "mpg",
        condition = ">=",
        value = "21"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 14)
  })

  test_that("apply_fitlers works with <=", {
    filter <- list(
      list(
        column = "mpg",
        condition = "<=",
        value = "21"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 20)
  })

  test_that("apply_filters works with compound filters", {
    filter <- list(
      list(
        column = "mpg",
        condition = "==",
        value = "21"
      ),
      list(
        column = "wt",
        condition = "==",
        value = "2.875"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 1)
  })

  test_that("apply_filters works with named filters", {
    filter <- list(
      filter1 = list(
        column = "mpg",
        condition = ">",
        value = "21"
      ),
      filter2 = list(
        column = "disp",
        condition = "<",
        value = "100"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 5)
  })

  test_that("apply_filters works with !=", {
    filter <- list(
      list(
        column = "mpg",
        condition = "!=",
        value = "21"
      )
    )

    expect_equal(nrow(apply_filters(mtcars, filter)), 30)
  })

  test_that("apply_filters works with min-max", {
    filter <- list(
      filter1  <- list(
        column = "mpg",
        condition = "min-max",
        value = c("15", "20")
      )
    )
    expect_equal(nrow(apply_filters(mtcars, filter)), 13)
  })

  test_that("apply_filters errors with min-max argument NA", {
    filter <- list(
      filter2  <- list(
        column = "mpg",
        condition = "min-max",
        value = c("15", NA)
      )
    )
    expect_error(apply_filters(mtcars, filter), "Invalid min-max value")
  })

  test_that("apply_filters skips NULL filters", {
    filters_with_null <- list(
      list(
        column = "mpg",
        condition = ">",
        value = "21"
      ),
      NULL
    )

    expect_equal(nrow(apply_filters(mtcars, filters_with_null)), 12)
  })

})

describe("apply_filters handles (NA) sentinel", {
  it("== with (NA) includes rows where column is NA", {
    df <- data.frame(x = c("a", "b", NA, NA, "a"), stringsAsFactors = FALSE)
    filter <- list(list(column = "x", condition = "==", value = "(NA)"))
    result <- apply_filters(df, filter)
    expect_equal(nrow(result), 2)
    expect_true(all(is.na(result$x)))
  })

  it("== with (NA) and other values includes both", {
    df <- data.frame(x = c("a", "b", NA, NA, "a"), stringsAsFactors = FALSE)
    filter <- list(list(column = "x", condition = "==", value = c("a", "(NA)")))
    result <- apply_filters(df, filter)
    expect_equal(nrow(result), 4)
  })

  it("!= with (NA) keeps NA rows and excludes matched values", {
    df <- data.frame(x = c("a", "b", NA, "a"), stringsAsFactors = FALSE)
    filter <- list(list(column = "x", condition = "!=", value = c("a", "(NA)")))
    result <- apply_filters(df, filter)
    expect_equal(nrow(result), 2)
    expect_equal(result$x, c("b", NA))
  })

  it("!= with only (NA) returns all rows", {
    df <- data.frame(x = c("a", "b", NA), stringsAsFactors = FALSE)
    filter <- list(list(column = "x", condition = "!=", value = "(NA)"))
    result <- apply_filters(df, filter)
    expect_equal(nrow(result), 3)
  })

  it("== without (NA) excludes NA rows", {
    df <- data.frame(x = c("a", "b", NA), stringsAsFactors = FALSE)
    filter <- list(list(column = "x", condition = "==", value = "a"))
    result <- apply_filters(df, filter)
    expect_equal(nrow(result), 1)
    expect_equal(result$x, "a")
  })
})

describe("apply_filters throws errors for invalid input", {
  test_that("apply_filters throws error when filter misses required fields", {
    filter <- list(
      filter1 = list(
        column = "mpg",
        condition = ">"
      )
    )

    expect_error(
      apply_filters(mtcars, filter),
      "Missing required filter fields: value"
    )
  })

  test_that("apply_filter throws error if column does not exist in data", {
    filter <- list(
      filter1 = list(
        column = "nonexistant",
        condition = ">",
        value = "21"
      )
    )

    expect_error(
      apply_filters(mtcars, filter),
      "Data is missing filtered column: nonexistant"
    )
  })
})

describe("set exclusion and signed numeric ranges", {
  it("excludes a set without recycling and keeps the existing NA policy", {
    data <- data.frame(id = 1:7, x = c(1, 2, 3, 1, 2, 3, NA))
    filter <- list(list(column = "x", condition = "!=", value = c("1", "2")))
    expect_identical(apply_filters(data, filter)$id, c(3L, 6L))
    filter[[1]]$value <- c("1", "2", "(NA)")
    expect_identical(apply_filters(data, filter)$id, c(3L, 6L, 7L))
    filter[[1]]$value <- "1"
    expect_identical(apply_filters(data, filter)$id, c(2L, 3L, 5L, 6L))
  })

  it("parses negative, decimal and exponent bounds", {
    data <- data.frame(x = c(-3, -2, -1, 0, 0.001, 0.002, 1, 2, 3))
    selected <- function(value) {
      apply_filters(data, list(list(column = "x", condition = "min-max", value = value)))$x
    }
    expect_equal(selected("-2-2"), data$x[2:8])
    expect_equal(selected("-3--1"), c(-3, -2, -1))
    expect_equal(selected(" 1e-3 - 2e-3 "), c(0.001, 0.002))
    expect_equal(selected(c("-2", "-1")), c(-2, -1))
    expect_equal(selected("-.5-+.5"), c(0, 0.001, 0.002))
    for (invalid in list("2-1", "Inf-2", "bad", c("1", NA), "1-2-3")) {
      expect_error(selected(invalid), "Invalid min-max value")
    }
  })

  it("rejects invalid operators and comparison values explicitly", {
    data <- data.frame(x = 1:3)
    expect_error(apply_filters(data, list(list(column = "x", condition = "unknown", value = 1))),
                 "Unsupported filter condition")
    expect_error(apply_filters(data, list(list(column = "x", condition = ">", value = c(1, 2)))),
                 "Comparison filters require one value")
    expect_error(apply_filters(data, list(list(column = "x", condition = "==", value = "bad"))),
                 "Invalid numeric filter value")
  })
})
