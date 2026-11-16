# Sample concentration data for testing with three doses
conc_data <- data.frame(
  USUBJID = c("001", "001", "001", "001", "001", "001", "001", "001", "001", "001"),
  AVAL = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  DOSNOA = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
  PROFILE = c("A", "A", "A", "B", "B", "B", "B", "C", "C", "C"),
  ARRLT = c(-1, 0, 1, -1, 0, 1, 2, 0, 1, 2),
  AFRLT = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8),
  NRRLT = c(-1, 0, 1, -1, 0, 1, 2, 0, 1, 2),
  NFRLT = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8)
)

result <- dose_profile_duplicates(conc_data, groups = c("USUBJID", "DOSNOA"), dosno = "DOSNOA")

describe("dose_profile_duplicates", {

  it("should maintain the original columns", {
    expect_true(all(colnames(conc_data) %in% colnames(result)))
  })

  it("should add predose duplicates correctly", {
    predose_rows <- result[result$DOSNOA == 1, ]
    expect_equal(nrow(predose_rows), 4)
  })

  it("should add last dose values when predose is missing", {
    last_dose_rows <- result[result$DOSNOA == 3, ]
    expect_equal(nrow(last_dose_rows), 4)
  })

  it("should return data with IX column if only one dose is present", {
    single_dose_data <- conc_data %>%
      filter(DOSNOA == 1)

    result_single_dose <- dose_profile_duplicates(single_dose_data,
                                                  groups = c("USUBJID", "DOSNOA"),
                                                  dosno = "DOSNOA")
    expect_true("IX" %in% colnames(result_single_dose))
  })

  it("handles character or numeric dosno correctly", {
    result_char_dosno <- dose_profile_duplicates(conc_data,
                                                 groups = c("USUBJID", "PROFILE"),
                                                 dosno = "PROFILE")
    expect_equal(result_char_dosno$PROFILE[1:4], c("A", "A", "A", "A"))
  })
})

describe("physical-time predose selection", {
  profile <- data.frame(
    USUBJID = "S1", PARAM = "DrugA", DOSNOA = c(1, 1, 2),
    AFRLT = c(7, 8, 10), ARRLT = c(7, 8, 0), NFRLT = c(7, 8, 10),
    NRRLT = c(7, 8, 0), AVAL = c(70, 80, 100), RRLTU = "hour"
  )
  predose <- function(x, ...) {
    result <- dose_profile_duplicates(x, ...)
    result[result$DOSNOA == 2 & result$ARRLT < 0, , drop = FALSE]
  }

  it("is invariant to row order and equivalent hour/minute/day inputs", {
    for (unit in c("hour", "min", "day")) {
      multiplier <- c(hour = 1, min = 60, day = 1 / 24)[[unit]]
      x <- profile
      x[c("AFRLT", "ARRLT", "NFRLT", "NRRLT")] <-
        x[c("AFRLT", "ARRLT", "NFRLT", "NRRLT")] * multiplier
      x$RRLTU <- unit
      for (order in list(1:3, c(2, 1, 3), c(3, 2, 1))) {
        actual <- predose(x[order, ])
        expect_equal(actual$AVAL, 80)
        expect_equal(actual$ARRLT / multiplier, -2)
      }
    }
  })

  it("includes the four-hour boundary and excludes earlier and postdose samples", {
    x <- profile[c(1, 3), ]
    x$AFRLT[1] <- x$ARRLT[1] <- x$NFRLT[1] <- x$NRRLT[1] <- 6
    expect_equal(nrow(predose(x)), 1)
    x$AFRLT[1] <- x$ARRLT[1] <- x$NFRLT[1] <- x$NRRLT[1] <- 5.99
    expect_equal(nrow(predose(x)), 0)
    x$AFRLT[1] <- x$ARRLT[1] <- x$NFRLT[1] <- x$NRRLT[1] <- 11
    expect_equal(nrow(predose(x)), 0)
  })

  it("supports legacy inputs and validates explicit or inconsistent units", {
    legacy <- profile
    legacy$RRLTU <- NULL
    expect_equal(predose(legacy)$AVAL, 80)
    expect_equal(predose(legacy, time_unit = "hour")$AVAL, 80)
    expect_error(predose(profile, time_unit = "ng"), "recognized time units")
    expect_error(predose(profile, time_unit = c("hour", "min")), "one unit per row")
    profile$RRLTU[2] <- "min"
    expect_error(predose(profile), "consistent")
  })
})
