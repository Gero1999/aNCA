# Settings literals must reconstruct the values, not merely look like R code.
describe("clean_deparse()", {
  roundtrip <- function(x) eval(parse(text = clean_deparse(x)), envir = baseenv())

  it("preserves escaped text, typed missing values, and exact numeric values", {
    values <- list(
      c('Study "A"', "C:\\temp\\data.csv", "line1\nline2", NA_character_),
      c(1.2345678901234567, NA_real_, NaN, Inf, -Inf),
      c(1L, NA_integer_), c(TRUE, FALSE, NA),
      character(), numeric(), integer(), logical(), list(), NULL
    )
    for (x in values) expect_identical(roundtrip(x), x)
  })

  it("preserves names, classes, and data-frame attributes", {
    df <- data.frame(
      `column with spaces` = c('quoted "value"', NA_character_),
      `back\\slash` = factor(c("A", "B")),
      date = as.Date(c("2026-01-01", NA)), check.names = FALSE
    )
    rownames(df) <- c("first row", "second row")
    attr(df[[1]], "label") <- "A label"
    x <- list(df, `a "quoted" name` = c(key = 1.2345678901234567),
              `empty table` = df[0, ], `empty frame` = data.frame())
    expect_identical(roundtrip(x), x)
    expect_identical(roundtrip(setNames(list(1, 2), c("", "NA"))),
                     setNames(list(1, 2), c("", "NA")))
  })

  it("represents tibbles as data frames and preserves long vectors", {
    df <- dplyr::tibble(x = 1:25, y = as.character(1:25))
    expect_identical(roundtrip(df), as.data.frame(df))
    expect_identical(roundtrip(dplyr::tibble()), data.frame())
    expect_identical(roundtrip(rep(c("a", "b", NA_character_), 20)),
                     rep(c("a", "b", NA_character_), 20))
  })
})

describe("get_settings_code: ", {
  data_file <- testthat::test_path("data/test-multispec-ADNCA.csv")
  placeholder <- "settings_list"

  it("reads the requested dataset path instead of a hardcoded relative location", {
    script_file <- withr::local_tempfile(fileext = ".R")
    get_settings_code(test_path("data/test-settings.yaml"), data_file, script_file)
    code <- readLines(script_file)
    env <- new.env(parent = environment())
    eval(parse(text = grep("^(data_path|adnca_data) <-", code, value = TRUE)), env)
    expect_identical(env$data_path, data_file)
    expect_equal(env$adnca_data, read_pk(data_file))
  })

  it("runs the ratio step when an exported settings file contains no ratio rows", {
    settings <- read_settings(test_path("data/test-settings.yaml"))
    settings$settings$ratio_table <- data.frame(TestParameter = character())
    settings$tlg_order <- list(t_pkct01 = list(Selection = FALSE, Comment = "Saved order"))
    settings_file <- withr::local_tempfile(fileext = ".yaml")
    script_file <- withr::local_tempfile(fileext = ".R")
    yaml::write_yaml(settings, settings_file)
    expect_null(read_settings(settings_file)$settings$ratio_table)
    get_settings_code(settings_file, data_file, output_path = script_file)

    # Execute the generated assignment and its downstream calculation, rather
    # than only checking that placeholders disappeared from the script.
    code <- parse(script_file)
    ratio_assignment <- Filter(function(expr) {
      is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("ratio_table"))
    }, as.list(code))
    env <- new.env(parent = environment())
    eval(ratio_assignment[[1]], envir = env)
    expect_identical(nrow(env$ratio_table), 0L)
    result <- calculate_table_ratios(FIXTURE_PKNCA_RES, env$ratio_table)
    expect_equal(result$result$PPORRES, FIXTURE_PKNCA_RES$result$PPORRES)
    expect_equal(nrow(result$result), nrow(FIXTURE_PKNCA_RES$result))
  })

  describe("with legacy settings (no mapping in YAML)", {
    setts_file <- testthat::test_path("data/test-settings.yaml")
    output_file <- tempfile(fileext = ".R")

    get_settings_code(
      settings_file_path = setts_file,
      data_path = data_file,
      output_path = output_file
    )

    it("writes a script R file output and substitutes placeholders", {
      expect_true(file.exists(output_file))
    })

    it("substitutes placeholders in the output file", {
      file_content <- readLines(output_file)
      expect_false(any(grepl(placeholder, file_content)))
    })
  })

  describe("with mapping, filters, and ratio_table in YAML", {
    setts_file <- testthat::test_path("data/test-settings-with-mapping.yaml")
    output_file <- tempfile(fileext = ".R")

    get_settings_code(
      settings_file_path = setts_file,
      data_path = data_file,
      output_path = output_file
    )

    it("writes a script R file output and substitutes placeholders", {
      expect_true(file.exists(output_file))
    })

    it("substitutes placeholders in the output file", {
      file_content <- readLines(output_file)
      expect_false(any(grepl(placeholder, file_content)))
    })

    it("reads mapping from YAML instead of using defaults", {
      file_content <- paste(readLines(output_file), collapse = "\n")
      # Custom mapping uses DOSE_AMT for select_DOSEA
      expect_match(file_content, '"DOSE_AMT"')
    })

    it("reads filters from YAML", {
      file_content <- paste(readLines(output_file), collapse = "\n")
      expect_match(file_content, '"Analyte01"')
      expect_match(file_content, '"Analyte02"')
    })

    it("reads ratio_table from YAML", {
      file_content <- paste(readLines(output_file), collapse = "\n")
      expect_match(file_content, "MRCMAX")
    })

    it("reads units_table from YAML", {
      file_content <- paste(readLines(output_file), collapse = "\n")
      expect_match(file_content, "ng/mL")
    })

    it("passes blq_imputation_rule to PKNCA_update_data_object (#1445)", {
      file_content <- readLines(output_file)
      # The BLQ rule must also be applied at interval setup, matching the app, so
      # the exported script reproduces the app's BLQ handling. The named argument
      # `blq_imputation_rule =` is unique to that call (the calculation step uses
      # `blq_rule =`), so a single occurrence confirms the fix.
      update_line <- grep("blq_imputation_rule =", file_content, value = TRUE)
      expect_length(update_line, 1)
      code <- as.list(parse(output_file))
      assignment <- Filter(function(x) {
        is.call(x) && identical(x[[1]], as.name("<-")) &&
          identical(x[[2]], as.name("pknca_obj"))
      }, code)
      update_call <- assignment[[1]][[3]][[3]]
      expect_identical(update_call[[1]], as.name("PKNCA_update_data_object"))
      expect_identical(eval(update_call$blq_imputation_rule, baseenv()),
                       read_settings(setts_file)$settings$data_imputation$blq_imputation_rule)
    })

    it("derives extra_vars_to_keep from YAML mapping", {
      file_content <- paste(readLines(output_file), collapse = "\n")
      # Custom grouping vars are COHORT, SEX (not TRT01A, RACE from default)
      expect_match(file_content, '"COHORT"')
      expect_match(file_content, '"SEX"')
      expect_match(file_content, '"DOSEA"')
      expect_match(file_content, '"ATPTREF"')
      expect_match(file_content, '"ROUTE"')
    })
  })
})

describe("get_session_code: ", {
  setts_file <- testthat::test_path("data/test-settings.yaml")
  data_file <- testthat::test_path("data/test-multispec-ADNCA.csv")
  output_file <- tempfile(fileext = ".R")
  placeholder <- "settings_list"

  it("reads the input dataset beside the ZIP session script", {
    export_dir <- withr::local_tempdir()
    saveRDS(data.frame(USUBJID = "S1", AVAL = 3.597), file.path(export_dir, "input_data.rds"))
    script <- file.path(export_dir, "session_code.R")
    get_session_code(list(userData = list()), script)
    code <- readLines(script)
    withr::local_dir(export_dir)
    env <- new.env(parent = environment())
    eval(parse(text = grep("^(data_path|adnca_data) <-", code, value = TRUE)), env)
    expect_identical(env$data_path, "input_data.rds")
    expect_equal(env$adnca_data$AVAL, 3.597)
  })

  default_mapping <- list(
    select_STUDYID = "STUDYID",
    select_USUBJID = "USUBJID",
    select_DOSEA = "DOSEA",
    select_DOSEU = "DOSEU",
    select_DOSETRT = "DOSETRT",
    select_PARAM = "PARAM",
    select_Metabolites = "Metab-DrugA",
    select_ARRLT = "ARRLT",
    select_NRRLT = "NRRLT",
    select_AFRLT = "AFRLT",
    select_NCAwXRS = c("NCA1XRS", "NCA2XRS"),
    select_NFRLT = "NFRLT",
    select_PCSPEC = "PCSPEC",
    select_ROUTE = "ROUTE",
    select_TRTRINT = "TRTRINT",
    select_ADOSEDUR = "ADOSEDUR",
    select_Grouping_Variables = c("TRT01A", "RACE", "SEX"),
    select_RRLTU = "RRLTU",
    select_VOLUME = "VOLUME",
    select_VOLUMEU = "VOLUMEU",
    select_AVAL = "AVAL",
    select_AVALU = "AVALU",
    select_ATPTREF = "ATPTREF"
  )

  # Create a mock session object
  setts <- read_settings(setts_file)
  session <- list(
    userData = list(
      settings_list = setts,
      data_path = data_file,
      mapping = default_mapping,
      ratio_table = data.frame()
    )
  )

  # Call the function and generate the output file
  get_session_code(
    session = session$userData,
    output_path = output_file
  )

  it("writes a script R file output", {
    expect_true(file.exists(output_file))
  })

  it("substitutes placeholders in the output file", {
    file_content <- readLines(output_file)
    expect_false(any(grepl(placeholder, file_content)))
  })
})


describe("generated script dependency checks", {
  it("does not execute package installation commands", {
    template <- system.file("www/templates/script_template.R", package = "aNCA")
    code <- parse(template)
    calls_install <- function(x) {
      if (!is.call(x) && !is.expression(x)) return(FALSE)
      if (is.call(x) && identical(x[[1]], as.name("install.packages"))) return(TRUE)
      any(vapply(as.list(x), calls_install, logical(1)))
    }
    expect_false(calls_install(code))
  })
})


describe("generated script literal roundtrips", {
  it("writes parseable scripts for quoted and Windows dataset paths", {
    for (path in c('/tmp/study "A".csv', "C:\\temp\\data.csv")) {
      out <- withr::local_tempfile(fileext = ".R")
      get_settings_code(test_path("data/test-settings.yaml"), path, out)
      code <- parse(out)
      assignment <- Filter(function(x) {
        is.call(x) && identical(x[[1]], as.name("<-")) &&
          identical(x[[2]], as.name("data_path"))
      }, as.list(code))
      env <- new.env(parent = baseenv())
      eval(assignment[[1]], env)
      expect_identical(env$data_path, path)
    }
  })

  it("preserves settings values through actual template substitution", {
    template <- withr::local_tempfile(fileext = ".R")
    out <- withr::local_tempfile(fileext = ".R")
    writeLines("x <- settings_list$test", template)
    values <- list(text = c('A "quote"', NA_character_),
                   number = 1.2345678901234567,
                   table = data.frame(`a b` = 1L, check.names = FALSE))
    get_code(list(test = values), out, template)
    env <- new.env(parent = baseenv())
    sys.source(out, env)
    expect_identical(env$x, values)
  })
})
