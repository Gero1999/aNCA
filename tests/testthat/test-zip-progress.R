local({
  library(shiny)
  env <- new.env(parent = environment())
  shiny_dir <- system.file("shiny", package = "aNCA")
  source(file.path(shiny_dir, "functions/zip-utils.R"), local = env)
  source(file.path(shiny_dir, "modules/tab_nca/zip.R"), local = env)
  env$prepare_export_files <- function(target_dir, ...) {
    writeLines("example", file.path(target_dir, "example.txt"))
  }

  describe("ZIP export progress cleanup", {
    it("closes progress after successful export and failure, allowing a retry", {
      state <- new.env()
      state$messages <- character()
      testServer(
        env$zip_server,
        args = list(res_nca = reactive(NULL), adnca_data = reactive(NULL),
                    settings = reactive(list()), grouping_vars = reactive(character())),
        {
          root_session <- session$rootScope()
          root_session$sendProgress <- function(type, message) {
            state$messages <- c(state$messages, type)
          }
          target <- tempfile(fileext = ".zip")
          on.exit(unlink(target), add = TRUE)
          original_wd <- getwd()
          .run_export(target)
          expect_true(file.exists(target))
          expect_equal(unzip(target, list = TRUE)$Name, "example.txt")
          expect_identical(tail(state$messages, 1), "close")
          expect_equal(sum(state$messages == "close"), 1)
          expect_identical(getwd(), original_wd)

          writer <- env$prepare_export_files
          env$prepare_export_files <- function(...) stop("simulated export failure")
          state$messages <- character()
          expect_error(.run_export(target), "simulated export failure")
          expect_false(file.exists(target))
          expect_identical(tail(state$messages, 1), "close")
          expect_equal(sum(state$messages == "close"), 1)
          expect_identical(getwd(), original_wd)

          env$prepare_export_files <- writer
          state$messages <- character()
          .run_export(target)
          expect_true(file.exists(target))
          expect_identical(tail(state$messages, 1), "close")
          expect_false(session$isClosed())
        }
      )
    })
  })
})
