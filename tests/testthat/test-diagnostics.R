local({
  diagnostics <- new.env(parent = environment())
  source(system.file("shiny/functions/diagnostics.R", package = "aNCA"),
         local = diagnostics)

  describe("session diagnostic reports", {
    it("does not change the random stream used by scientific calculations", {
      withr::local_seed(93)
      before <- .Random.seed
      store <- diagnostics$diagnostics_store()
      store$record("NCA calculation")
      expect_identical(.Random.seed, before)
    })

    it("retains only a bounded recent history with distinct references", {
      store <- diagnostics$diagnostics_store(limit = 3L)
      records <- lapply(seq_len(8), function(i) store$record("NCA calculation"))
      entries <- store$snapshot()$entries

      expect_length(entries, 3)
      expect_identical(vapply(entries, `[[`, character(1), "id"),
                       vapply(tail(records, 3), `[[`, character(1), "id"))
      expect_length(unique(vapply(records, `[[`, character(1), "id")), 8)
      expect_true(all(vapply(entries, function(x) nzchar(x$time), logical(1))))
    })

    it("does not accept arbitrary identifiers, HTML, paths, or credentials as operation text", {
      store <- diagnostics$diagnostics_store()
      private_text <- paste(
        "STUDY-PRIVATE SUBJECT-002 /srv/private/study.csv",
        "Bearer secret-token <img src=x onerror=alert(1)>"
      )
      store$record(private_text)
      entry <- store$snapshot()$entries[[1]]
      report <- paste(capture.output(str(store$snapshot())), collapse = "\n")

      expect_identical(entry$operation, "Application operation")
      expect_false(grepl("STUDY-PRIVATE|SUBJECT-002|/srv/private|secret-token|<img", report))
      expect_true(nzchar(entry$summary))
    })

    it("isolates concurrent sessions and releases records on close", {
      first <- diagnostics$diagnostics_store()
      second <- diagnostics$diagnostics_store()
      first_record <- first$record("NCA calculation")
      second_record <- second$record("Data upload")

      expect_length(first$snapshot()$entries, 1)
      expect_length(second$snapshot()$entries, 1)
      expect_false(identical(first_record$id, second_record$id))
      expect_identical(first$snapshot()$entries[[1]]$operation, "NCA calculation")
      expect_identical(second$snapshot()$entries[[1]]$operation, "Data upload")

      first$close()
      first$record("Data mapping")
      expect_length(first$snapshot()$entries, 0)
      expect_length(second$snapshot()$entries, 1)
    })
  })

  describe("recoverable operation boundaries", {
    it("stops failed work, clears stale state, cleans up and permits a retry", {
      session <- shiny::MockShinySession$new()
      on.exit(session$close(), add = TRUE)
      session$userData$diagnostics <- diagnostics$diagnostics_store()
      cleaned <- 0L
      invalidated <- 0L
      stale_result <- 42L
      continued <- FALSE
      notifications <- list()
      logged <- character()
      testthat::local_mocked_bindings(
        log_error = function(...) {
          logged <<- c(logged, as.character(glue::glue(..., .envir = parent.frame())))
        },
        .package = "logger"
      )
      testthat::local_mocked_bindings(
        showNotification = function(ui, ...) {
          notifications[[length(notifications) + 1L]] <<- list(ui = ui, ...)
        },
        .package = "shiny"
      )

      result <- diagnostics$diagnostics_run({
        stop("PRIVATE-SUBJECT /srv/secret.csv token=private-token <script>alert(1)</script>")
        continued <- TRUE
      }, "NCA calculation", session,
      cleanup = function() cleaned <<- cleaned + 1L,
      invalidate = function() {
        invalidated <<- invalidated + 1L
        stale_result <<- NULL
      })

      expect_null(result)
      expect_false(continued)
      expect_null(stale_result)
      expect_identical(cleaned, 1L)
      expect_identical(invalidated, 1L)
      expect_length(session$userData$diagnostics$snapshot()$entries, 1)
      expect_length(notifications, 1)
      expect_length(logged, 1)
      reference <- session$userData$diagnostics$snapshot()$entries[[1]]$id
      expect_true(grepl(reference, logged, fixed = TRUE))
      exposed <- paste(capture.output(str(list(
        notifications, session$userData$diagnostics$snapshot()
      ))), collapse = "\n")
      expect_false(grepl("PRIVATE-SUBJECT|/srv/secret|private-token|<script>", exposed))
      expect_true(grepl("PRIVATE-SUBJECT", logged, fixed = TRUE))

      retry <- diagnostics$diagnostics_run({
        stale_result <- 84L
        stale_result
      }, "NCA calculation", session,
      cleanup = function() cleaned <<- cleaned + 1L,
      invalidate = function() invalidated <<- invalidated + 1L)
      expect_identical(retry, 84L)
      expect_identical(stale_result, 84L)
      expect_identical(cleaned, 2L)
      expect_identical(invalidated, 1L)
      expect_length(session$userData$diagnostics$snapshot()$entries, 1)
      expect_length(logged, 1)
    })

    it("preserves req and validate silent errors without diagnostics or invalidation", {
      session <- shiny::MockShinySession$new()
      on.exit(session$close(), add = TRUE)
      session$userData$diagnostics <- diagnostics$diagnostics_store()
      cleaned <- 0L
      invalidated <- 0L
      cleanup <- function() cleaned <<- cleaned + 1L
      invalidate <- function() invalidated <<- invalidated + 1L

      expect_error(diagnostics$diagnostics_run(
        shiny::req(FALSE), "NCA calculation", session, cleanup, invalidate
      ), class = "shiny.silent.error")
      expect_error(diagnostics$diagnostics_run(
        shiny::validate(shiny::need(FALSE, "Select a dataset first")),
        "NCA calculation", session, cleanup, invalidate
      ), class = "shiny.silent.error")
      expect_identical(cleaned, 2L)
      expect_identical(invalidated, 0L)
      expect_length(session$userData$diagnostics$snapshot()$entries, 0)
    })

    it("still recovers and cleans up when the server log sink fails", {
      session <- shiny::MockShinySession$new()
      on.exit(session$close(), add = TRUE)
      session$userData$diagnostics <- diagnostics$diagnostics_store()
      cleaned <- FALSE
      testthat::local_mocked_bindings(
        log_error = function(...) stop("synthetic full log disk"), .package = "logger"
      )
      testthat::local_mocked_bindings(showNotification = function(...) NULL, .package = "shiny")
      result <- diagnostics$diagnostics_run(stop("synthetic operation error"),
                                            "Data mapping", session,
                                            cleanup = function() cleaned <<- TRUE)
      expect_null(result)
      expect_true(cleaned)
      expect_length(session$userData$diagnostics$snapshot()$entries, 1)
    })
  })

  describe("logger adapter and session lifecycle", {
    it("preserves rotating server-file logging while sending only session-safe metadata", {
      namespace <- "anca-diagnostics-test"
      on.exit(rm(list = namespace, envir = get("namespaces", asNamespace("logger"))),
              add = TRUE)
      setup <- new.env(parent = diagnostics)
      list2env(mget(getNamespaceExports("logger"), asNamespace("logger"), inherits = TRUE),
               envir = setup)
      source(system.file("shiny/functions/utils.R", package = "aNCA"), local = setup)
      setup$log_layout <- function(layout = NULL, index = 1) {
        logger::log_layout(layout, namespace = namespace, index = index)
      }
      setup$log_threshold <- function(level = NULL, index = 1) {
        logger::log_threshold(level, namespace = namespace, index = index)
      }
      setup$log_formatter <- function(formatter = NULL, index = 1) {
        logger::log_formatter(formatter, namespace = namespace, index = index)
      }
      setup$log_appender <- function(appender = NULL, index = 1) {
        logger::log_appender(appender, namespace = namespace, index = index)
      }
      setup$appender_console <- function(lines) invisible(NULL)
      withr::local_dir(withr::local_tempdir())
      dir.create("log")
      old_files <- file.path("log", sprintf("aNCA_app_%04d.log", seq_len(7)))
      file.create(old_files)
      setup$setup_logger()

      session <- shiny::MockShinySession$new()
      on.exit(session$close(), add = TRUE)
      session$userData$diagnostics <- diagnostics$diagnostics_store()
      private_text <- "PRIVATE-SUBJECT /srv/dataset.csv <img src=x>"
      shiny::withReactiveDomain(session, logger::log_warn(private_text, namespace = namespace))
      files <- list.files("log", full.names = TRUE)
      expect_length(files, 5)
      expect_false(any(file.exists(head(old_files, 3))))
      server_text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
      expect_true(grepl(private_text, server_text, fixed = TRUE))
      entries <- session$userData$diagnostics$snapshot()$entries
      expect_length(entries, 1)
      expect_identical(entries[[1]]$level, "WARN")
      exposed <- paste(capture.output(str(entries)), collapse = "\n")
      expect_false(grepl("PRIVATE-SUBJECT|/srv/dataset|<img", exposed))
    })

    it("bounds rotated server logs by bytes and count, including multibyte writes", {
      directory <- withr::local_tempdir()
      file <- file.path(directory, "server.log")
      append <- diagnostics$diagnostics_file_appender(file, max_bytes = 128, max_files = 3)
      append("normal log text")
      expect_identical(readLines(file), "normal log text")
      for (i in seq_len(20)) {
        append(paste0(sprintf("record-%02d ", i), paste(rep("\u00e9", 80), collapse = "")))
      }
      files <- list.files(directory, full.names = TRUE)
      expect_length(files, 3)
      expect_true(all(file.info(files)$size <= 128))
      expect_true(grepl("record-20", paste(readLines(file), collapse = "\n"), fixed = TRUE))
    })

    it("routes safe metadata to the current session without sharing messages or references", {
      first <- shiny::MockShinySession$new()
      second <- shiny::MockShinySession$new()
      on.exit(first$close(), add = TRUE)
      on.exit(second$close(), add = TRUE)
      previous_options <- options("shiny.error", "shiny.sanitize.errors")
      diagnostics$diagnostics_session(first)
      diagnostics$diagnostics_session(second)
      expect_identical(options("shiny.error", "shiny.sanitize.errors"), previous_options)

      private_message <- "PRIVATE-ID /srv/patient.csv password=secret <script>alert(1)</script>"
      safe_lines <- diagnostics$diagnostics_log_layout(logger::WARN, private_message)
      expect_identical(safe_lines, "WARN")
      shiny::withReactiveDomain(first, diagnostics$diagnostics_log_appender(safe_lines))
      shiny::withReactiveDomain(second, diagnostics$diagnostics_log_appender("INFO"))

      first_entries <- first$userData$diagnostics$snapshot()$entries
      second_entries <- second$userData$diagnostics$snapshot()$entries
      expect_length(first_entries, 1)
      expect_length(second_entries, 1)
      expect_identical(first_entries[[1]]$level, "WARN")
      expect_identical(second_entries[[1]]$level, "INFO")
      expect_false(identical(first_entries[[1]]$id, second_entries[[1]]$id))
      rendered <- paste(capture.output(str(first_entries)), collapse = "\n")
      expect_false(grepl("PRIVATE-ID|/srv/patient|password|<script>", rendered))

      duplicate <- diagnostics$diagnostics_log_layout(
        logger::ERROR, "Diagnostic reference example: NCA calculation failed"
      )
      expect_length(duplicate, 0)
      shiny::withReactiveDomain(first, diagnostics$diagnostics_log_appender(duplicate))
      expect_length(first$userData$diagnostics$snapshot()$entries, 1)
      diagnostics$diagnostics_log_appender("WARN")
      expect_length(first$userData$diagnostics$snapshot()$entries, 1)
    })

    it("keeps unexpected-error hooks session-local and releases data on session end", {
      new_session <- function() {
        session <- new.env(parent = emptyenv())
        session$userData <- new.env(parent = emptyenv())
        closed <- FALSE
        session$isClosed <- function() closed
        session$sendCustomMessage <- function(...) NULL
        session$onFlushed <- function(...) NULL
        session$onSessionEnded <- function(callback) session$userData$on_end <- callback
        session$onUnhandledError <- function(callback) session$userData$on_error <- callback
        session$close <- function() {
          closed <<- TRUE
          session$userData$on_end()
        }
        session
      }
      # A small session double keeps callback registration instance-local;
      # MockShinySession shares its error callback registry in some Shiny versions.
      first <- new_session()
      second <- new_session()
      on.exit(second$close(), add = TRUE)
      diagnostics$diagnostics_session(first)
      diagnostics$diagnostics_session(second)
      notifications <- 0L
      testthat::local_mocked_bindings(
        showNotification = function(...) notifications <<- notifications + 1L,
        .package = "shiny"
      )
      testthat::local_mocked_bindings(log_error = function(...) NULL, .package = "logger")

      first$userData$on_error(simpleError("private unexpected failure"))
      expect_length(first$userData$diagnostics$snapshot()$entries, 1)
      expect_length(second$userData$diagnostics$snapshot()$entries, 0)
      expect_identical(notifications, 1L)

      silent <- structure(list(message = "not ready", call = NULL),
                          class = c("shiny.silent.error", "error", "condition"))
      first$userData$on_error(silent)
      expect_length(first$userData$diagnostics$snapshot()$entries, 1)
      expect_identical(notifications, 1L)
      first$close()
      expect_length(first$userData$diagnostics$snapshot()$entries, 0)
      second$userData$on_error(simpleError("another private failure"))
      expect_length(second$userData$diagnostics$snapshot()$entries, 1)
    })
  })
})
