# Browser reports deliberately contain only allowlisted labels, never condition messages,
# logger message bodies, dataset contents, session tokens, filenames or server paths.
diagnostics_operation <- function(operation) {
  allowed <- c("NCA calculation", "Data preparation", "Data mapping", "Data upload",
               "Data filtering", "Settings upload", "Application operation",
               "Session started", "Application event")
  if (length(operation) == 1L && !is.na(operation) && operation %in% allowed) {
    operation
  } else {
    "Application operation"
  }
}

diagnostics_store <- function(limit = 100L) {
  stopifnot(length(limit) == 1L, is.finite(limit), limit >= 1, limit <= 100)
  entries <- list()
  closed <- FALSE
  # A local reference independent of authentication tokens and R's scientific RNG.
  # Only the generated basename is used; no server directory is exposed.
  prefix <- basename(tempfile("diagnostic-"))
  counter <- 0L
  record <- function(operation, level = "ERROR") {
    if (closed) return(invisible(NULL))
    operation <- diagnostics_operation(operation)
    if (!level %in% c("INFO", "WARN", "ERROR", "SUCCESS")) level <- "INFO"
    counter <<- counter + 1L
    entry <- list(
      id = paste0(prefix, "-", counter),
      time = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      level = level, operation = operation,
      summary = if (level == "ERROR") {
        paste0(operation, " failed. Check your inputs and settings, then retry. ",
               "If it continues, share this diagnostic report with support.")
      } else {
        paste0(operation, ": ", tolower(level), " event recorded.")
      }
    )
    entries <<- tail(c(entries, list(entry)), limit)
    entry
  }
  list(
    record = record,
    snapshot = function() {
      list(package_version = as.character(utils::packageVersion("aNCA")), entries = entries)
    },
    close = function() {
      closed <<- TRUE
      entries <<- list()
    }
  )
}

diagnostics_ui <- function() {
  tagList(
    tags$button("Diagnostics", id = "diagnostics-open", type = "button",
                class = "btn btn-outline-secondary"),
    tags$dialog(
      id = "diagnostics-panel", `aria-labelledby` = "diagnostics-title",
      tags$h2("Session diagnostics", id = "diagnostics-title"),
      tags$p("No raw datasets, error messages, identifiers or server paths are included. ",
             "Detailed logs remain on the server for authorized support."),
      tags$p(id = "diagnostics-status", role = "status"),
      tags$h3("Latest error"), tags$pre(id = "diagnostics-summary"),
      tags$h3("Recent session events"), tags$pre(id = "diagnostics-events"),
      tags$button("Download diagnostic report", id = "diagnostics-download", type = "button",
                  class = "btn btn-primary"),
      tags$button("Close", id = "diagnostics-close", type = "button", class = "btn btn-secondary")
    )
  )
}

# Installed once with the normal logger configuration, never once per session.
# Layout discards message bodies before the appender sees them. No session lookup
# or buffer writes take reactive dependencies, so logging cannot retrigger an output.
diagnostics_log_layout <- function(level, msg, ...) {
  if (any(startsWith(as.character(msg), "Diagnostic reference "))) return(character())
  label <- attr(level, "level")
  if (is.null(label)) as.character(level) else label
}

diagnostics_log_appender <- function(lines) {
  session <- shiny::getDefaultReactiveDomain()
  if (is.null(session) || is.null(session$userData$diagnostics)) return(invisible(NULL))
  for (level in lines) {
    tryCatch(session$userData$diagnostics$record("Application event", level),
             error = function(error) NULL)
  }
  invisible(NULL)
}

# logger's built-in rotation requires an optional package. Keep its file sink,
# but rotate with base R so diagnostics add no runtime dependency.
diagnostics_file_appender <- function(file, max_bytes = 5 * 1024^2, max_files = 5L) {
  append <- logger::appender_file(file)
  function(lines) {
    # Bound an unusually large individual log write as well as the file count.
    lines <- substr(paste(lines, collapse = "\n"), 1L, floor(max_bytes / 4) - 1L)
    size <- if (file.exists(file)) file.info(file)$size else 0
    if (is.finite(size) && size + nchar(lines, type = "bytes") + 1 > max_bytes) {
      unlink(paste0(file, ".", max_files - 1L))
      if (max_files > 2L) {
        for (i in seq.int(max_files - 2L, 1L)) {
          from <- paste0(file, ".", i)
          if (file.exists(from)) file.rename(from, paste0(file, ".", i + 1L))
        }
      }
      if (file.exists(file)) file.rename(file, paste0(file, ".1"))
    }
    append(lines)
  }
}

diagnostics_session <- function(session) {
  store <- diagnostics_store()
  publish <- function() {
    if (!session$isClosed()) {
      tryCatch(session$sendCustomMessage("anca-diagnostics", store$snapshot()),
               error = function(error) NULL)
    }
  }
  record <- function(operation, level = "ERROR") {
    entry <- store$record(operation, level)
    publish()
    entry
  }
  session$userData$diagnostics <- list(
    record = record, snapshot = store$snapshot, close = store$close
  )
  session$onSessionEnded(store$close)
  session$onFlushed(function() record("Session started", "INFO"), once = TRUE)
  # This hook records remaining unexpected failures; recovery belongs at the
  # operation boundary. It does not replace Shiny's own error handling.
  shiny::onUnhandledError(function(error) {
    if (!inherits(error, "shiny.silent.error")) {
      diagnostics_failure("Application operation", session, error)
    }
  }, session = session)
  invisible(session$userData$diagnostics)
}

diagnostics_failure <- function(operation, session, error = NULL) {
  entry <- if (!is.null(session) && !is.null(session$userData$diagnostics)) {
    session$userData$diagnostics$record(operation)
  } else {
    diagnostics_store()$record(operation)
  }
  if (is.null(entry)) return(invisible(NULL))
  # Preserve the cause in the existing admin-only server sinks, never the report.
  # A full disk or broken sink must not turn recovery into a second error.
  detail <- if (is.null(error)) "" else substr(conditionMessage(error), 1L, 8192L)
  tryCatch(logger::log_error("Diagnostic reference {entry$id}: {entry$operation} failed: {detail}"),
           error = function(error) NULL)
  if (!is.null(session) && !session$isClosed()) {
    shiny::showNotification(
      paste(entry$summary, "Reference:", entry$id), type = "error", duration = NULL,
      action = tags$button("Open diagnostics", type = "button", class = "btn btn-link",
                           `data-anca-diagnostics` = "open"), session = session
    )
  }
  invisible(entry)
}

diagnostics_run <- function(expr, operation, session = shiny::getDefaultReactiveDomain(),
                            cleanup = function() {}, invalidate = function() {}) {
  on.exit(cleanup(), add = TRUE)
  tryCatch(force(expr), error = function(error) {
    # req()/validate() are control flow, including silent subclasses.
    if (inherits(error, "shiny.silent.error")) stop(error)
    invalidate()
    diagnostics_failure(operation, session, error)
    NULL
  })
}
