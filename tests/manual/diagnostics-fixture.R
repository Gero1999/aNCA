# Synthetic integration app; run from the checkout with:
# Rscript tests/manual/diagnostics-fixture.R
# No clinical data or production fault-injection controls are used.
library(shiny)
devtools::load_all(".", quiet = TRUE)

shiny_dir <- normalizePath("inst/shiny", mustWork = TRUE)
source(file.path(shiny_dir, "functions", "diagnostics.R"), local = TRUE)
addResourcePath("diagnostics-test", file.path(shiny_dir, "www"))
logger::log_appender(diagnostics_file_appender(
  Sys.getenv("DIAGNOSTICS_TEST_LOG", tempfile("diagnostics-fixture-")),
  max_bytes = 1024 * 1024, max_files = 2L
))

ui <- fluidPage(
  includeCSS(file.path(shiny_dir, "www", "main.css")),
  tags$script(src = "diagnostics-test/diagnostics.js"),
  diagnostics_ui(),
  h1("Synthetic diagnostics fixture"),
  checkboxInput("fixture_invalid", "Fail the next operation", FALSE),
  shinyjs::useShinyjs(),
  actionButton("fixture_run", "Run synthetic operation"),
  textOutput("fixture_result"),
  textOutput("fixture_loading"),
  actionButton("fixture_literal", "Send literal rendering fixture"),
  actionButton("fixture_flood", "Send bounded buffer fixture"),
  actionButton("fixture_disconnect", "Disconnect this synthetic session")
)

server <- function(input, output, session) {
  diagnostics_session(session)
  result <- reactiveVal(NULL)
  loading <- reactiveVal(FALSE)
  output$fixture_result <- renderText(if (is.null(result())) "No current result" else result())
  output$fixture_loading <- renderText(if (loading()) "Loading" else "Idle")

  observeEvent(input$fixture_run, {
    loading(TRUE)
    shinyjs::disable("fixture_run")
    value <- diagnostics_run({
      if (isTRUE(input$fixture_invalid)) {
        stop("Raw synthetic secret=do-not-export; /srv/private/study.csv; subject ALPHA-42")
      }
      "Synthetic result: 42"
    }, "NCA calculation", session = session, cleanup = function() {
      loading(FALSE)
      shinyjs::enable("fixture_run")
    }, invalidate = function() result(NULL))
    result(value)
  })

  observeEvent(input$fixture_literal, {
    session$sendCustomMessage("anca-diagnostics", list(entries = list(list(
      operation = "<b>literal operation</b>",
      summary = '<img src=x onerror="window.diagnosticsInjected=true">',
      id = "literal-fixture", time = "2026-01-01T00:00:00Z", level = "ERROR"
    ))))
  })

  observeEvent(input$fixture_flood, {
    entries <- lapply(seq_len(250), function(i) {
      list(
        operation = paste(rep("Synthetic operation", 100), collapse = " "),
        summary = paste(rep("Synthetic bounded message", 100), collapse = " "),
        id = paste0("flood-", i), time = "2026-01-01T00:00:00Z", level = "INFO"
      )
    })
    session$sendCustomMessage("anca-diagnostics", list(entries = entries))
  })

  observeEvent(input$fixture_disconnect, session$close())
}

runApp(shinyApp(ui, server), host = "127.0.0.1",
       port = as.integer(Sys.getenv("DIAGNOSTICS_TEST_PORT", "3849")), launch.browser = FALSE)
