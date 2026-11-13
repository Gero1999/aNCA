local({
  library(shiny)
  library(shinyWidgets)
  env <- new.env(parent = environment())
  shiny_dir <- system.file("shiny", package = "aNCA")
  source(file.path(shiny_dir, "modules/tab_nca/setup/page_and_searcher.R"), local = env)

  describe("slope plot subject search", {
    it("stores the original subject ID separately from plot names", {
      result <- suppressWarnings(get_halflife_plots(FIXTURE_PKNCA_DATA))
      expect_gt(length(result$plots), 0)
      subjects <- vapply(result$plots, attr, character(1), which = "anca_subject")
      expected <- vapply(result$data, function(x) as.character(x$USUBJID[1]), character(1))
      expect_identical(subjects, expected)
    })

    it("matches complete IDs literally and survives special characters and no matches", {
      ids <- c("1", "10", "SUB[01", "A.B", "A_B", "A+B", "A(B)")
      plots <- lapply(ids, function(id) structure(list(), anca_subject = id))
      names(plots) <- paste0("STUDYID=S_USUBJID=", ids,
                             "_GROUP_WITH_UNDERSCORES=Test_PARAM=A_start=0_end=24")
      selected <- reactiveVal(NULL)
      testServer(
        env$page_and_searcher_server,
        args = list(search_subject = selected, plot_outputs = reactive(plots),
                    plots_per_page = reactive(3)),
        {
          session$flushReact()
          expect_true(all(is_plot_searched()))
          for (id in ids) {
            selected(id)
            session$flushReact()
            expect_identical(ids[is_plot_searched()], id)
            expect_false(session$isClosed())
          }
          selected(c("1", "SUB[01"))
          session$flushReact()
          expect_identical(ids[is_plot_searched()], c("1", "SUB[01"))
          selected("missing")
          session$flushReact()
          expect_false(any(is_plot_searched()))
          expect_equal(num_pages(), 0)
          expect_false(session$isClosed())
          selected(NULL)
          session$flushReact()
          expect_true(all(is_plot_searched()))
          expect_equal(num_pages(), 3)
        }
      )
    })
  })
})
