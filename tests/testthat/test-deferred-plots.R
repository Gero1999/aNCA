describe("deferred plot entries", {
  it("builds only requested entries and evicts old widgets from its bounded cache", {
    calls <- integer()
    entries <- .map_tlg_plots(1:5, function(id) {
      calls <<- c(calls, id)
      id * 10
    }, defer = TRUE, cache_size = 2)
    names(entries) <- letters[1:5]
    expect_length(calls, 0)
    expect_equal(.realize_tlg_items(entries[1:2]), list(a = 10, b = 20))
    expect_equal(calls, 1:2)
    expect_equal(.realize_tlg_items(entries[2]), list(b = 20))
    expect_equal(calls, 1:2)
    expect_equal(.realize_tlg_items(entries[3]), list(c = 30))
    expect_equal(.realize_tlg_items(entries[1]), list(a = 10))
    expect_equal(calls, c(1L, 2L, 3L, 1L))
    expect_equal(.realize_tlg_items(list(a = 42)), list(a = 42))
    expect_identical(.realize_tlg_items("not available"), "not available")
  })

  it("defers actual individual and combined plot construction without changing output", {
    calls <- 0L
    original <- .pkcg_plotly
    testthat::local_mocked_bindings(.pkcg_plotly = function(...) {
      calls <<- calls + 1L
      original(...)
    }, .package = "aNCA")
    data <- subset(adnca_example, PARAM == "DrugA" & PCSPEC == "SERUM")
    for (render in list(g_pkcg01_lin, g_pkcg02_lin)) {
      before <- calls
      plan <- render(data, .lazy = TRUE)
      expect_gt(length(plan), 1)
      expect_equal(calls, before)
      page <- .realize_tlg_items(plan[1:2])
      expect_equal(calls, before + 2L)
      expect_equal(.realize_tlg_items(plan[1:2]), page)
      expect_equal(calls, before + 2L)
      eager <- render(data)
      expect_equal(names(plan), names(eager))
      expect_equal(page[[1]]$x$data, eager[[1]]$x$data)
      expect_equal(page[[2]]$x$layout, eager[[2]]$x$layout)
      expect_equal(attr(page[[1]], "ggplot")$data,
                   attr(eager[[1]], "ggplot")$data)
      expect_equal(names(.realize_tlg_items(plan)), names(eager))
    }
  })
})
