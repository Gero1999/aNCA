#' Build named plot entries eagerly or defer them for app pagination.
#'
#' The shared cache belongs to one immutable data/options generation and retains
#' at most twenty rendered widgets. A new generation discards the old cache.
#' @param ids Ordered group identifiers.
#' @param render Function rendering one group identifier.
#' @param defer Whether to return deferred entries instead of rendered plots.
#' @param cache_size Maximum number of rendered plots retained by a deferred list.
#' @returns A list of plot objects or callable deferred entries.
#' @noRd
.map_tlg_plots <- function(ids, render, defer = FALSE, cache_size = 20L) {
  if (!isTRUE(defer)) return(lapply(ids, render))
  stopifnot(length(cache_size) == 1L, is.finite(cache_size), cache_size >= 1L)
  cache <- new.env(parent = emptyenv())
  cache$plots <- list()
  cache$order <- character()
  make_entry <- function(index) {
    force(index)
    structure(function() {
      key <- as.character(index)
      if (!key %in% cache$order) {
        value <- render(ids[[index]])
        cache$plots[key] <- list(value)
      }
      cache$order <- c(setdiff(cache$order, key), key)
      if (length(cache$order) > cache_size) {
        oldest <- cache$order[[1]]
        cache$plots[[oldest]] <- NULL
        cache$order <- cache$order[-1]
      }
      cache$plots[[key]]
    }, class = "anca_deferred_plot")
  }
  lapply(seq_along(ids), make_entry)
}

#' Realize only selected deferred entries, retaining names and eager inputs.
#' @param items Named list of plots or deferred entries.
#' @returns A list of rendered plot objects.
#' @noRd
.realize_tlg_items <- function(items) {
  if (!is.list(items)) return(items)
  lapply(items, function(item) {
    if (inherits(item, "anca_deferred_plot")) item() else item
  })
}
