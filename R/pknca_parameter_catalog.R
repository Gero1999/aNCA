# PKNCA is the stable settings key. runtime_name is the name that can be
# requested from this installation, including parameters registered by aNCA.
.pknca_parameter_catalog <- function(metadata = aNCA::metadata_nca_parameters,
                                     registered = names(PKNCA::get.interval.cols()),
                                     version = utils::packageVersion("PKNCA")) {
  runtime <- metadata$PKNCA
  preferred <- metadata$PKNCA_preferred %in% registered
  runtime[preferred] <- metadata$PKNCA_preferred[preferred]
  runtime[!runtime %in% registered] <- NA_character_

  minimum <- metadata$min_pknca_version
  gated <- which(nzchar(minimum, keepNA = FALSE))
  runtime[gated[numeric_version(minimum[gated]) > version]] <- NA_character_
  metadata$runtime_name <- runtime
  metadata
}

.pknca_settings_keys <- function(params, metadata = aNCA::metadata_nca_parameters) {
  aliases <- metadata[nzchar(metadata$PKNCA_preferred, keepNA = FALSE), , drop = FALSE]
  index <- match(params, aliases$PKNCA_preferred)
  replace(params, !is.na(index), aliases$PKNCA[index[!is.na(index)]])
}

.resolve_pknca_parameters <- function(params, catalog = .pknca_parameter_catalog(),
                                      registered = names(PKNCA::get.interval.cols())) {
  index <- match(.pknca_settings_keys(params, catalog), catalog$PKNCA)
  # Preserve parameters registered by callers outside aNCA's catalog.
  resolved <- replace(params, !is.na(index), catalog$runtime_name[index[!is.na(index)]])
  unavailable <- unique(params[!resolved %in% registered])
  if (length(unavailable)) {
    stop("Parameters unavailable with PKNCA ", utils::packageVersion("PKNCA"), ": ",
         paste(unavailable, collapse = ", "), call. = FALSE)
  }
  resolved
}

# Accept interval objects prepared on another installation, including objects
# with both flags present. Combine requests per row.
.normalize_pknca_intervals <- function(intervals, catalog) {
  keys <- intersect(names(intervals), c(catalog$PKNCA, catalog$PKNCA_preferred))
  keys <- keys[vapply(intervals[keys], is.logical, logical(1))]
  runtime <- catalog$runtime_name[match(.pknca_settings_keys(keys, catalog), catalog$PKNCA)]
  unavailable <- keys[is.na(runtime)]
  requested <- unavailable[vapply(intervals[unavailable], any, logical(1), na.rm = TRUE)]
  .resolve_pknca_parameters(requested, catalog)
  registered <- setdiff(names(PKNCA::get.interval.cols()), c("start", "end"))
  intervals[setdiff(unavailable, registered)] <- NULL

  for (i in which(!is.na(runtime) & keys != runtime)) {
    old <- intervals[[keys[i]]]
    target <- runtime[i]
    intervals[[target]] <- if (target %in% names(intervals)) intervals[[target]] | old else old
    if (keys[i] %in% registered) {
      intervals[[keys[i]]] <- rep(FALSE, nrow(intervals))
    } else {
      intervals[[keys[i]]] <- NULL
    }
  }
  # pk.nca expects the full set of flags on an already constructed PKNCAdata.
  for (key in setdiff(registered, names(intervals))) {
    intervals[[key]] <- rep(FALSE, nrow(intervals))
  }
  intervals
}
