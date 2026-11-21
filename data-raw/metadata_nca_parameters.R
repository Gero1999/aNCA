# Load NCA parameters dataset
metadata_nca_parameters <- read.csv(
  "data-raw/metadata_nca_parameters.csv",
  colClasses = c(PKNCA_preferred = "character", min_pknca_version = "character")
)
usethis::use_data(metadata_nca_parameters, overwrite = TRUE)
