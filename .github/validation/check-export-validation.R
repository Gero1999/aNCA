checks <- yaml::read_yaml(".github/validation/export-validation.yml")$checks

expected <- c(
  "standard_output_structure", "cdisc_schema", "settings_structure",
  "script_structure", "slides_structure", "pre_specs_structure",
  "session_info_structure"
)

ids <- vapply(checks, `[[`, "", "id")
blocks_export <- vapply(checks, function(check) isTRUE(check$blocks_export), logical(1))
tests <- vapply(checks, `[[`, "", "test")

stopifnot(
  "Each export check appears exactly once" = identical(sort(ids), sort(expected)),
  "Each check documents a blocking error" = all(blocks_export),
  "Each referenced regression test exists" = all(file.exists(tests))
)
