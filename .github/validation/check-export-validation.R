rules <- yaml::read_yaml(".github/validation/export-validation.yml")$requirements

expected <- c(
  "standard_output_structure", "cdisc_schema", "settings_structure",
  "script_structure", "slides_structure", "pre_specs_structure",
  "session_info_structure"
)

ids <- vapply(rules, `[[`, "", "id")
validators <- vapply(rules, `[[`, "", "validator")
requirements <- vapply(rules, `[[`, "", "requirement")
blocks_export <- vapply(rules, function(rule) isTRUE(rule$blocks_export), logical(1))
tests <- vapply(rules, `[[`, "", "test")

stopifnot(
  "Each export check appears exactly once" = identical(sort(ids), sort(expected)),
  "Each requirement names its validator" = identical(ids, validators),
  "Each check cites the applicable requirement" = all(requirements == "21 CFR \u00a711.10(a)"),
  "Each check documents a blocking error" = all(blocks_export),
  "Each referenced regression test exists" = all(file.exists(tests))
)
