# Preserve CSV identifiers and restore filters safely

## How to test

1. Start this branch's app (`devtools::load_all(export_all = FALSE); aNCA::run_app()`).
2. Upload [csv-identifiers.csv](tests/manual/data/csv-identifiers.csv). In **Mapping**, map **USUBJID → participant**; retain the matching-name mappings for other fields. Continue to **Filtering**. Both `001` and `1` must remain distinct; the unfiltered dataset has **50 rows**.
3. Add an **AFRLT / min-max / -2-2** filter and submit: **34 rows** remain. Change its value to **-3--1**: **18 rows** remain. Enter `bad`: an explicit validation message appears, and the session stays connected; a corrected range works again.
4. Save settings with a valid filter. Remove the filter, submit, and save settings without filters. Restore the filtered, unfiltered, and filtered settings in the same session. The row counts must match each restore, without a stuck modal or disconnection.
5. Supplemental API check for multi-value exclusion: `aNCA::apply_filters(data.frame(x = c(1,2,3,1,2,3)), list(list(column = "x", condition = "!=", value = c("1","2"))))` must return only the two `3` rows. Numeric-range exponent notation such as `1e-3-2e-3` is covered by automated tests.

The CSV fixture is needed to exercise lexical subject IDs through upload. Other tests create small synthetic data inline; no additional study data are included.

## Summary

Addresses review findings **01, 07, and 12**, which meet at the input/mapping/filtering boundary.

- CSV inference now retains text if numeric/logical inference would change its representation, including leading zeros, long identifiers, and custom source columns.
- Mapping converts known numeric measurement fields only after their roles are assigned, and rejects invalid numeric measurements instead of manufacturing missing observations.
- Filter removal destroys obsolete input observers. Restore uses monotonically increasing IDs and discards cleared modules before invoking their reactives; empty restores clear the saved applied-filter state too.
- `!=` excludes the complete selected set. Numeric ranges support signs, decimals, and exponents. Invalid conditions or ranges fail explicitly instead of silently returning unfiltered data, and app validation remains recoverable.

## Validation

- `devtools::document()` with roxygen2 **7.3.3** completed.
- Full `devtools::test(stop_on_failure = TRUE)`: **4,370 passed**, **0 failures**, **0 warnings**, **3 explicit skips** (installed-only namespace test, opt-in browser runner, full-session slide integration).
- Targeted reader/mapping/filter/restore tests: **203 passed**; the final CSV/RDS parity test run passed **114 reader expectations**.
- `lintr::lint_package()`: **0 findings**; final reader-file lint also passed.
- Final `devtools::check(document = FALSE, cran = FALSE, manual = FALSE, args = c("--no-tests", "--no-manual"))`: **0 errors, 0 warnings, 0 notes**. Tests were run separately above; this is not a remote CRAN incoming check.
- Installed-app Chromium smoke verified upload/mapping (**50 rows**, IDs `001` and `1`), signed range (**34 rows**), an invalid-range notification without disconnection, correction (**18 rows**), and removal (**50 rows**).

Browser YAML replay remains a reviewer step above; repeated imported-filter state is covered by production-module `testServer()` regressions.

## Notes for reviewers

The existing `(NA)` sentinel policy is preserved: selecting the sentinel includes missing rows in the retained result. Invalid filter inputs now error through the public API instead of warning and returning the original data; the app presents that error without terminating the session.

Lossless CSV inference intentionally preserves numeric-looking text until mapping establishes its role. This also prevents values beyond exact double precision from collapsing before identifier mapping. Package version advances to `0.2.0.9001`; the sibling branches use the same development version and remain independent.
