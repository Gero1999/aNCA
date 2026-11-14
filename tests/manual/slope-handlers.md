# Manual check: Bound manual slope edit handler registration

Base: `chore/cran-review-fixes`. Fix: `fix/cran-slope-handlers`.

Compare `chore/cran-review-fixes` with this fix branch in **fresh R sessions**. Keep downloaded CSVs/probes outside the checkout so switching branches does not remove them. Use CRAN PKNCA 0.12.1 unless this guide explicitly requests development PKNCA.

For app checks, run from the selected checkout:

```r
devtools::load_all(export_all = FALSE)
aNCA::run_app()
```

Use Data → Upload, then Next → Mapping → Next → Filtering → Next → Preview. Keep same-name column mappings and apply no filters. For the supplied CSVs, check the mappings named below before continuing.

**Data:** default `adnca_example`; no CSV needed.

**UI smoke check:** map the default data, open NCA → Setup → Slope Selector, add a manual exclusion/selection, edit its range/reason repeatedly, remove it, then add another row. On the fix, edits continue to apply and stale events do not recreate removed rows.

The original leak does **not** guarantee a visible freeze in a short manual session. To compare branches deterministically, run the supplied `test-manual_slopes_table.R` probe below. It creates a synthetic rule table internally and counts actual event-observer registration across 20 edits, clearing, and a new grouping column.

**CRAN branch:** observer registrations grow after edits; the bounded-count assertion fails, and a stale edit can error after clearing the table.

**Fix branch:** all 27 assertions pass; registration stays constant during edits, stale events are ignored, and one new column receives one handler.
Stop the app first. In a fresh R session in the selected checkout, run this command and choose the supplied probe file from the downloaded bundle:

```r
devtools::load_all(export_all = FALSE)
testthat::test_file(
  file.choose(),
  env = new.env(parent = asNamespace("aNCA")),
  stop_on_failure = TRUE
)
```

The same probe loads production code from whichever branch is checked out. Failures on the CRAN branch are expected.

Files: [test-manual_slopes_table.R](../testthat/test-manual_slopes_table.R).
