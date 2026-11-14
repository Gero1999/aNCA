# Manual check: Close ZIP export progress after success and failure

Base: `chore/cran-review-fixes`. Fix: `fix/cran-export-progress`.

Compare `chore/cran-review-fixes` with this fix branch in **fresh R sessions**. Keep downloaded CSVs/probes outside the checkout so switching branches does not remove them. Use CRAN PKNCA 0.12.1 unless this guide explicitly requests development PKNCA.

For app checks, run from the selected checkout:

```r
devtools::load_all(export_all = FALSE)
aNCA::run_app()
```

Use Data → Upload, then Next → Mapping → Next → Filtering → Next → Preview. Keep same-name column mappings and apply no filters. For the supplied CSVs, check the mappings named below before continuing.

**Data:** default `adnca_example`; no CSV needed.

1. Map the default data and open Export as ZIP. Select only `extras` → `session_info`; deselect slides, plots, tables and other items.
2. Click Export ZIP with Results, then Export in the Ready to Export dialog. Repeat once.

**CRAN branch:** the ZIP can download, but a progress notification remains after each export.

**Fix branch:** each ZIP downloads and its progress notification disappears.

**Failure and retry:** use the supplied `test-zip-progress.R` probe below. It creates a tiny synthetic text file, deliberately fails the export writer, then retries. On the CRAN branch its progress-close assertions fail; on the fix all 13 assertions pass. This controlled failure occurs after progress starts; an invalid CSV may be rejected earlier and would not reliably test this path.
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

Files: [test-zip-progress.R](../testthat/test-zip-progress.R).
