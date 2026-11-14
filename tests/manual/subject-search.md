# Manual check: Match complete subject IDs in the slope selector

Base: `chore/cran-review-fixes`. Fix: `fix/cran-subject-search`.

Compare `chore/cran-review-fixes` with this fix branch in **fresh R sessions**. Keep downloaded CSVs/probes outside the checkout so switching branches does not remove them. Use CRAN PKNCA 0.12.1 unless this guide explicitly requests development PKNCA.

For app checks, run from the selected checkout:

```r
devtools::load_all(export_all = FALSE)
aNCA::run_app()
```

Use Data → Upload, then Next → Mapping → Next → Filtering → Next → Preview. Keep same-name column mappings and apply no filters. For the supplied CSVs, check the mappings named below before continuing.

**Data:** `subject-search.csv`: 18 invented serum samples for subjects `1`, `10`, and `SUB[01`. The default data do not contain these edge-case IDs.

1. Upload the CSV; check that USUBJID retains those exact strings. Map the data without filters.
2. In NCA Setup select FakeDrug, SERUM and DOSE 1. In Parameter Selection include LAMZHL (half-life), then open Slope Selector. Each subject has a valid 2-hour terminal half-life.
3. Select subject `1`, then select `SUB[01` alone. Repeat on the fix branch in a fresh R session.

**CRAN branch:** selecting `1` also matches `10`; selecting `SUB[01` raises a regular-expression error and can disconnect the session.

**Fix branch:** each selection shows exactly that subject, including `SUB[01`. Clear the subject selector to restore all three subjects' plots.

Files: [subject-search.csv](data/subject-search.csv).
