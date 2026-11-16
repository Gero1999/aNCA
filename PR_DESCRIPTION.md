# FIX: Restore current slope rules and invalidate stale excretion results

## How to test

Use the bundled example; no external dataset is needed. Start with
`devtools::load_all(export_all = FALSE); aNCA::run_app()`.

1. **Repeated slope restore:** Data → Next through Mapping, Filtering, Preview.
   In NCA → Setup select DrugA/SERUM, then select LAMZHL in Parameter Selection.
   Open Slope Selector, find S1-01,
   double-click a point to exclude it, and download settings as **A**.
2. Change the excluded point and save settings as **B**. Upload A, then B in
   the same app session. The manual rule table and selected points must match
   the last uploaded file. Run NCA and compare with B loaded in a fresh session.
3. Remove all manual rules, save **C**, and restore A then C. C must clear A's
   rules. Ordinary settings changes after a manual edit must preserve that edit.
4. **Excretion retry:** Include URINE in NCA Setup and run NCA once to enable
   NCA exports. Open Additional Analysis → Excretion. Uncheck body-weight adjustment; select URINE, AEFRLT, RCAMINT,
   and Samples; Submit. A result table and excretion export must be available.
5. Change the end-time column to AFRLT and Submit. A warning must appear;
   the old table and excretion entry in Export as ZIP must disappear. Restore
   AEFRLT and Submit: current results return without reconnecting.
6. Change the NCA specimen/profile selection or restore another session after a
   successful excretion run. Until Excretion is submitted again, previous
   excretion results must not be displayed or offered for ZIP/slide export.

## Description

Addresses branch-review findings **06** and **08**. Repeated session restores
previously ignored imported slope rules once any rules existed, and failed or
outdated excretion calculations could leave successful older results available
for export.

- Each imported slope override is treated as an explicit restore: clear the
  previous rules, then apply the import once matching data are available.
  Empty imports clear rules; ordinary data/settings refreshes preserve edits.
- Excretion calculations retain the submitted dataset generation. Changed or
  removed input data invalidate display/export state, and unsuccessful submissions
  clear the previous export table. A successful resubmission restores both.
- Add production-module regressions covering repeated/empty/deferred imports,
  edit preservation, failure/retry, changed/removed data, and downstream cache state.
- Increment the development version and exclude review/PR documents from source
  package builds.

## Validation

- Pinned roxygen2 7.3.3 documentation generation passed.
- Focused production-module tests: 78 passing expectations, no failures/warnings/skips.
- Lint for changed source/test files passed.
- Source-app Chromium startup smoke: upload UI, live Shiny connection, no output errors.
- Intermediate broad/check runs were stopped before completion because the related
  table-editor fixes will be integrated and the combined branch checked again.

The browser steps above are reviewer instructions, not a claim that every
interactive step was executed automatically.

## Notes to reviewer

A nonmatching imported rule set remains pending until matching data arrive;
previous rules remain cleared. An intervening user edit takes precedence over
that pending import. Excretion selections remain submit-driven: changing a
selection alone does not run a calculation.

This is a local branch based on `docs/cran-review-fixes-summary` (`ff86d834`).
It does not include the independent observer-lifecycle fixes in the table-editor
branch, which also changes the manual-slope module. Nothing has been pushed.
