# FIX: Keep restored rules, edited tables, and exported results current

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
   Repeat a restore after typing into a rule's Reason cell; the old typed value
   must not overwrite the imported rule.
4. **Excretion retry:** Include URINE in NCA Setup and run NCA once to enable
   NCA exports. Open Additional Analysis → Excretion. Uncheck body-weight
   adjustment; select URINE, AEFRLT, RCAMINT, and Samples; Submit. A result table
   and excretion export must be available.
5. Change the end-time column to AFRLT and Submit. A warning must appear;
   the old table and excretion entry in Export as ZIP must disappear. Restore
   AEFRLT and Submit: current results return without reconnecting.
6. Change the NCA specimen/profile selection or restore another session after a
   successful excretion run. Until Excretion is submitted again, previous
   excretion results must not be displayed or offered for ZIP/slide export.
7. **TLG Reset:** Open TLG, add `pkcl01 – Concentration Listing (first dose)`,
   submit, and open Variables formatting. Change the USUBJID label, confirm,
   Reset to defaults, reopen, and confirm without editing. The default
   **Unique Subject Identifier** must remain; old edits must not return.
8. Repeat cell edits, confirmations, Cancel, row additions/removals, and manual
   slope edits. Each edit applies once and the app remains responsive. Cancel
   discards pending changes; Reset discards them even after reopening.

## Description

Addresses branch-review findings **06, 08, 14, and 15** together because restored
manual rules and editable-table observer lifecycles share state. Repeated session
restores previously ignored imported slope rules once any rules existed, and failed or
outdated excretion calculations could leave successful older results available
for export.

- Each imported slope override is treated as an explicit restore: clear the
  previous rules, then apply the import once matching data are available.
  Empty imports clear rules; ordinary data/settings refreshes preserve edits.
- Excretion calculations retain the submitted dataset generation. Changed or
  removed input data invalidate display/export state, and unsuccessful submissions
  clear the previous export table. A successful resubmission restores both.
- TLG Reset synchronizes committed and pending table values. Shared cell-edit
  handling retains one observer per current column, destroys removed/session-ended
  handlers, ignores malformed/stale edits, and avoids replaying retained browser
  input after a schema is recreated.
- Add production-module regressions covering repeated/empty/deferred imports,
  edit preservation, failure/retry, changed/removed data, and downstream cache state.
- Increment the development version and exclude review/PR documents from source
  package builds.

## Validation

- Documentation generation passed with pinned roxygen2 7.3.3.
- Combined focused production-module tests: **134 passing expectations**, no
  failures, warnings, or skips, including restore after a typed cell edit.
- Chromium source-app startup smoke passed: upload UI, live connection, no output
  errors. A production-module browser harness also passed Reset/Cancel, ten typed
  slope edits, and subsequent restores of different, empty, and original rules.
- Package check passed: **0 errors, 0 warnings, 0 notes** (`cran = FALSE`,
  `--no-tests --no-build-vignettes --no-manual`); tests are run separately.
- Final full development suite: **4,401 passing expectations**, no failures or
  warnings; three skips (opt-in browser smoke, installed-only namespace regression,
  and full-session export integration test).
- Final `lintr::lint_package()`: **0 findings**.

The browser steps above are reviewer instructions; the automated browser coverage
uses a production-module harness, not the complete saved-YAML replay in steps 1–3.
Full presentation rendering was not checked: the existing Quarto CLI installation
reports a missing Deno helper.

## Notes to reviewer

A nonmatching imported rule set remains pending until matching data arrive;
previous rules remain cleared. An intervening user edit takes precedence over
that pending import. Excretion selections remain submit-driven: changing a
selection alone does not run a calculation.

This is a local branch based on `docs/cran-review-fixes-summary` (`ff86d834`).
The related table-editor precursor is included here; no sibling branch is
required for these four findings. Nothing has been pushed.
