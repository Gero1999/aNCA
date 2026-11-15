# Prepared pull request descriptions: CRAN app stability fixes

This document prepares six separate pull requests for existing, pushed fixes.
It does not apply or combine those fixes. No pull request was created or updated
while preparing this document: the GitHub API was blocked by the environment
proxy (`CONNECT tunnel failed, response 403`). Existing matching PRs could not
be checked; check for the same fork/head/base before publishing and update a
matching description rather than creating a duplicate.

**For every PR:**

- Base repository: `Gero1999/aNCA`
- Base branch: `chore/cran-review-fixes`
- Head repository: `Gero1999/aNCA`
- Head branch: the corresponding `fix/cran-*` branch below

Do not target `main` or an upstream repository. Do not merge these branches as
part of PR preparation, and do not include diagnostics-feature work.

The freshly fetched base remains
[`63822b1d55ce3ddd4f09b9d3902574ea7462dc9f`](https://github.com/Gero1999/aNCA/commit/63822b1d55ce3ddd4f09b9d3902574ea7462dc9f).
Every listed fix has that exact merge base. If the base later incorporates a
fix, use this historical commit for its before comparison.

Branch diffs, commits, regression-test source, NEWS and manual guides were
inspected for this document. Local historical validation JSON/logs were not
available. Historical test/check outcomes below are explicitly attributed to
the supplied evidence and checked-in guides; no fresh tests or CI success are
claimed. Manual GUI follow-ups remain listed for reviewers. Unchecked checklist
items distinguish unverified/not-applicable requirements from source coverage.

Synthetic CSV/probe links point to the corresponding branch in this fork.
Download them outside the checkout before changing revisions, and use fresh R
sessions. CRAN PKNCA **0.12.1** is the default test dependency; only the strict
`VOLUME` rejection specifically needs development PKNCA **0.12.1.9003**, commit
`be353efaed885402f651238a13de34db03e41c9c`.

## Branch index

| Topic | Head branch | Inspected head commit |
|---|---|---|
| restore | `fix/cran-restore` | [`63ed793a`](https://github.com/Gero1999/aNCA/commit/63ed793af598f79416e556eff73eda1e43b504d6) |
| excretion | `fix/cran-excretion` | [`768e27b3`](https://github.com/Gero1999/aNCA/commit/768e27b32216f9a427bce5ab7d34dfb02289ad4f) |
| subject-search | `fix/cran-subject-search` | [`9030622c`](https://github.com/Gero1999/aNCA/commit/9030622c7639cc82bed01dacb701796bcb1d7401) |
| interval-metadata | `fix/cran-interval-metadata` | [`70100cb6`](https://github.com/Gero1999/aNCA/commit/70100cb6e3a062a9c0170e120c2309af51b8d90f) |
| export-progress | `fix/cran-export-progress` | [`a5f5ecfc`](https://github.com/Gero1999/aNCA/commit/a5f5ecfc216a3ed4af993be5dd9b361831bdc631) |
| slope-handlers | `fix/cran-slope-handlers` | [`283e4716`](https://github.com/Gero1999/aNCA/commit/283e4716d79347488fa56adfc7060ad623156d5f) |

## 1. FIX: Complete repeated settings restores

**Title:** `FIX: Complete repeated settings restores`

**Head:** `Gero1999/aNCA:fix/cran-restore`

**Base:** `Gero1999/aNCA:chore/cran-review-fixes`

### Issue

No issue number supplied.

### Description

A second settings restore in the same app session can leave “Restoring session…” on screen: readiness remains TRUE, so completion is not signalled again. Reset readiness when each restore starts so every completed restore dismisses the loading modal.

### Definition of Done

Repeated restores signal completion independently, including a third restore; existing restore behavior is preserved. This PR contains only the existing restore fix, regression test, manual guide and folded saved-session NEWS entry.

### How to test

Use **CRAN PKNCA 0.12.1** and the default `adnca_example`; no CSV is required. Compare the base and fix in fresh R sessions. The current base is the tested baseline `63822b1d`; if the base later incorporates this fix, use that commit for the before comparison.

On each checkout, run `devtools::load_all(export_all = FALSE); aNCA::run_app()`.

1. From Data → Upload, proceed through Mapping → Filtering → Preview with same-name mappings and no filters.
2. In NCA → Setup → Download settings → Save, download settings **without running NCA**. Save two copies outside the checkout as `restore-a.yaml` and `restore-b.yaml`.
3. In the same app session, upload A via Data → Upload and choose Restore if prompted. Wait for completion, then upload and restore B.
4. **Before (`63822b1d`):** the first restore completes; the second can leave the loading modal visible. **After (`fix/cran-restore`):** both restores complete and close it. Restore a third time to check repeatability.

Optional manual follow-up: repeat with saved settings created after Run NCA to exercise repeated NCA auto-run. Full [manual guide](https://github.com/Gero1999/aNCA/blob/fix/cran-restore/tests/manual/restore.md).

### Contributor checklist

- [ ] Code passes lintr checks — no new lint run in PR preparation
- [ ] Code passes all unit tests — historical validation reported; no new run claimed
- [x] New logic covered by unit tests (`test-tab_data.R` exercises three restores and connection state)
- [x] New logic is documented in the manual guide
- [x] App changes are reflected in the existing saved-session NEWS entry
- [ ] Package version is incremented — unchanged in this existing fix
- [ ] R script works with the new implementation — no new script validation
- [ ] Settings upload works with the new implementation — GUI steps above remain a manual reviewer check
- [ ] If any `.scss` change was done, run `data-raw/compile_css.R` — not applicable
- [ ] If a package dependency was added/changed, run `data-raw/test_suggests_hidden.R` — not applicable

### Notes to reviewer

Fresh inspection confirms the fetched branch diff against `Gero1999/aNCA:chore/cran-review-fixes` is limited to this fix; no diagnostics-feature work is included. Historical validation supplied with this work reported full checks with **zero errors/warnings**. Remaining NOTEs concerned time verification and, in parallel runs, example timings; isolated example reruns completed below three seconds. Local validation artifacts were unavailable, so these historical results are preserved as reported, not represented as newly verified tests or CI. Manual GUI checks above have not been rerun during PR preparation.

## 2. FIX: Recover from invalid excretion inputs and calculation errors

**Title:** `FIX: Recover from invalid excretion inputs and calculation errors`

**Head:** `Gero1999/aNCA:fix/cran-excretion`

**Base:** `Gero1999/aNCA:chore/cran-review-fixes`

### Issue

No issue number supplied.

### Description

Invalid collection end times or body-weight units can cause excretion calculation errors and disconnect the session. Validate collection intervals and finite positive body weights/recognized mass units before calculation, and contain unexpected calculation errors so users can correct inputs and retry.

### Definition of Done

Invalid inputs produce a recoverable message; correction and resubmission return results without disconnecting or altering the original input data. The existing excretion calculation formulas are preserved.

### How to test

Use **CRAN PKNCA 0.12.1** in fresh R sessions for the base and fix. The current base is the tested baseline `63822b1d`; if it later incorporates this fix, compare against that historical commit.

Download the synthetic fixtures from this fork/branch and keep them outside the checkout:

- [excretion-invalid-weight-unit.csv](https://raw.githubusercontent.com/Gero1999/aNCA/refs/heads/fix/cran-excretion/tests/manual/data/excretion-invalid-weight-unit.csv)
- [excretion-valid.csv](https://raw.githubusercontent.com/Gero1999/aNCA/refs/heads/fix/cran-excretion/tests/manual/data/excretion-valid.csv)

On each checkout, run `devtools::load_all(export_all = FALSE); aNCA::run_app()`. Proceed Data → Upload → Mapping → Filtering → Preview, retaining same-name mappings and applying no filters.

**Collection-time failure and retry:** with default `adnca_example` mapped, open NCA → Additional Analysis → Excretion. Select URINE, uncheck Adjust for Body Weight, select RCAMINT and Samples, set the end-time column to AFRLT, and Submit. **Before:** calculation can error and disconnect. **After:** a collection-time warning appears and the session remains connected; change the end-time column to AEFRLT and resubmit to obtain results.

**Body-weight failure and retry:** upload the invalid CSV and confirm AEFRLT, VOLUME, VOLUMEU, WTBL and WTBLU mappings. In Excretion select URINE, AEFRLT, RCAMINT and FREXINT, **Profiles only**, with **Adjust for Body Weight enabled**, then Submit. **Before:** the weight-unit calculation errors and can disconnect. **After:** a recognized-mass-unit warning appears without disconnection. Upload the valid CSV, repeat these selections and Submit: expect **100 µg recovered and 1% excreted**. Four samples contain 1, 2, 3, 4 µg/mL × 10 mL; dose 5 mg/kg × weight 2 kg = 10 mg. The fixtures differ only in WTBLU (`unknown` versus `kg`); the default data lack these weight fields.

Full [manual guide](https://github.com/Gero1999/aNCA/blob/fix/cran-excretion/tests/manual/excretion.md). These GUI checks remain for reviewers; they were not rerun during PR preparation.

### Contributor checklist

- [ ] Code passes lintr checks — no new lint run in PR preparation
- [ ] Code passes all unit tests — historical validation reported; no new run claimed
- [x] New logic covered by unit tests (invalid end times; unknown, non-mass and missing weight units; injected PKNCA error; correction/retry and connection state)
- [x] New logic is documented in the manual guide and fixtures
- [x] App changes are reflected in the existing urine-output NEWS entry
- [ ] Package version is incremented — unchanged in this existing fix
- [ ] R script works with the new implementation — no new script validation
- [ ] Settings upload works with the new implementation — not changed
- [ ] If any `.scss` change was done, run `data-raw/compile_css.R` — not applicable
- [ ] If a package dependency was added/changed, run `data-raw/test_suggests_hidden.R` — not applicable

### Notes to reviewer

Fresh inspection confirms the fetched diff against `Gero1999/aNCA:chore/cran-review-fixes` is limited to excretion validation/recovery, its tests/fixtures, manual guide and NEWS; it contains no diagnostics-feature work. Historical validation supplied with this work reported full checks with **zero errors/warnings**. Remaining NOTEs concerned time verification and, in parallel runs, example timings; isolated example reruns completed below three seconds. Local validation artifacts were unavailable, so these are preserved historical results, not newly verified tests or CI.

## 3. FIX: Match complete subject IDs in the slope selector

**Title:** `FIX: Match complete subject IDs in the slope selector`

**Head:** `Gero1999/aNCA:fix/cran-subject-search`

**Base:** `Gero1999/aNCA:chore/cran-review-fixes`

### Issue

No issue number supplied.

### Description

Slope-selector subject search builds a regular expression from IDs: selecting `1` also matches `10`, and `SUB[01` can raise a regex error and disconnect. Store each plot's original subject ID as metadata and match complete IDs literally; handle empty matches without creating invalid page choices.

### Definition of Done

Each ID selects only its own plots, regex punctuation is treated literally, no matches remain safe, and clearing the selector restores all plots. Scientific slope/half-life formulas are unchanged.

### How to test

Use **CRAN PKNCA 0.12.1** in fresh R sessions for the base and fix. The current base is the tested baseline `63822b1d`; if it later incorporates this fix, use that commit for the before comparison.

Download [subject-search.csv](https://raw.githubusercontent.com/Gero1999/aNCA/refs/heads/fix/cran-subject-search/tests/manual/data/subject-search.csv) from this fork/branch and keep it outside the checkout. It contains 18 invented serum samples for IDs **`1`, `10`, and `SUB[01`**, each with a valid 2-hour terminal half-life. The default dataset does not contain these edge cases.

1. Run `devtools::load_all(export_all = FALSE); aNCA::run_app()` on the selected checkout.
2. Data → Upload the CSV, confirm USUBJID retains those exact strings, then proceed through Mapping → Filtering → Preview with same-name mappings and no filters.
3. In NCA Setup select **FakeDrug**, **SERUM**, and **DOSE 1**. Include **LAMZHL (half-life)** in Parameter Selection, then open Slope Selector.
4. Select subject `1`, then select `SUB[01` alone. **Before (`63822b1d`):** `1` also matches `10`; `SUB[01` raises a regex error and can disconnect. **After (`fix/cran-subject-search`):** each selection shows exactly that subject and the session remains usable.
5. Clear the selector: all three subjects' plots return.

Full [manual guide](https://github.com/Gero1999/aNCA/blob/fix/cran-subject-search/tests/manual/subject-search.md). GUI checks above remain for reviewers; they were not rerun during PR preparation.

### Contributor checklist

- [ ] Code passes lintr checks — no new lint run in PR preparation
- [ ] Code passes all unit tests — historical validation reported; no new run claimed
- [x] New logic covered by unit tests (original-ID metadata, complete/literal matching, punctuation, multiple selections, no matches, clearing, pagination and connection state)
- [x] New logic is documented in the manual guide and synthetic fixture
- [x] App changes are reflected in the existing half-life-selection NEWS entry
- [ ] Package version is incremented — unchanged in this existing fix
- [ ] R script works with the new implementation — no new script validation
- [ ] Settings upload works with the new implementation — not changed
- [ ] If any `.scss` change was done, run `data-raw/compile_css.R` — not applicable
- [ ] If a package dependency was added/changed, run `data-raw/test_suggests_hidden.R` — not applicable

### Notes to reviewer

Fresh inspection confirms the fetched diff against `Gero1999/aNCA:chore/cran-review-fixes` contains only subject-ID plot metadata/filtering, safe page choices, tests/fixture, manual guide and NEWS; no diagnostics-feature work is included. Historical validation supplied with this work reported full checks with **zero errors/warnings**. Remaining NOTEs concerned time verification and, in parallel runs, example timings; isolated example reruns completed below three seconds. Local validation artifacts were unavailable, so these historical results are preserved as reported, not represented as newly verified tests or CI.

## 4. FIX: Preserve NCA compatibility and recover from slope plot errors

**Title:** `FIX: Preserve NCA compatibility and recover from slope plot errors`

**Head:** `Gero1999/aNCA:fix/cran-interval-metadata`

**Base:** `Gero1999/aNCA:chore/cran-review-fixes`

### Issue

No issue number supplied.

### Description

Incidental `VOLUME` columns in the dose-interval table cause strict development PKNCA to reject NCA/slope calculations; an uncontained slope-generation error can disconnect the session and leave old plots. CRAN PKNCA also needs the existing non-dotted callback aliases exported for run-time lookup in installed packages.

Omit incidental interval `VOLUME` unless explicitly retained or used for concentration grouping, while preserving sample volumes. Contain slope update failures, clear plots/cached slope state and permit correction/retry. Export and document the existing volume, excretion-rate and fraction-excreted callback aliases without changing formulas or native PKNCA definitions.

### Definition of Done

- Serum/urine values match an explicitly retained-metadata reference in shared result columns; explicit retention remains supported.
- Slope failures clear stale plots without closing the session, and corrected inputs regenerate them.
- Callback lookup works with namespace-restricted development loading and a real installed package.

### How to test

**Data:** built-in `adnca_example`; no CSV needed. **Strict rejection requires development PKNCA 0.12.1.9003**, commit `be353efaed885402f651238a13de34db03e41c9c`. **CRAN PKNCA 0.12.1 does not reproduce that specific VOLUME rejection.** Restart R when changing PKNCA versions.

1. On before baseline `63822b1d`, run `devtools::load_all(export_all = FALSE); aNCA::run_app()`. Use default data: Data → Upload → Next → Mapping → Next → Filtering → Next → Preview; keep same-name mappings and apply no filters.
2. In NCA Setup choose DrugA and SERUM; include LAMZHL in Parameter Selection. Open Slope Selector, then Run NCA if connected. With the stated development PKNCA, slope generation or NCA reports that VOLUME is not allowed; the slope observer can disconnect.
3. Repeat on `fix/cran-interval-metadata` with identical data/dependency. Slopes and NCA calculate without that rejection; also check URINE parameters and verify concentration sample volumes remain present.
4. For deterministic forced-error/retry, download [test-interval-metadata.R](https://raw.githubusercontent.com/Gero1999/aNCA/refs/heads/fix/cran-interval-metadata/tests/testthat/test-interval-metadata.R) and follow the probe steps below. The guide records 24 passing assertions on the fix, including failed-plot clearing, corrected retry, sample volumes and numerical equivalence; the before revision fails the metadata and forced-error checks.
5. Repeat the probe with **CRAN PKNCA 0.12.1** for callback lookup. Install each selected aNCA revision, restart R, attach it with `library(aNCA)` and check urine volume/excretion-rate callbacks. Before can report `anca_calc_volpk` missing; the fix resolves callbacks. Default `load_all()` exporting internals can conceal this installed-package problem.

Stop the app, keep the downloaded probe outside the checkout, and use fresh R sessions on the before and fix revisions. Run:

```r
devtools::load_all(export_all = FALSE)
testthat::test_file(
  file.choose(),
  env = new.env(parent = asNamespace("aNCA")),
  stop_on_failure = TRUE
)
```

Choose the downloaded probe. It loads production code from the selected checkout, so failures on the before revision are expected.

Detailed [manual guide](https://github.com/Gero1999/aNCA/blob/fix/cran-interval-metadata/tests/manual/interval-metadata.md).

### Contributor checklist

- [ ] Code passes lintr checks — historical results only; not rerun here
- [ ] Code passes all unit tests — historical results only; not rerun here
- [x] New logic covered by the included regression probe
- [x] New logic is documented in the manual guide
- [x] App or package changes are reflected in the existing NEWS entry
- [ ] Package version is incremented — unchanged in this existing fix
- [ ] R script works with the new implementation — no new script validation
- [ ] Settings upload works with the new implementation — not changed by this fix
- [ ] If any `.scss` change was done, run `data-raw/compile_css.R` — not applicable
- [ ] If a package dependency was added/changed, run `data-raw/test_suggests_hidden.R` — not applicable

### Notes to reviewer

Historical validation supplied with this work reported full checks with **zero errors and zero warnings**. Remaining NOTEs concerned time verification and, in parallel runs, example timings; isolated example reruns finished below three seconds. These results are historical, not new test or CI runs. This PR preparation re-inspected the fetched commits, diffs, tests, NEWS and manual guide; local validation JSON/log files were unavailable.

Both repositories are `Gero1999/aNCA`; base is `chore/cran-review-fixes`. The tested before baseline is [`63822b1d`](https://github.com/Gero1999/aNCA/commit/63822b1d55ce3ddd4f09b9d3902574ea7462dc9f); use that historical commit if the base later incorporates this fix. No diagnostics-feature work or history rewrite is included. No new CI success is claimed.

## 5. FIX: Close ZIP export progress after success and failure

**Title:** `FIX: Close ZIP export progress after success and failure`

**Head:** `Gero1999/aNCA:fix/cran-export-progress`

**Base:** `Gero1999/aNCA:chore/cran-review-fixes`

### Issue

No issue number supplied.

### Description

ZIP creation leaves its progress notification open after success or an export-writer failure, so repeated exports can leave stale progress UI. Register `progress$close()` with `on.exit(..., add = TRUE)` immediately after creating the progress object, covering success, failure and retry without changing exported content.

### Definition of Done

- Every started ZIP export closes its progress indicator on exit.
- A writer failure leaves no partial ZIP, restores the working directory and allows a successful retry in the same session.

### How to test

**Data:** built-in `adnca_example`; no CSV needed. **PKNCA:** CRAN 0.12.1; no development-specific behavior is required.

1. On before baseline `63822b1d`, in a fresh R session run `devtools::load_all(export_all = FALSE); aNCA::run_app()`. Map default data through Data → Upload → Next → Mapping → Next → Filtering → Next → Preview, keeping same-name mappings and applying no filters.
2. Open Export as ZIP. Select only `extras` → `session_info`; deselect slides, plots, tables and other items. Click Export ZIP with Results, then Export in Ready to Export. Repeat once. Before can download ZIPs but leaves progress notifications open.
3. Repeat on `fix/cran-export-progress`: both ZIPs download and each progress indicator disappears.
4. For **failure/retry**, download the existing [test-zip-progress.R controlled probe](https://raw.githubusercontent.com/Gero1999/aNCA/refs/heads/fix/cran-export-progress/tests/testthat/test-zip-progress.R) and run it on both revisions using the commands below. It writes a tiny synthetic text file, forces export-writer failure after progress starts, then retries. An invalid CSV can fail too early to exercise this path. Before fails progress-close assertions; the guide records all 13 assertions passing on the fix, including one close per attempt and no partial ZIP after failure.

Stop the app, keep the downloaded probe outside the checkout, and use fresh R sessions on the before and fix revisions. Run:

```r
devtools::load_all(export_all = FALSE)
testthat::test_file(
  file.choose(),
  env = new.env(parent = asNamespace("aNCA")),
  stop_on_failure = TRUE
)
```

Choose the downloaded probe. It loads production code from the selected checkout, so failures on the before revision are expected.

Detailed [manual guide](https://github.com/Gero1999/aNCA/blob/fix/cran-export-progress/tests/manual/export-progress.md).

### Contributor checklist

- [ ] Code passes lintr checks — historical results only; not rerun here
- [ ] Code passes all unit tests — historical results only; not rerun here
- [x] New logic covered by the included regression probe
- [x] New logic is documented in the manual guide
- [x] App or package changes are reflected in the existing NEWS entry
- [ ] Package version is incremented — unchanged in this existing fix
- [ ] R script works with the new implementation — no new script validation
- [ ] Settings upload works with the new implementation — not changed by this fix
- [ ] If any `.scss` change was done, run `data-raw/compile_css.R` — not applicable
- [ ] If a package dependency was added/changed, run `data-raw/test_suggests_hidden.R` — not applicable

### Notes to reviewer

Historical validation supplied with this work reported full checks with **zero errors and zero warnings**. Remaining NOTEs concerned time verification and, in parallel runs, example timings; isolated example reruns finished below three seconds. These results are historical, not new test or CI runs. This PR preparation re-inspected the fetched commits, diffs, tests, NEWS and manual guide; local validation JSON/log files were unavailable.

Both repositories are `Gero1999/aNCA`; base is `chore/cran-review-fixes`. The tested before baseline is [`63822b1d`](https://github.com/Gero1999/aNCA/commit/63822b1d55ce3ddd4f09b9d3902574ea7462dc9f); use that historical commit if the base later incorporates this fix. No diagnostics-feature work or history rewrite is included. No new CI success is claimed.

## 6. FIX: Bound manual slope edit handler registration

**Title:** `FIX: Bound manual slope edit handler registration`

**Head:** `Gero1999/aNCA:fix/cran-slope-handlers`

**Base:** `Gero1999/aNCA:chore/cran-review-fixes`

### Issue

No issue number supplied.

### Description

Every manual slope-table change re-registers edit observers, accumulating duplicate handlers. Delayed edits can also index rows that have already been removed. Register once per editable column, add a handler only for a newly introduced column, and ignore edits whose table/column/row is no longer valid.

### Definition of Done

- Repeated edits keep handler registration bounded.
- Removed-row edits do not recreate data or fail; a new grouping column receives one handler.

### How to test

**Data:** built-in `adnca_example` for GUI smoke; the controlled probe builds its own tiny synthetic rule table. **PKNCA:** CRAN 0.12.1; no development-specific behavior is required.

1. Download the existing [test-manual_slopes_table.R controlled probe](https://raw.githubusercontent.com/Gero1999/aNCA/refs/heads/fix/cran-slope-handlers/tests/testthat/test-manual_slopes_table.R). Compare before baseline `63822b1d` with `fix/cran-slope-handlers` in fresh R sessions using the commands below.
2. The probe counts actual observer registrations across 20 range edits, clears the table, delivers a stale edit and introduces a grouping column. Before registrations grow, its bounded-count assertion fails, and a stale edit after clearing can error. On the fix, the guide records all 27 assertions passing: count stays constant during edits, stale edits are ignored, and the new column gets exactly one handler.
3. For GUI smoke on each revision, run `devtools::load_all(export_all = FALSE); aNCA::run_app()`. Map default data through Data → Upload → Next → Mapping → Next → Filtering → Next → Preview, retaining same-name mappings and applying no filters. Open NCA → Setup → Slope Selector, add a manual exclusion/selection, repeatedly edit its range/reason, remove it, and add another row. On the fix, edits continue to apply and stale events do not restore removed rows.

A short GUI session does **not** necessarily reproduce a visible freeze. Use the controlled registration-count probe for the deterministic before/after comparison.

Stop the app, keep the downloaded probe outside the checkout, and use fresh R sessions on the before and fix revisions. Run:

```r
devtools::load_all(export_all = FALSE)
testthat::test_file(
  file.choose(),
  env = new.env(parent = asNamespace("aNCA")),
  stop_on_failure = TRUE
)
```

Choose the downloaded probe. It loads production code from the selected checkout, so failures on the before revision are expected.

Detailed [manual guide](https://github.com/Gero1999/aNCA/blob/fix/cran-slope-handlers/tests/manual/slope-handlers.md).

### Contributor checklist

- [ ] Code passes lintr checks — historical results only; not rerun here
- [ ] Code passes all unit tests — historical results only; not rerun here
- [x] New logic covered by the included regression probe
- [x] New logic is documented in the manual guide
- [ ] App or package changes are reflected in NEWS — this existing branch has no NEWS change
- [ ] Package version is incremented — unchanged in this existing fix
- [ ] R script works with the new implementation — no new script validation
- [ ] Settings upload works with the new implementation — not changed by this fix
- [ ] If any `.scss` change was done, run `data-raw/compile_css.R` — not applicable
- [ ] If a package dependency was added/changed, run `data-raw/test_suggests_hidden.R` — not applicable

### Notes to reviewer

Historical validation supplied with this work reported full checks with **zero errors and zero warnings**. Remaining NOTEs concerned time verification and, in parallel runs, example timings; isolated example reruns finished below three seconds. These results are historical, not new test or CI runs. This PR preparation re-inspected the fetched commits, diffs, tests, NEWS and manual guide; local validation JSON/log files were unavailable.

Both repositories are `Gero1999/aNCA`; base is `chore/cran-review-fixes`. The tested before baseline is [`63822b1d`](https://github.com/Gero1999/aNCA/commit/63822b1d55ce3ddd4f09b9d3902574ea7462dc9f); use that historical commit if the base later incorporates this fix. No diagnostics-feature work or history rewrite is included. No new CI success is claimed.
