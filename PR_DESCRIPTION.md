# FIX: Preserve literal text and reproducible export values

## How to test

Use the bundled example data; no new dataset is needed.

1. Start the app with `devtools::load_all(); aNCA::run_app()`. Continue through Data → Mapping → Filtering → Preview.
2. In TLG, add an individual concentration plot. Enter `Study {A} {{B}} — $PARAM` as its title. Braces must display literally; `$PARAM` must still resolve to the analyte. Check `!PARAM` label references and multiline footnotes too.
3. Set the project name to `Study "A"\Site`. Select DrugA/SERUM and dose 1 in NCA Setup, run NCA, then export the presentation/QMD. Its YAML header must parse with that same project title; the existing newline-to-space rule remains.
4. Add a manual slope rule with a reason containing quotes/backslashes, save settings, and export the R script. Open/parse the script: the reason must remain unchanged and no syntax error should occur.
5. API check for dataset paths (no data file is required to generate/parse the script): call `get_settings_code()` with the existing `tests/testthat/data/test-settings.yaml` and a path such as `/tmp/study "A".csv` or `C:\temp\data.csv`. Parse the result and verify its `data_path` assignment reproduces the exact input string.

Automated regressions cover annotation text that previously evaluated R, typed missing values, named objects, non-syntactic column names, factors/dates, numeric precision, public script generation, and YAML/setup-path roundtrips.

## Description

Addresses review findings **09, 10, and 13** together because all concern user text crossing into generated expressions or documents.

- Escape literal annotation braces before introducing the supported column/label lookups, so entered text cannot become an arbitrary glue expression. Substituted values remain literal too.
- Replace handwritten R serialization with base deparsing, preserving quotes, backslashes, missing-value types, names, classes, and numeric precision. Disable `niceNames` because it can emit unsafe backslashes in argument names.
- Serialize presentation YAML through `yaml`, retain portable true/false booleans and the existing newline policy, and serialize setup paths/library names as R literals.
- Replace internal formatting-only serializer tests with reconstruct-and-compare assertions; retain the generated BLQ-rule test as a value-level check.
- Bump the development version to 0.2.0.9001 and record the user-visible changes in NEWS. Review/PR Markdown files stay outside built source packages.

## Validation

- Full development suite: 4,335 passing assertions, 0 failures, 0 warnings, 3 explicit skips (opt-in browser runner, installed-only namespace test, full-session export integration).
- Focused generated-script regression suite: 46 passing assertions, no failures or warnings.
- `lintr::lint_package()`: 0 findings; documentation regenerated with roxygen2 7.3.3.
- `R CMD check --no-tests --no-manual`: 0 errors, 0 warnings, 1 environment time-verification NOTE. Tests execute separately in the full development suite.
- Isolated installed-package API checks passed for quoted/Windows script paths, literal annotations, and quoted/backslash YAML titles.
- Installed app connected in Chromium, displayed upload controls, and had no Shiny output errors.

The complete manual app flows above have not been replayed in the browser. A full presentation render was not run: this environment's Quarto wrapper lacks its bundled `deno` executable; YAML headers and generated R setup paths were verified directly.

## Notes to reviewer

Generated scripts use standard `structure(...)` expressions for some named values instead of the previous custom pretty-print format. Their values are preserved; exact formatting is not a public contract. No additional runtime dependencies or test datasets were introduced. This is an independent branch based on `docs/cran-review-fixes-summary`; nothing was pushed.
