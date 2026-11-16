# Fix ratios between different partial-AUC intervals

## How to test

1. Start the app from this branch with `devtools::load_all(); aNCA::run_app()`.
2. Data → Upload: choose `tests/manual/data/partial-auc-ratios.csv`. Keep same-name mappings, mark B as the metabolite, and apply no filters.
3. NCA → Setup: select A and B / SERUM. In Partial Interval Calculations add AUCINT for **0–2** and **2–4 hours**. Use linear trapezoidal integration.
4. In Ratios choose numerator **B / AUCINT_0-2**, denominator **A / AUCINT_2-4**, subject aggregation **no**, factor **1**, and output code `INTRATIO`. Run NCA.
5. Expect exactly one `INTRATIO` row, value **0.4** (= 4 / 10). A's partial AUCs are 2 and 10; B's are 4 and 40.
6. Repeat with both ranges 0–2: expect **2**. Save settings, export the R script, and verify the same result after replay.

## Description

Addresses review finding **03**. Numerator and denominator are now selected independently before reference averaging, so two ranges sharing `AUCINT` cannot pool together. Ratio descriptions retain both requested range labels. The public ratio API and ordinary non-interval calculation behavior are preserved.

The branch includes the small dataset needed to make the app result hand-checkable. Tests cover exact values, same/different ranges, multiple doses, dose-relative versus absolute boundaries, and subject aggregation.

## Validation

- Full development suite: **4,334 passed**, 0 failures, 0 warnings, 3 explicit skips.
- Focused ratio tests: **93 passed**, no failures/warnings/skips.
- Full-precision synthetic CSV passed the real PKNCA → CDISC → interval-ratio pipeline: partial AUCs **2/10** and **4/40**, final ratio **0.4**.
- Lint: **0 findings**. Documentation regenerated with roxygen2 **7.3.3**.
- Package check with `--no-tests --no-manual`, tests run separately: **0 errors, 0 warnings, 1 NOTE** (environment cannot verify current time).

Manual app steps above are acceptance instructions; a full browser replay was not performed. The environment's Quarto launcher lacks Deno, so slide rendering was not verified. This branch does not change slide rendering.
