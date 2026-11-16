# Correct geometric CV in descriptive statistics

## How to test

1. Start this branch's app: `devtools::load_all(); aNCA::run_app()`.
2. Data → Upload: choose `tests/manual/data/geometric-cv.csv`; keep same-name mappings and apply no filters.
3. NCA → Setup: select DrugA / SERUM, select CMAX, and run NCA. The three subjects' CMAX values must be **1, 2, 4 ng/mL**.
4. NCA → Descriptive Statistics: pool subjects (do not group by USUBJID), display Geomean and Geocv. Expect **2** and **78.537%**.
5. Create a TLG PK summary for the same group: expect geometric mean **2** and geometric CV **78.5%** (one-decimal TLG formatting).
6. Download the summary/export ZIP and check the same values. Save/replay settings if comparing a previously reported analysis.

## Description

Addresses review finding **05**. Both summary paths use one unrounded log-variance definition: `100 * sqrt(exp(var(log(x))) - 1)`. Previously the descriptive table divided arithmetic SD by geometric mean, which can misstate variability. Previously reported descriptive geometric CVs should be recalculated.

Geometric statistics use finite positive observations. No eligible observations give NA; a singleton gives its geometric mean and NA geometric CV. Each output keeps its existing display precision. Other arithmetic statistics retain their existing behavior.

A small three-subject dataset is included because it makes the app result independently hand-checkable. Regression tests cover known values, cross-output agreement, missing/nonpositive/nonfinite inputs, singleton behavior, and the dataset's actual CMAX values.

## Validation

- Full development suite: **4,334 passed**, 0 failures, 0 warnings, 3 explicit skips.
- Focused summary/TLG tests: **102 passed**, no failures/warnings/skips.
- The manual CSV passed the actual PKNCA calculation and summary pipeline: CMAX **1/2/4**, geometric CV **78.537%**.
- Lint: **0 findings**. Documentation regenerated with roxygen2 **7.3.3**.
- Package check with `--no-tests --no-manual`, tests run separately: **0 errors, 0 warnings, 1 NOTE** (environment cannot verify current time).

Manual app steps above are acceptance instructions, not a completed browser replay. The environment's Quarto launcher lacks Deno; full slide rendering was not verified and is unchanged by this branch.
