# Preserve predose samples across time units and input ordering

## How to test

1. Start this branch's app: `devtools::load_all(); aNCA::run_app()`.
2. Upload `tests/manual/data/predose-hours.csv`, keep same-name mappings, and apply no filters.
3. Exploration → select DrugA / SERUM, Timescale **By Dose Profile**, profile **DOSE2**. The preceding concentration **80 ng/mL** must be present at **−2 hours**.
4. Upload `tests/manual/data/predose-minutes.csv` with the same choices. Expect the same preceding point at **−120 minutes**, not a missing point. The physical profile must be unchanged.
5. Reorder the rows within DOSE1 and repeat: the latest eligible sample remains **80**, not **70**. This order case is also tested directly because normal app preprocessing sorts input.

## Description

Addresses review finding **11**. Predose backfill selects the latest eligible sample by actual time within four physical hours of the next dose, instead of taking the last input row within four raw units. Samples after the next dose are excluded. Tied latest samples are retained rather than arbitrarily discarding a measurement.

The helper automatically uses RRLTU; legacy inputs without units continue to mean hours. An optional `time_unit` argument supports explicit units. Unknown/inconsistent time units are rejected instead of silently changing the physical window. Existing original observations and previous-dose predose duplication are preserved.

Exploration now performs backfill before converting absolute times to time since dose, and reassigned samples receive the correct dose-time metadata. This fixes the actual app plotting path as well as the helper. The plotting boundary retains data.frame indexing semantics for attribute-bearing PKNCA unit-column metadata, including mixed-unit dose-escalation exports.

The two small equivalent-unit datasets are needed for the app comparison. Regression tests cover hours/minutes/days, row order, the four-hour boundary, excluded earlier/postdose samples, legacy inputs, and invalid/mixed units.

## Validation

- Full development suite: **4,365 passed**, 0 failures, 0 warnings, 3 explicit skips.
- Focused helper/exploration/diagnostic/ZIP suite: **247 passed**, no failures or warnings, 1 explicit full-session integration skip.
- After the final explicit column-lookup qualification, **31 helper assertions** and touched-file lint passed; full package lint had **0 findings**.
- Documentation regenerated with roxygen2 **7.3.3**.
- Final package check with `--no-tests --no-manual`, full tests run separately: **0 errors, 0 warnings, 1 NOTE** (environment cannot verify current time).
- Chromium production-function harness rendered actual points **x = [-2, 0, 1]**, **y = [80, 100, 90]** in hours, then **x = [-120, 0, 60]** with identical concentrations in minutes. Shiny stayed connected with **0 output errors**.

The browser harness exercised the branch CSVs through `read_pk()` → PKNCA creation → production `exploration_individualplot()` → Plotly. It did not automate full-app upload/navigation; the terse steps above remain that reviewer acceptance procedure.

Full slide rendering was not verified because the existing Quarto launcher lacks its Deno helper. ZIP utility regressions and package examples passed.
