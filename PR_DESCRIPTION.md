# Keep matrix ratios and collection amounts invariant across equivalent units

## How to test

1. Launch this checkout with `devtools::load_all(); aNCA::run_app()`.
2. **Data → Upload:** choose `tests/manual/data/matrix-unit-equivalence.csv`.
   Keep the same-name mappings, apply no filters, and continue to Preview.
3. Open **NCA → Setup** and select **DrugA**, both specimens, and a dose
   profile if the controls are not already populated. Open
   **Additional Analysis → Concentration Ratios** and add a ratio pair:
   **BLOOD / PLASMA**. Every sample ratio and geometric mean must be **1**.
   Source values stay **1 ug/mL** and **1000 ng/mL** in the result table.
4. Reverse the pair: ratios must still be **1**. Download the table and verify
   the same values. Compare with the base branch, which reports **0.001** for
   BLOOD / PLASMA.
5. To check recovery, edit a copy of the CSV so PLASMA has `hour` as AVALU,
   a blank AVALU, or a zero concentration. Upload it and recreate the pair: a clear error
   notification appears, the app stays connected, and no previous ratio table
   remains available for export. Correct the data and retry.

For the collection-density correction, the smallest exact check is in R:

```r
d <- data.frame(AVAL = c(1, 1000), AVALU = c("mg/g", "mg/kg"),
                VOLUME = 1, VOLUMEU = "mL")
v <- aNCA::convert_volume_units(d)
v$AVAL * v$VOLUME  # both 1 mg, rather than 1 mg and 1000 mg
```

No separate density dataset is needed. The automatic fallback remains an
**assumption of 1 g/mL**, not a measured specimen density. Supply measured
collection mass/volume data when that assumption is inappropriate.

## Description

Addresses review findings **02, 04, and 16** together because all use the unit
conversion layer.

- Matrix ratios convert denominator concentrations to numerator units before
  division, retain original values/units, and reject missing/blank or incompatible units,
  infinite concentrations, and zero denominators. Missing pairs remain omitted.
- Ratio errors become recoverable app notifications and clear the stored export
  result. Duplicate grouping names are normalized before the pair join.
- Collection conversion uses one physical density, **1 g/mL**, consistently in
  both mass/volume directions and reports the assumption. Other incompatible
  dimensions no longer receive an invented conversion.
- Conversion factors are computed once per distinct unit pair within a call,
  then mapped back without an unbounded global cache. Integer-coded pair keys
  avoid collisions from unit text. The unused first summary conversion is removed.

The existing 3-significant-digit ratio presentation is preserved. This branch
changes computed values where previous results depended on unit notation; those
outputs should be recalculated.

## Validation

- Final focused scientific/API/module suite: **252 passed**, no failures or warnings.
- Chromium production-module check: ratios **1**, visible invalid-unit notification,
  export cache cleared, then successful recovery without disconnection. Module
  regressions also cover NA and blank unit labels with notification/cache recovery.
- CSV roundtrip builds a PKNCA object and produces **88 matching ratios of 1**.
- Full-package lint: **0 findings**; final touched-file lint after the unit-label
  guard: **0 findings**. Documentation generated with roxygen2 7.3.3.
- `devtools::check()` with `--no-tests --no-manual`, CRAN mode off:
  **0 errors, 0 warnings, 0 notes**, rerun after the final unit-label guard.
  Tests are validated separately below.
- 10,000 repeated unit pairs: baseline **15.162 s**, optimized **0.011 s**,
  identical output. These local timings were recorded under concurrent load.
- Full development suite before the final missing-unit guard: **4,364 passed**,
  **0 failures/warnings**, **3 skips** (browser opt-in, installed-only namespace,
  full-session slide integration). Final guard regressions are covered in the
  focused suite and final package check rather than rerunning that broad suite.

Quarto deck rendering was not validated. The environment's Quarto launcher
reports a missing Deno executable; this did not affect the completed package check.

## Contributor checklist

- [x] Value-level regressions cover equivalent and incompatible units, density,
  recycling, missing/identity conversions, and matrix calculation recovery.
- [x] Public documentation regenerated with roxygen2 7.3.3.
- [x] User-visible scientific corrections recorded in NEWS; version 0.2.0.9001.
- [x] Full-suite baseline and final focused/lint/package checks recorded above.

## Notes to reviewer

The CSV fixture is a synthetic matched-unit variant of bundled example data,
created only for the app check above. No source data or analysis is pushed.
Full-app upload steps are supplied for review. The production Concentration
Ratios module was exercised in Chromium: equivalent units produced ratios of 1,
invalid units showed a notification and cleared export state, and corrected units
recovered without disconnecting. Automated fixture/calculation and module tests
provide additional directly executed evidence. The performance claim is
limited to the recorded conversion benchmark, not an end-to-end app guarantee.
