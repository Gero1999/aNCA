# Render concentration TLG pages on demand without truncating exports

## How to test

1. Start the app with the bundled example data. In Data, keep the default mapping,
   filter PCSPEC to SERUM and PARAM to DrugA, submit the filters, then open Preview.
   No custom dataset is needed.
2. Open **TLG → Order details → Submit Order Details**, then **Graphs → pkcg01 -
   Individual Linear**. Set **Entries per page = 2** and move Next → Previous.
   The eight subject groups must keep their order; revisiting the page should be
   immediate. Repeat with Individual Log / Side-by-Side and Combined Linear.
3. Open the graph options and change the title or grouping. The current page
   must refresh; navigating back must not show the old title or old grouping.
4. Choose **Export as ZIP**, select the concentration TLG graphs, and export
   HTML/PDF/PNG as available. Every selected group must be included, including
   pages never opened. Titles, points, and colors must agree with the app.

## Description

Addresses review finding **17**. Individual and combined concentration plots
(`pkcg01`, `pkcg02`, and their six linear/log/side-by-side catalog wrappers)
prepare their group list without constructing every Plotly widget. The module
resolves only the visible page and reuses a bounded cache of up to twenty group
widgets per data/options generation. Changing data, labels, or options replaces
that generation, so outdated widgets are not reused.

Normal programmatic calls remain eager. The additive internal `.lazy` flag is
used by the app only for functions explicitly declaring support. Group order,
plot templates, labels, and color inputs are still derived from the complete
input. Summary plots (`pkcg03`) and unrelated TLG types retain their current
rendering behavior.

Export uses an explicit all-group callable with the same active options. It
materializes every deferred group only when requested and does not keep an
unbounded full-study reactive cache after download. Deferred rendering retains
notification/error handling rather than letting a failed plot close the session.

## Validation

- Focused deferred-plot and TLG module regression suite: **238 expectations
  passed**, no failures, warnings, or skips. Covers construction counts, page
  reuse, cache eviction, data/title invalidation, eager/lazy plot equivalence,
  error recovery, and complete exports.
- Local 1,440-row / twenty-subject probe: eager twenty plots **18.866 s**;
  deferred group plan plus first five plots **2.362 s**; repeated cached page
  below the timer's resolution. Shared-machine load affects these timings;
  they are not a universal end-to-end speedup claim.
- Installed-package production-module smoke passed: two real Plotly entries on
  each of the first two pages, then all eight ordered groups in explicit export.
- Documentation regenerated with roxygen2 7.3.3; package lint found no issues.
- Package check with `--no-tests` passed installation, namespace loading,
  examples, code analysis, and documentation: **0 errors, 0 warnings, 1 NOTE**
  (environment could not verify current time).
- Complete development suite: **4,370 passing expectations, 0 failures,
  0 warnings, 3 skips** (opt-in browser test, installed-only namespace test,
  and the existing full-session export integration placeholder).
- Full browser interaction and multi-format exported files were not manually
  exercised. Quarto's shared CLI also reports a missing Deno executable after
  check; this branch does not change or validate slide rendering.

## Scope and review notes

No custom testing data were added: the bundled study supplies eight subjects.
The cache retains at most twenty widgets, plus currently referenced outputs;
requesting All or a complete export necessarily constructs the full requested
set temporarily. Other TLG families can adopt the deferred protocol separately.
Version is incremented to `0.2.0.9001`, and NEWS records the user-visible change.
