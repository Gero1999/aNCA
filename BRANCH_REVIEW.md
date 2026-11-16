# Branch review: bugs and substantial improvement opportunities

Reviewed `docs/cran-review-fixes-summary`, based on
`chore/cran-review-fixes`, at commit `aea6b5834ef1839626af9fd6f2f2f40b084154b5`.
This reviews the code present on this branch, including inherited code; it does
not imply that every finding was introduced by the CRAN fixes.

No application code was changed. The implementation sections below describe
proposed work, not completed fixes.

## Findings at a glance

| ID | Short description | Priority | Evidence |
| --- | --- | --- | --- |
| 01 | CSV inference can collapse distinct subject identifiers | High | Public API reproduction |
| 02 | Assumed density changes calculated amount when equivalent units change | High | Value-level reproduction |
| 03 | Partial-AUC ratios pool the wrong numerator and denominator intervals | High | Calculation-path reproduction |
| 04 | Matrix ratios divide values without converting concentration units | High | Public API reproduction |
| 05 | Descriptive statistics use the wrong geometric CV formula | High | Value-level reproduction |
| 06 | A second restore can retain the first session's manual slope rules | High | Production module reproduction |
| 07 | Clearing restored filters can call NULL modules and close the session | High | Production module reproduction |
| 08 | Failed or outdated excretion calculations leave old export results available | High | Production module reproduction |
| 09 | Annotation text is evaluated as R code | High on a shared server | Harmless execution probe |
| 10 | Generated scripts lose escaping, missing-value types, and precision | High | Installed-package reproduction |
| 11 | Predose duplication depends on time units and input row order | Medium | Public API reproduction |
| 12 | Exclusion filters recycle vectors; signed ranges are ignored | Medium | Public API reproduction |
| 13 | Quoted project titles produce invalid Quarto YAML | Medium | Header parsing reproduction |
| 14 | Resetting a TLG table leaves stale edits ready to reapply | Medium | Production module reproduction |
| 15 | Table edits repeatedly register additional observers | Medium | Measured observer accumulation |
| 16 | Repeated unit parsing is a large avoidable calculation cost | High-value optimization | Equivalent-output microbenchmark |
| 17 | TLG pagination builds all plots before selecting the visible page | High-value optimization | Source inspection and plot benchmark |
| 18 | Core NCA users must install the entire Shiny stack | Dependency opportunity | Source and dependency-graph analysis |

“High” prioritizes incorrect scientific results, stale results, broken restore
workflows, or execution of user text. Module reproductions used real production
server functions with `shiny::testServer()`; browser and full-study confirmation
are explicitly requested where still needed. Optimization timings are local
measurements, not promises about total application speed.

## 01. Preserve subject identifiers before CSV type inference can merge them

**Location:** [R/readers.R](R/readers.R#L54),
[R/apply_mapping.R](R/apply_mapping.R#L100).

`read.csv()` infers numeric columns. Subject identifiers `001` and `1` both
become integer `1` through the public `read_pk()` function. Mapping cannot
recover the original strings afterward. This can merge subjects during grouping,
duplicate detection, and NCA, rather than merely changing their display labels.

**How to test:** Write a CSV containing the following three lines and import it
with `aNCA::read_pk()`:

```text
USUBJID,AVAL
001,1
1,2
```

The current result has IDs `1, 1`; require the distinct character values `001`
and `1`. Compare equivalent CSV and RDS inputs, then include a custom source
column mapped to USUBJID and verify that two subject profiles remain distinct.

**How to implement:** Preserve identifier text at ingestion, before inference
destroys it. Account for custom mapping names, not just columns already named
USUBJID. Convert mapped numeric analysis fields deliberately and validate them.
Simply applying `as.character()` after reading the current CSV is insufficient;
forcing every column to character without adapting numeric conversion would
also break the existing workflow.

## 02. Replace the unit-dependent density assumption with an explicit physical conversion

**Location:** [R/conversions.R](R/conversions.R#L258).

When direct volume conversion fails, `convert_volume_units()` invents a density
of one in the requested target unit divided by the original unit. Equivalent
concentrations of `1 mg/g` and `1000 mg/kg`, each with `1 mL` of specimen, are
therefore assigned `1 g` and `1 kg` respectively. Their derived amounts become
`1 mg` and `1000 mg`. The documented density assumption does not justify
changing the physical result when only unit notation changes.

**How to test:** Run:

```r
d <- data.frame(PCSPEC = "FECES", AVAL = c(1, 1000),
                AVALU = c("mg/g", "mg/kg"), VOLUME = 1, VOLUMEU = "mL")
v <- aNCA::convert_volume_units(d)
v$AVAL * v$VOLUME  # currently 1 and 1000, both expressed in mg
```

Require equal physical amounts for equivalent input units. Add mass/volume
conversions, incompatible dimensions, missing densities, and an excretion
calculation that checks both recovered amount and fraction excreted.

**How to implement:** Require a measured/configured density, with explicit
physical units, when crossing volume and mass dimensions. If a default density
is scientifically approved, express it once in a canonical unit and convert it
properly. Otherwise reject the calculation with a useful message. Preserve the
assumption in exported settings and reproducible scripts for auditability.

## 03. Select numerator and reference intervals separately for partial-AUC ratios

**Location:** [R/ratio_calculations.R](R/ratio_calculations.R#L366),
[ratio application](R/ratio_calculations.R#L501),
[reference averaging](R/ratio_calculations.R#L170).

When both requested parameters have the same base code but different intervals,
the helper retains both intervals, removes their range identity from parameter
selection, and allows both into the ratio calculation. Reference values are
then averaged across the retained intervals.

**How to test:** Create reference analyte A with AUC0–2 = 2 and AUC2–4 = 10,
and test analyte B with AUC0–2 = 4 and AUC2–4 = 40, for one subject/dose/specimen.
Request `B AUCINT_0-2 / A AUCINT_2-4` through `calculate_ratio_app()`. The current
path produces two values, approximately `0.6667` and `6.6667`, instead of the
single requested value `4 / 10 = 0.4`. Extend the existing different-interval
test to assert exact values and row counts, not just that some rows exist.

**How to implement:** Filter numerator and reference rows independently by
their complete parameter-and-interval identities before joining. Carry role
information or separate tables through the calculation instead of pooling
identically named base parameters. Preserve interval provenance in the output.
Cover same-range, different-range, multiple-dose, and subject-aggregation cases.

## 04. Convert compatible concentration units before calculating matrix ratios

**Location:** [R/ratio_calculations.R](R/ratio_calculations.R#L58),
[app caller](inst/shiny/modules/tab_nca/non_nca_ratios.R#L167).

`multiple_matrix_ratios()` divides raw numeric concentrations even though it
also receives their unit columns. Equal concentrations can produce ratios that
are wrong by several orders of magnitude; downstream summaries use that Ratio.

**How to test:**

```r
d <- data.frame(USUBJID = "A", NFRLT = 1,
                MATRIX = c("BLOOD", "PLASMA"), CONC = c(1, 1000),
                UNITS = c("ug/mL", "ng/mL"))
aNCA::multiple_matrix_ratios(d, "MATRIX", "CONC", "UNITS",
                            spec1 = "BLOOD", spec2 = "PLASMA")
# Ratio is currently 0.001; expected 1.
```

Also test incompatible units, missing values, zero denominators, and the actual
additional-analysis UI.

**How to implement:** Reuse `get_conversion_factor()` to express the denominator
in numerator units before division. Reject or clearly flag incompatible
dimensions and undefined divisions. Keep original values/units for traceability
and round the displayed result rather than the calculation input.

## 05. Use the log-variance formula for geometric CV consistently

**Location:** [R/calculate_summary_stats.R](R/calculate_summary_stats.R#L101).

The descriptive-statistics path uses arithmetic SD divided by geometric mean.
The conventional lognormal geometric CV is
`100 * sqrt(exp(var(log(x))) - 1)`, which is already implemented in the TLG
summary helper in [R/utils-tlg.R](R/utils-tlg.R#L251). Consequently, two output
paths can report different geometric CVs for the same observations.

**How to test:** Summarize values `c(1, 2, 4)` as one CMAX group with consistent
units. `calculate_summary_stats()` reports `76.376`; the log-variance formula
gives `78.537` after rounding. Assert agreement between descriptive statistics,
TLG summaries, and exported tables before rounding, or normalize their display
precision first: the TLG helper currently rounds GeoCV to one decimal, while
this path uses three. Include NA, nonpositive, and singleton cases.

**How to implement:** Share one geometric-statistics definition between these
paths. Define which observations are eligible, handle insufficient positive
observations explicitly, and apply rounding at presentation. Add value-level
scientific assertions rather than only checking table structure.

## 06. Apply the new session's manual slope rules on every explicit restore

**Location:**
[manual_slopes_table.R](inst/shiny/modules/tab_nca/setup/manual_slopes_table.R#L59),
[NCA consumer](inst/shiny/modules/tab_nca/nca_setup.R#L125).

Imported rules are applied only while `manual_slopes()` is NULL. After Session
A supplies rules, a later Session B override is ignored. A production module
probe retained A's exclusion at time `4.48` after supplying B's exclusion at
`6.011` and changed PKNCA data. These rules flow into
`PKNCA_update_data_object(hl_adj_rules = ...)`, so a completed restore can retain
the wrong analysis choices.

**How to test:** Restore two sessions with deliberately different manual slope
rules in one app session. Compare the table, selected slope points, and rerun
NCA with a fresh session containing only B. Also restore B with no rules and
verify A's rules disappear. The state defect was reproduced with `testServer()`;
full browser replay and numerical half-life/AUC effects still need verification.

**How to implement:** Treat an explicit restore as a separate event with its
own generation/identity. Replace or clear imported rules on that event, while
preserving intentional edits on routine reactive data refreshes. Reconcile
rules with the new dataset and ensure selectors and settings exports use the
same committed state.

## 07. Remove cleared filter modules before invoking their reactive functions

**Location:**
[clearing filters](inst/shiny/modules/tab_data/data_filtering.R#L77),
[reminder](inst/shiny/modules/tab_data/data_filtering.R#L171),
[filter execution](inst/shiny/modules/tab_data/data_filtering.R#L226),
[exported filter list](inst/shiny/modules/tab_data/data_filtering.R#L254).

Clearing assigns NULL into existing `reactiveValues` keys. Consumers then call
every value as `x()` before discarding NULL results. Restoring an empty filter
list after an existing `SEX == M` filter produced “could not find function
'x'” in the reminder observer and closed the mock session. The same pattern
appears in multiple consumers.

**How to test:** Exercise filtered → unfiltered and filtered → differently
filtered restores, including fewer filters than before. Check that the browser
remains connected, the modal completes, row counts match the new settings, and
saved/exported filters contain only current entries. The production module
probe reproduced the crash with UI insertion/rendering stubbed; confirm browser
behavior as part of the fix.

**How to implement:** Drop NULL entries before invoking module reactives, or
store modules in a reactive named list that supports actual key removal. Reuse
one collector in all three consumers. Ensure removal also manages obsolete
observers and clears any previous applied-filter cache when no filters remain.

## 08. Invalidate excretion export results after failure or changed input data

**Location:**
[excretion_analysis.R](inst/shiny/modules/tab_nca/excretion_analysis.R#L264),
[export-cache observer](inst/shiny/modules/tab_nca/excretion_analysis.R#L291).

The calculation is bound only to Submit. Its export cache updates only on
successful non-NULL output. After a successful 32-row calculation, submitting
invalid collection times made `analysis_result()` NULL but left the old export
table unchanged. Changing upstream concentrations also left the prior
calculation cached until resubmission. ZIP/slide export can consume that table.

**How to test:** Calculate successfully, then submit invalid end times; assert
old results are removed from export availability. Next calculate dataset A,
load dataset B with doubled concentrations, and verify results become explicitly
unavailable/stale before resubmission and recompute correctly afterward.

**How to implement:** Associate results with dataset and committed calculation
settings. Invalidate both displayed and exported results when upstream data
changes or a submitted calculation fails. If retaining prior results for user
convenience, mark them as stale and prevent them being exported as current.
Keep displayed tables and export availability driven by one result state.

## 09. Interpret annotation tokens without evaluating arbitrary R expressions

**Location:** [R/utils.R](R/utils.R#L88),
[TLG text inputs](inst/shiny/modules/tab_tlg/tlg_option_text.R#L27),
[rendering](inst/shiny/modules/tab_tlg/tlg_module.R#L333).

`parse_annotation()` replaces `$COLUMN`/`!COLUMN` tokens, then passes the entire
user string to `glue()`. Existing `{...}` text is evaluated as an R expression.
A harmless options-only probe confirmed execution. Even an ordinary title
`Group {A}` errors because R tries to look up A. On a shared deployment where
users can edit or restore annotation text, this crosses from text input into
server-side code execution under the application account.

**How to test:** Require `parse_annotation(data.frame(x = 1), "Group {A}")` to
return literal text. A regression probe using
`{options(anca_review_annotation_probe=TRUE); "text"}` must not change that
option; restore its previous value after the test. Also retain the documented
`$COLUMN`, `!COLUMN`, newline, missing-column, and label behavior. Exercise text
entered in the UI and restored from settings.

**How to implement:** Resolve only the supported tokens with direct string and
column/attribute lookup. Do not send remaining text to a general expression
evaluator. Define literal escaping and HTML handling separately. Preserve useful
ERR behavior for unavailable tokens without granting expression execution.

## 10. Serialize generated R code with faithful literals, names, and missing values

**Location:** [R/get_session_code.R](R/get_session_code.R#L103),
[character serialization](R/get_session_code.R#L143).

Character values are wrapped with `sprintf('"%s"', obj)` without escaping.
The installed public `get_settings_code()` writes unparsable R for paths such
as `/tmp/study "A".csv` or Windows paths with backslashes. The same serializer
turns character NA into the literal string `"NA"`, and numeric formatting loses
precision. List and data-frame names also use incomplete manual quoting.

**How to test:** Generate a script from
`tests/testthat/data/test-settings.yaml` with each quoted/backslash path, then
parse it. Add literal roundtrip tests for `c("A", NA_character_)`, control
characters, non-syntactic names, named lists, realistic settings tables, and
high-precision numbers. The current numeric probe changed
`1.2345678901234567` to `1.23456789012346`; that small difference demonstrates
serialization loss, not a separately established clinical error.

**How to implement:** Use base R's deparsing/serialization facilities for
literals and names instead of maintaining custom escaping rules. Preserve NA
types, numeric precision, names, and supported classes. Keep readability as a
formatting step after correctness. Validate against an installed package as
well as development loading, and execute a representative generated analysis.

## 11. Make predose backfill invariant to time units and row ordering

**Location:** [R/dose_profile_duplicates.R](R/dose_profile_duplicates.R#L109).

The documented four-hour window is implemented as a raw `>= -4` comparison;
RRLTU is not consulted. An eligible sample two hours before the next dose is
retained in hours but rejected when the same times are expressed in minutes.
Separately, `slice_tail(n = 1)` selects the last input row, not necessarily the
latest sample time.

**How to test:** For dose 1 samples at 7 h and 8 h, and dose 2 at 10 h, require
the 8 h sample to appear as a -2 h predose sample for dose 2. Multiply all four
time columns by 60 and set RRLTU to minutes; the physical result must remain
the same. Shuffle rows so 8 h precedes 7 h; the chosen sample must still be 8 h.
Test days and the four-hour boundary too.

**How to implement:** Accept/derive time units and convert the physical window
into each profile's units, with validation for inconsistent units. Select the
maximum eligible sample time with a deterministic tie rule, rather than the last
physical row. The main PKNCA creation path currently sorts samples, limiting the
row-order defect there; the exported helper still lacks that guarantee. The
unit defect affects sorted inputs too, including profile/diagnostic consumers.

## 12. Make exclusion and range filters obey their declared selection

**Location:** [R/apply_filters.R](R/apply_filters.R#L104).

There are two independently reproduced errors. Without the NA sentinel, `!=`
compares vectors element by element and recycles the exclusion values. For
`x = c(1,2,3,1,2,3)`, excluding `c(1,2)` retains `3,1,2,3` instead of only
`3,3`. The range parser splits on every hyphen, so a signed range such as
`-2-2` is rejected and the original unfiltered data are returned with a warning.
Negative predose times are a normal domain input.

**How to test:** Assert exact retained row IDs for multi-value exclusions,
with/without NA and with vector lengths that do not divide the row count.
Filter `x = -3:3` to `-2-2`; require `-2:2`. Include negative-to-negative
limits, decimals, exponents, reversed limits, and malformed input. Cover the
public API, restored YAML filters, and displayed filtering feedback.

**How to implement:** Implement exclusion as set non-membership with an explicit
NA policy. Parse two signed numeric bounds with an unambiguous grammar, or
store the bounds separately and adapt the legacy text representation. Validate
unsupported conditions and malformed values before committing a filter; do not
allow users to mistake an ignored filter for an applied one.

## 13. Quote report titles using YAML serialization

**Location:** [R/quarto-utils.R](R/quarto-utils.R#L24),
[slide caller](inst/shiny/functions/zip-utils.R#L914).

`create_qmd_doc()` inserts the project/title directly between YAML double
quotes. A title such as `Study "A"` creates an invalid header, preventing a
valid presentation from being generated. Newline replacement alone does not
handle quotes and backslashes.

**How to test:** Generate QMD files with quotes, backslashes, colons, Unicode,
and newlines in the title. Parse the YAML header and assert the intended title
roundtrips, then render one representative slide deck with the Quarto toolchain.
Header parsing was reproduced here; rendering every variant was not performed.

**How to implement:** Serialize header values with the already-required `yaml`
package and preserve the existing newline policy deliberately. Check interpolated
file paths in setup code using the R literal strategy from finding 10. Keep
output format and execution settings unchanged while adding escaping.

## 14. Reset both the committed TLG table and its editable copy

**Location:**
[tlg_option_table.R](inst/shiny/modules/tab_tlg/tlg_option_table.R#L143).

Reset changes `output_table` but leaves `edits_table` unchanged. The editor
renders and confirms the latter. In a production module probe, reset restored
default A while the editable copy remained EDIT4; confirming without another
edit restored EDIT4 and silently undid the reset.

**How to test:** Edit → confirm → reset → reopen → confirm without editing.
Defaults must remain in both the editor and generated output. Repeat with
added/removed rows and cancellation.

**How to implement:** Reset both representations atomically, clear editor
selection/memory, and trigger a refresh. Define one committed/default/editing
state transition so Reset and Cancel cannot disagree about which table is
authoritative. Avoid reusing stale pending edits after data or schema changes.

## 15. Register editable-column observers once rather than on every table change

**Location:**
[manual_slopes_table.R](inst/shiny/modules/tab_nca/setup/manual_slopes_table.R#L195),
[tlg_option_table.R](inst/shiny/modules/tab_tlg/tlg_option_table.R#L110).

Both modules create inner `observeEvent()` handlers inside an observer that
depends on the table being edited. Old handlers remain registered. Reactlog
counts in the slope fixture grew `13 → 21 → 29 → 37` across three edits in
an eight-column table. A two-column TLG fixture grew `10 → 12 → 14 → 16 → 18`
across four confirmations. The retained handlers all listen to subsequent
edits, creating duplicate work; callback counts were not separately measured.

**How to test:** Count observer registrations and callback invocations while
performing 100 edits, confirmations, resets, and restores. Counts should remain
bounded by the current schema, and one edit should produce one state update.
Also test a changed grouping-column schema and session teardown.

**How to implement:** Register stable column handlers once during module setup.
When schemas genuinely change, retain observer handles and explicitly destroy
obsolete ones before adding replacements. Do not recreate observers simply
because cell values changed. This is a confirmed lifecycle defect; its overall
latency and memory effect on long sessions should be measured after correction.

## 16. Convert repeated unit pairs once per call instead of reparsing every row

**Location:** [R/conversions.R](R/conversions.R#L14),
[summary callers](R/calculate_summary_stats.R#L70),
[ratio callers](R/ratio_calculations.R#L153).

`get_conversion_factor()` parses units repeatedly for every row, even when
only a few unit pairs occur. A 10,000-row probe with five repeated pairs took
**5.181 seconds** in the current helper and **0.005 seconds** in a prototype
that converted distinct pairs once and mapped factors back; `identical()`
passed. The summary also calculates an initial `conv_factor` column that is
unused downstream before calculating its final conversion factors.

**How to test:** Compare outputs for valid, invalid, missing, identity, recycled,
named, and zero-length inputs. Benchmark full descriptive-statistics and ratio
workflows on representative studies, including mostly-unique unit pairs, not
just the repeated-pair microbenchmark.

**How to implement:** Build a collision-safe table of distinct source/target
unit pairs per call, evaluate each once, and restore original order and length.
Preserve fallback/error behavior. Remove the unused summary conversion after
checking its callers. Prefer bounded per-call reuse over an unbounded global
cache. The local Debian 13/R 4.5.3 timing demonstrates helper-level potential;
it is not an application-wide speedup estimate.

## 17. Construct only the current page's TLG plots, with an explicit full-export path

**Location:**
[tlg_module.R rendering](inst/shiny/modules/tab_tlg/tlg_module.R#L333),
[pagination](inst/shiny/modules/tab_tlg/tlg_module.R#L382),
[plot construction](R/g_pkcg.R#L235).

The complete TLG list is built before pagination slices the visible entries.
Thus displaying five subject plots can still construct hundreds of Plotly
widgets. A probe cloning the first example subject's SERUM/DrugA records into 20 subject IDs
(1,440 rows) took **3.109 seconds** for all 20 plots versus **0.902 seconds**
for the first five. No incremental-memory/RSS reduction was established.

**How to test:** Instrument plot-conversion calls for 100 groups with five
entries per page. Only the visible groups should construct widgets; revisiting
a cached page should reuse valid results. Change data/options and verify cache
invalidation. Assert complete group counts, ordering, titles, and values in
all-group ZIP/PDF/HTML exports.

**How to implement:** Separate cheap group indexing/counting from expensive
plot construction. Render/cache pages by dataset identity, grouping, options,
and labels, preserving consistent palettes and ordering across pages. Exports
must explicitly iterate over every group: the current module returns the full
list for export, so simply truncating it would lose outputs. Profile on realistic
studies before deciding cache size and invalidation granularity.

## 18. Make the Shiny stack optional for programmatic NCA users

**Location:** [R/imports-shiny.R](R/imports-shiny.R#L7),
[DESCRIPTION](DESCRIPTION), [startup checks](R/run_app.R#L95).

Dummy namespace imports make application-only packages mandatory even for
programmatic NCA. Candidates for an optional application dependency layer are
`bslib`, `htmltools`, `htmlwidgets`, `logger`, `reactable`, `reactable.extras`,
`shiny`, `shinycssloaders`, `shinyjs`, `shinyjqui`, `shinyWidgets`, and `zip`.
Some remain transitively required through plotting dependencies, so moving all
twelve direct declarations does not remove twelve installed packages.

In this environment's installed dependency graph, removing those required
edges reduced the recursive Depends/Imports/LinkingTo closure from **94 to 79**
packages. This includes base/recommended dependencies and depends on installed
versions; it is not a universal CRAN package count or a measured download-size
or startup-time improvement.

**How to test:** Install/load/use core NCA functions with the optional UI stack
absent, and run checks with and without suggested packages. Then install the
full app dependencies and test installed-app startup, upload, restore, and
exports. Audit transitive requirements and core plotting APIs before changing
any declaration.

**How to implement:** Move truly app-only dependencies to Suggests and remove
their dummy namespace imports. Simultaneously replace `run_app()`'s current
Imports-derived check with an explicit application requirement list and clear
installation guidance. Publish a convenient full-app installation command.
Keep actual core dependencies required. This is a product/packaging choice:
programmatic users benefit, while app users still need the complete stack.

## Validation performed and suggested implementation order

- The current development suite completed with **4,327 passing expectations,
  0 failures, 0 warnings, and 3 skips**. Skips were the opt-in browser smoke,
  installed-only namespace regression, and full-session export integration test.
- `lintr::lint_package()` completed with **0 findings**.
- Focused value/API probes covered identifiers, density, ratios, geometric CV,
  predose duplication, filtering, and annotation evaluation. Production module
  probes covered restores, excretion cache, reset behavior, and observer counts.
- Generated-script escaping was checked using an isolated installation of the
  reviewed commit. Quarto headers were parsed; full report rendering and full
  browser replay of the newly found restore defects were not performed.
- Local probe artifacts and available logs are retained under `/workspace/.anca-setup/review/`.
  That directory is outside this repository and is not a portable test suite;
  convert the relevant reproductions into repository tests when implementing.

Start with subject identity, units, ratio interval selection, and geometric
statistics because they can silently change scientific output. Address restore
state and stale exports next, alongside annotation evaluation for any shared
deployment. Then fix serialization, filters, and table lifecycle behavior.
Benchmark unit conversion and page rendering before and after their changes.
Consider dependency restructuring as separate packaging work.

A passing general suite does not invalidate the targeted counterexamples. Each
fix should add a regression asserting the correct values/state, then rerun the
affected installed-package, settings replay, generated-script, and export paths.
Scientific acceptance of changed calculations remains necessary.
