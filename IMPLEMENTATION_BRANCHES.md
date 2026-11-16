# Local implementation branches

All nine solution branches start from the reviewed branch at `ff86d834`.
They are local commits; nothing was pushed and no remote PR was opened.

## How to test

Choose a branch below and start with its **PR_DESCRIPTION.md → How to test** section.
The guides put short app steps first, followed by API checks where a direct value test is clearer.
They distinguish executed validation from manual acceptance steps.

From `/workspace/aNCA`:

```bash
git switch fix/review-input-filters
```

Use the corresponding branch name for another solution. Restart R when switching branches. Each branch contains its own PR_DESCRIPTION.md.

`BRANCH_REVIEW.md` records the historical review; this index and each branch's PR guide describe the implemented fixes.

## Branch map

| Branch | Review findings | Scope | Workspace PR guide |
| --- | --- | --- | --- |
| `fix/review-input-filters` | 01, 07, 12 | Lossless identifiers and filter correctness/lifecycle | [PR guide](/workspace/anca-fixes/input-filters/PR_DESCRIPTION.md) |
| `fix/review-unit-ratios` | 02, 04, 16 | Physical unit conversion, matrix ratios, and repeated-pair reuse | [PR guide](/workspace/anca-fixes/unit-ratios/PR_DESCRIPTION.md) |
| `fix/review-interval-ratios` | 03 | Independent partial-AUC numerator/reference intervals | [PR guide](/workspace/anca-fixes/interval-ratios/PR_DESCRIPTION.md) |
| `fix/review-geometric-cv` | 05 | Consistent log-variance geometric CV | [PR guide](/workspace/anca-fixes/geometric-cv/PR_DESCRIPTION.md) |
| `fix/review-session-results` | 06, 08, 14, 15 | Restored analysis state, stale exports, and editor lifecycle | [PR guide](/workspace/anca-fixes/session-results/PR_DESCRIPTION.md) |
| `fix/review-literal-serialization` | 09, 10, 13 | Literal annotations and faithful R/YAML serialization | [PR guide](/workspace/anca-fixes/literal-serialization/PR_DESCRIPTION.md) |
| `fix/review-predose` | 11 | Time-unit/order invariant backfill through plotting consumers | [PR guide](/workspace/anca-fixes/predose/PR_DESCRIPTION.md) |
| `fix/review-lazy-tlg` | 17 | Deferred concentration-plot pagination with complete exports | [PR guide](/workspace/anca-fixes/lazy-tlg/PR_DESCRIPTION.md) |
| `fix/review-optional-app` | 18 | Core NCA without mandatory application dependencies | [PR guide](/workspace/anca-fixes/optional-app/PR_DESCRIPTION.md) |

The restored-state and table-editor work is deliberately combined: stale table observers could replay old edits during restore. The temporary table-editors precursor is included in session-results rather than offered as a separate solution.

## Validation

All solution branches have version **0.2.0.9001** and an app-first PR guide. Package checks run with `--no-tests` where the full suite was executed separately; those checks do not count as another test run.

| Branch suffix | Full-suite passing assertions | Failures / warnings | Check errors / warnings / notes |
| --- | ---: | --- | --- |
| input-filters | 4,370 | 0 / 0 | 0 / 0 / 0 |
| unit-ratios | 4,364 | 0 / 0 | 0 / 0 / 0 |
| interval-ratios | 4,334 | 0 / 0 | 0 / 0 / 1 |
| geometric-cv | 4,334 | 0 / 0 | 0 / 0 / 1 |
| session-results | 4,401 | 0 / 0 | 0 / 0 / 0 |
| literal-serialization | 4,335 | 0 / 0 | 0 / 0 / 1 |
| predose | 4,365 | 0 / 0 | 0 / 0 / 1 |
| lazy-tlg | 4,370 | 0 / 0 | 0 / 0 / 1 |
| optional-app | 4,328 | 0 / 0 | Status OK in full and UI-absent libraries |

- Each normal full suite has **three explicit skips**, recorded in its PR guide. Passed counts vary because branches add different regressions.
- All branches passed package lint; final small follow-up changes also received focused/touched-file checks.
- The units branch's full suite preceded its final missing-label guard. That final guard and API/module behavior passed **252 focused assertions** and a final package check. Its CSV line endings were normalized afterward with identical parsed rows.
- The predose branch also passed **247 focused assertions**, including ZIP consumers, and a Chromium production-function harness with actual points **(-2 h, 80)** and **(-120 min, 80)**. A final explicit column lookup also passed 31 helper assertions and the final package check.
- Optional-app additionally passed **3,559 installed assertions** in a genuinely UI-absent test library, with **52 reported skips** for unavailable optional integrations. Nine UI packages remained absent; core NCA results matched the full library.
- Browser validation scopes vary: complete startup, selected upload/filter flows, or production-module/function harnesses. Each guide names the exact scope; manual acceptance steps are not presented as completed tests.
- One-NOTE checks report the environment's inability to verify the current time. The existing Quarto launcher lacks its Deno helper, so full slide rendering was not verified; YAML/setup serialization tests passed.

Validation ran on Debian 13 with R 4.5.3 and PKNCA 0.12.1; cross-platform release acceptance was not performed.

Logs and local validation artifacts are retained under `/workspace/.anca-setup/implementation/`, outside the repositories.

## Local commits

| Branch | Commit |
| --- | --- |
| `fix/review-geometric-cv` | `22ae78ae` |
| `fix/review-input-filters` | `74a17de8` |
| `fix/review-interval-ratios` | `c50d8df3` |
| `fix/review-lazy-tlg` | `7a42d4f4` |
| `fix/review-literal-serialization` | `3ed14fb7` |
| `fix/review-optional-app` | `095c8795` |
| `fix/review-predose` | `92e310ef` |
| `fix/review-session-results` | `19c259c0` |
| `fix/review-unit-ratios` | `867ef915` |
