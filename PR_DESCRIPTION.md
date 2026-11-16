# Make the Shiny stack optional for programmatic NCA

## How to test

1. In a normal full-app library, start this branch with `devtools::load_all(); aNCA::run_app()`. With the bundled data, complete Data mapping → run NCA → view results and export a ZIP. Existing app behavior should work.
2. In a clean R environment without the app stack, install this checkout with `pak::local_install(dependencies = NA)`. Run `library(aNCA)` and a bundled-data NCA calculation. The app-only stack (including `shiny`, `reactable`, and `logger`) should be absent while calculation still works. Some moved declarations remain transitive core plotting dependencies.
3. Call `aNCA::run_app()` there: expect one clear message naming the missing UI packages and the explicit installation command. It must not install anything automatically.
4. Install the listed packages, restart R, and repeat the app flow.

## Description

Addresses review finding **18**. Application-only dependencies move from Imports to Suggests; their dummy namespace imports are removed. App startup checks a dedicated requirement list, and README explains the full-app installation command. Core calculations and core plotting dependencies stay required. Tests guard only the affected optional app/reporting groups; the remaining core suite still runs.

No test dataset is needed: the bundled study exercises both core and app paths. The dependency reduction benefits programmatic use; running the app still requires its complete UI stack. Transitive package counts vary by dependency versions.

## Validation

- Full development suite: **4,328 passing assertions**, no failures/warnings,
  and three existing skips (opt-in browser, installed-only namespace test, and
  full-session slide integration).
- Restricted installed suite: **3,559 passing assertions**, no failures/warnings,
  and **52 reported skips**, including 13 app-fixture skips outside test cases.
  The returned test-case data reports 39 skips; neither count is a passing test.
- Package checks in full and UI-absent libraries: **Status: OK**. Both used
  `--no-tests --no-manual --no-vignettes`; suites ran separately as above.
  The restricted check used `_R_CHECK_FORCE_SUGGESTS_=false` and identified
  unavailable suggested packages as INFO. Remote CRAN/Bioconductor index access
  was blocked; incoming checks and a PDF manual were not verified.
- Documentation regenerated with roxygen2 **7.3.3**; package lint passed with
  zero findings, and the subsequently changed reporting tests also passed lint.
- Installed full-app Chromium smoke passed: Shiny connected, the Upload screen
  rendered the bundled **1,308 rows**, and no output errors appeared. The full
  manual NCA/ZIP sequence above remains a reviewer test plan.

The minimal core smoke used a private R runtime and a single restricted library,
without deleting or modifying shared installations. `shiny`, `logger`,
`reactable`, `reactable.extras`, `shinycssloaders`, `shinyjs`, `shinyjqui`,
`shinyWidgets`, and `zip` were genuinely unavailable. NCA produced **16 rows**
from the bundled study and `run_app()` correctly listed those nine missing
packages without attempting installation.

For the restricted suite, compatible non-UI test/reporting dependencies were
added only when their required dependency closure excluded those same nine
packages. Their actual absence was verified again afterward, and the NCA
result was identical to the full-library reference. Packages such as htmltools
that remain transitive plotting requirements were retained; this is not a claim
that every moved declaration disappears from the dependency graph.

No additional dataset was created. All validation artifacts are outside the
checkout under `/workspace/.anca-setup/implementation/optional-app/`.
