# automated Non-Compartmental Analysis (aNCA)
<img src='man/figures/aNCA_logo.png' align="right" alt="aNCA logo" width="200">

<!-- badges: start -->

[![Pharmaverse](https://pharmaverse.org/shields/aNCA.svg)](https://pharmaverse.org)
[![License](https://img.shields.io/badge/License-Apache_2.0-yellow.svg)](https://github.com/pharmaverse/aNCA?tab=Apache-2.0-1-ov-file)
[![CRAN status](https://www.r-pkg.org/badges/version/aNCA)](https://CRAN.R-project.org/package=aNCA)
[![R build status](https://github.com/pharmaverse/aNCA/actions/workflows/main.yml/badge.svg)](https://github.com/pharmaverse/aNCA/actions)
![Code coverage](https://raw.githubusercontent.com/pharmaverse/aNCA/badges/main/test-coverage.svg)
[![RStudio community](https://img.shields.io/badge/community-shiny-blue?style=social&logo=rstudio&logoColor=75AADB)](https://forum.posit.co/new-topic?category=shiny&tags=shiny)

<!-- badges: end -->

<br/>

> Our aim is to develop and share an open-source R Shiny application for performing Non-Compartmental Analysis (NCA) on clinical and non-clinical datasets worldwide and across pharmaceutical companies.
<br/>

## Description

This application enables users to upload their datasets and perform Non-Compartment Analysis (NCA) on both pre-clinical and clinical datasets, with the results being easily visualizable. Designed with user-friendliness in mind, this app aims to make NCA accessible and straightforward for all scientists. Among the features it currently possess, the App can:

- **Customize half life calculation**: Either by rule settings definitions or performing manual in-plot adjustments
- **Define AUC intervals of interest**: Providing by default last and to-infinite calculations
- **Visualize data and results** with interactive boxplots, summary statistic tables and scatter plots
- **Produce PP and ADPP** dataset formats of the resulting parameters
- **Save your analysis settings** and reupload them later to keep on analysing!
- **Reproduce your NCA outside of the app**, using the R script created

## Installation

### From CRAN

The stable version of the package can be easily downloaded:

```R
install.packages("aNCA")
```

---

### From GitHub

We recommend using [pak](https://github.com/r-lib/pak) for package installation, along with all system dependencies. If you do not have `pak` available, you will need to set it up first:

```R
install.packages("pak")
```

then you can install [aNCA](.) by running:

```R
pak::pak("pharmaverse/aNCA")
```

in your R console.

---

### Via cloning the repository (for contributors)

Alternatively, you can set up the package by cloning the repository through your terminal/shell:

```bash
git clone https://github.com/pharmaverse/aNCA.git
```

and then loading it directly using [devtools](https://github.com/r-lib/devtools) in your IDE (e.g. RStudio) console:

```R
devtools::load_all()
```

## Quick start

To run the application, simply invoke:

```R
aNCA::run_app()
```

The testing data will be automatically loaded upon application startup. You can provide your own dataset in the **data** tab. Here you can also specify pre-processing filters.

In the **NCA** tab, start off by loading the pre-processed data using _Submit_ button. You will also need to choose dose number in the _Settings_. Then, you will be able to run the NCA analysis. From there, you can also specify different analysis options, like applying flag rule sets and selecting slopes.

After the setup is done and analysis is performed, you are free to explore the results in the **Outputs** tab. Application supports various customizable plots, as well as report exporting.

For more detailed instructions, check out the [Get Started](https://pharmaverse.github.io/aNCA/articles/aNCA.html) page.

### Diagnostics and recovery

The **Diagnostics** button in the header is available from every tab. If a data
or NCA operation fails, follow **Open diagnostics** in the notification, check
your inputs and settings, and retry. A failed NCA attempt clears its previous
results and export caches; it does not change the NCA formulas.

The panel shows a reference ID and up to 100 recent session events. **Download
diagnostic report** creates a JSON file in the browser, including after a
disconnection. It contains only events received before disconnection and is
lost when the page is reloaded or closed. Reports use fixed summaries, excluding
raw error text, datasets, study/subject identifiers, filenames, server paths,
and credentials. Detailed shared server logs are never offered for download.
Support can correlate a reference with the server's normal, access-controlled
logs, which retain original error details and can contain sensitive information.
Restrict server-log access to authorized administrators. File logs rotate at
5 MB with five files per run; startup retains up to
four files from previous runs. Browser reports can omit events outside a live
session or below INFO severity.

Recovery covers the main upload, settings upload, mapping, filtering, data
preparation, and NCA operations. Other unexpected failures retain Shiny's normal
handling and get a session diagnostic when its error hook runs. `req()` and
`validate()` remain normal validation flow. Process termination, out-of-memory
failures, and native-code crashes can still require hosting-level recovery;
browser reports cannot reconstruct events they never received.

## Contributing

Check our [contributing guidelines](https://github.com/pharmaverse/aNCA/blob/main/CONTRIBUTING.md).


## Documentation and references

- Go to [our Website](https://pharmaverse.github.io/aNCA/) for further information on the **aNCA app**.
- The main package used by the App is `PKNCA`. You can find more of it on its [GitHub](https://github.com/humanpred/pknca)
- aNCA has different system validations, you can consult them [here](https://pharmaverse.github.io/aNCA/articles/validations.html)
