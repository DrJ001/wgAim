# Speed Optimal Experimental Design R Package Development

## Project Structure
- R package source directory.
- `R/`: Contains function definitions and source code.
- `tests/testthat/`: Contains unit tests.
- `man/`: Auto-generated documentation (do not edit manually). Use roxygen
  documentation if you have it.
- `inst/`: Package assets, examples and installed files.

## Tool Locations
- **Quarto executable:** `C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe`
  - Not on the system PATH; must be called with the full path from a shell.
  - To render a vignette: `& "C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe" render <file.qmd>`
  - Always run from inside the `vignettes/` directory (or pass the full file path).

## Common Development Commands
- **Load All Functions:** `Rscript -e "devtools::load_all()"`
- **Run All Tests:** `Rscript -e "devtools::test()"`
- **Run Specific Test File:** `Rscript -e "testthat::test_file('tests/testthat/test-filename.R')"`
- **Update Documentation:** `Rscript -e "devtools::document()"`
- **Check Package:** `Rscript -e "devtools::check(args = '--no-manual')"`
- **Install Dependencies:** `Rscript -e "pak::pkg_install('.')"` or `Rscript -e "install.packages(c('pkg1', 'pkg2'))"`

## Code Style & Rules
- **Formatting:** Follow the [Tidyverse Style Guide](https://style.tidyverse.org/). Use 2-space indentation.
- **Functions:** Use implicit returns (avoid `return()` unless necessary for early exits).
- **Namespacing:** Use `pkg::function()` for external calls; list dependencies in `DESCRIPTION`.
- **Documentation:** Use [roxygen2](https://r-lib.org) tags (e.g., `@param`, `@return`, `@export`).
- **Tests:** Prefer `expect_snapshot()` for checking messages or complex outputs over `expect_message()`.

## Workflow Preferences
- Always run `devtools::document()` after modifying R/ files.
- Before suggesting a fix, run tests to verify current failures.
- If a function is no longer needed, remove it entirely instead of deprecating
  (active dev phase).
- The Quarto executable is in the local RStudio bin directory 
  
## Locally available R packages 
- ASReml-R, ASExtras4
