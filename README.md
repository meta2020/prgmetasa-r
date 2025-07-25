# [Reproducible results with R codes] 

# Sensitivity analysis for PB in meta-analysis of prognostic studies

## Ciatation

Zhou, Y., Huang, A., & Hattori, S. (2025). Sensitivity analysis for reporting bias on the time-dependent summary receiver operating characteristics curve in meta-analysis of prognosis studies with time-to-event outcomes. Research Synthesis Methods, 16(3), 528–549. [https://doi.org/10.1017/rsm.2025.14](https://doi.org/10.1017/rsm.2025.14)

This folder contains reproducible R codes of simulation studies and re-analysis of the example data.

The following packages are used in the simulation or example data

- "mixmeta", "foreach", "parallel", "doSNOW", "doRNG", used in the simulations

- "latex2exp", "kableExtra"; 

If they are not installed, please install from R CRAN `install.packages("package_name")`.


## Folders

- [example](example/): `Ki67-2.Rmd` is the main analysis R codes for Application

- [exampledata](exampledata/): the data used in Application

- [simulation-new](simulation-new/): the simulation codes and results

	- [main-simulation](simulation-new/main-simulation/): the simulation codes

		-	Step 1. calculate the alpha using `alpha-calc.R` (can be omit)

		-	Step 2. do simulation using `main-simulation-comp-npt.R` when distribution of censoring is correctly specified; do simulation using `main-simulation-comp-npt-mis.R` when distribution of censoring is not correctly specified

	- [5npt](simulation-new/5npt/), [7npt](simulation-new/5npt/), [nptmis-LN](simulation-new/nptmis-LN/), [nptmis-UNIF](simulation-new/nptmis-UNIF/): results are saved in these folders

	-	Step 3. Create summarize table in Latex using `summary-kable.R` for summarizing SAUC; summarize table using `summary-kable-par.R` for summarizing parameters; or, print results in HTML using `summary-kable-print.Rmd` and `summary-kable-par-print.Rmd`.



