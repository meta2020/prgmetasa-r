# [Reproducible results with R codes] 

# Sensitivity analysis for PB in meta-analysis of prognostic studies


This folder contains reproducible R codes of simulation studies and re-analysis of the example data.

The following packages are used in the simulation or example data

- "mixmeta", "foreach", "parallel", "doSNOW", "doRNG", used in the simulations

- "latex2exp", "kableExtra"; 

If they are not installed, please install from R CRAN `install.packages("package_name")`.


## Folders

- [example](example/): `Ki67-2.Rmd` is the main analysis R codes for Application

- [exampledata](exampledata/): the data used in Application

- [simulation](simulation/): the simulation codes

	-	Step 1. calculate the alpha using `alpha-calc.R` (can be omit)

	-	Step 2. do simulation using `comp-npt.R` when distribution of censoring is correctly specified; do simulation using `comp-npt-mis.R` when distribution of censoring is not correctly specified

	-	Step 3. Create summarize table using `summary-kable.R` for summarizing SAUC; summarize table using `summary-kable-par.R` for summarizing parameters



