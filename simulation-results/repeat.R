rm(list = ls())

setwd("~/GitHub-Bios/prgmetasa-r/simulation-results")


##
## CALCULATE ALPHA
##

setwd("~/GitHub-Bios/prgmetasa-r/simulation-results/exp0.2/a-calc")
source("simu.repeat.a.R")

setwd("~/GitHub-Bios/prgmetasa-r/simulation-results/U14/a-calc")
source("simu.repeat.a.R")

##
## SIMULATE RESULTS
##

setwd("~/GitHub-Bios/prgmetasa-r/simulation-results/exp0.2")
source("1000-times-simulation.R")


setwd("~/GitHub-Bios/prgmetasa-r/simulation-results/U14")
source("1000-times-simulation.R")

