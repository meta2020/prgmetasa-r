##
## CALCULATE ALPHA
##
setwd("~/2023pb-tsroc/simulation-results0207-2")

setwd("~/2023pb-tsroc/simulation-results0207-2/exp0.2/a-calc")
source("~/2023pb-tsroc/simulation-results0207-2/exp0.2/a-calc/simu.repeat.a.R")

setwd("~/2023pb-tsroc/simulation-results0207-2/U14/a-calc")
source("~/2023pb-tsroc/simulation-results0207-2/U14/a-calc/simu.repeat.a.R")

##
## SIMULATE RESULTS
##

setwd("~/2023pb-tsroc/simulation-results0207-2/exp0.2")
source("~/2023pb-tsroc/simulation-results0207-2/exp0.2/1000-times-simulation.R")


setwd("~/2023pb-tsroc/simulation-results0207-2/U14")
source("~/2023pb-tsroc/simulation-results0207-2/U14/1000-times-simulation.R")

