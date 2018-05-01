if(!is.loaded("mpi_initialize")) {
  library(Rmpi)
}

library(MASS)
library(doParallel)
library(doMPI)
library(doRNG)
library(slackr)
library(xtable)

wd <- getwd()

source(paste0(wd, "/functions.R"))
source(paste0(wd, "/generateSMART.R"))
source(paste0(wd, "/simulateSMART.R"))

slackrSetup(config_file = paste0(wd, "/d3slack.dcf"))

### CLUSTER SETUP ###
# Check for Flux and create cluster appropriately
if (Sys.info()["sysname"] != "Windows"){
  clus <- startMPIcluster()
} else {
  ncore <- detectCores()
  clus <- makeCluster(ncore)
}

# Set up cluster for parallel computation

# clusterExport(clus, varlist = "wd")
# clusterEvalQ(clus, {
#   library(MASS)
#   library(doParallel)
#   library(doRNG)
#   library(slackr)
#   library(xtable)
#   
#   source(paste0(wd, "/functions.R"), echo = F)
#   source(paste0(wd, "/generateSMART.R"), echo = F)
#   source(paste0(wd, "/simulateSMART.R"), echo = F)
# })

if (Sys.info()["sysname"] != "Windows"){
  registerDoMPI(clus)
} else {
  registerDoParallel(clus)
}

if (Sys.info()["sysname"] != "Windows"){
  nWorkers <- clus$workerCount
} else {
  nWorkers <- length(clus)
}

