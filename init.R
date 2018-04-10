library(MASS)
library(doParallel)
library(doRNG)
library(slackr)
library(xtable)

wd <- getwd()

source(paste0(wd, "/functions.R"))
source(paste0(wd, "/generateSMART.R"))
source(paste0(wd, "/simulateSMART.R"))

slackrSetup(config_file = paste0(wd, "/d3slack.dcf"))

### CLUSTER SETUP ###
nodefile <- Sys.getenv("PBS_NODEFILE")
# Check for Flux and create cluster appropriately
if (nodefile != "") {
  hostlist <- read.table(nodefile, skip = 1, header = FALSE)
  clus <- makeCluster(c(as.character(hostlist$V1)))
} else {
  ncore <- detectCores()
  clus <- makeCluster(ncore)
}

# Set up cluster for parallel computation

clusterExport(clus, varlist = "wd")
clusterEvalQ(clus, {
  library(MASS)
  library(doParallel)
  library(doRNG)
  library(slackr)
  library(xtable)
  
  source(paste0(wd, "/functions.R"), echo = F)
  source(paste0(wd, "/generateSMART.R"), echo = F)
  source(paste0(wd, "/simulateSMART.R"), echo = F)
})
registerDoParallel(clus)