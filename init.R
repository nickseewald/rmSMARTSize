# init.R
# Initialization of simulation environment
#
# Copyright 2018 Nicholas J. Seewald
#
# This file is part of rmSMARTsize.
# 
# rmSMARTsize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rmSMARTsize is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with rmSMARTsize.  If not, see <https://www.gnu.org/licenses/>.

check.rmpi  <- require(Rmpi)
check.dompi <- require(doMPI)

require(MASS)
library(doParallel)
library(doRNG)
library(slackr)
library(xtable)
library(here)

source(here("functions.R"))
source(here("generativeFunctions.R"))
source(here("generateSMART.R"))
source(here("simulateSMART.R"))

# Configuration for slackr notifications -- replace with appropriate .dcf file
slackrSetup(config_file = here("njssimsslack.dcf"))


### CLUSTER SETUP ###
# Check for Rmpi and create cluster appropriately
if (check.dompi) {
  clus <- startMPIcluster(verbose = T, 
                          logdir = here("outFiles", "workers", ""))
} else {
  ncore <- detectCores()
  clus <- makeCluster(ncore, 
                      outfile = here("outFiles",
                                          paste0("log", Sys.Date(), ".txt")))
}

if (check.dompi) {
  Sys.sleep(15)
  registerDoMPI(clus)
  Sys.sleep(15)
  nWorkers <- clus$workerCount
} else {
  registerDoParallel(clus)
  nWorkers <- length(clus)
}

