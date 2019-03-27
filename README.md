# rmSMARTSize
Companion code for "Sample size considerations for comparing dynamic treatment regimens in a sequential multiple-assignment randomized trial with a continuous longitudinal outcome"

## File Descriptions
- [sampleSize.R](sampleSize.R): Contains function to compute sample size for a SMART with a longitudinal outcome in which the primary aim is to compare two embedded DTRs.
- [functions.R](functions.R): Helper functions for simulations and parameter estimation.
- [generativeFunctions.R](generativeFunctions.R): Helper functions for data generation.
- [generateSMART.R](generateSMART.R): Main function used to generate data from a SMART.
- [simulateSMART.R](simulateSMART.R): Wrapper to generateSMART and helper functions to enable iterative data generation and parameter estimation.
- [init.R](init.R): Script to initialize R environment for simulations
- [resultFunctions.R](resultFunctions.R): Helper functions for compiling results and tables

The [Simulations](Simulations) folder contains scripts to perform simulations under all scenarios compiled in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). Simulations are designed to be run in parallel in an environment where all available cores can be dedicated to this task; as such, we recommend appropriately modifying [simulateSMART.R](simulateSMART.R) before running simulation scripts on, say, a laptop rather than a high-performance cluster.
