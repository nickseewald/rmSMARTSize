# rmSMARTSize
Companion code for "Sample size considerations for comparing dynamic treatment regimens in a sequential multiple-assignment randomized trial with a continuous longitudinal outcome"

## File Descriptions
- [functions.R](functions.R): Helper functions for simulations. Contains sample size function, as well as functions used for parameter estimation.
- [generativeFunctions.R](generativeFunctions.R): Helper functions for data generation.
- [generateSMART.R](generateSMART.R): Main function used to generate data from a SMART.
- [simulateSMART.R](simulateSMART.R): Wrapper to generateSMART and helper functions to enable iterative data generation and parameter estimation.
- [init.R](init.R): Script to initialize R environment for simulations
- [resultFunctions.R](resultFunctions.R): Helper functions for compiling results and tables

The [Simulations](Simulations) folder contains scripts to perform simulations under all scenarios compiled in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)).
