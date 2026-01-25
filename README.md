# selective-inference-cox

This repository contains the code accompanying the manuscript

**“Statistical inference after variable selection in Cox models: A simulation study”**.

It provides fully reproducible implementations of the simulation designs and evaluation
pipelines used to study post-selection inference in Cox proportional hazards models after
Lasso-type variable selection.

## Scope

The code covers inference procedures including:
- sample splitting,
- exact post-selection inference (PSI),
- debiased Lasso–based inference.

Performance is evaluated in terms of:
- selective coverage,
- selective confidence interval width,
- selective power and type I error,
- predictive performance measured by the integrated Brier score (IBS) and the
  concordance index (C-index),
- selection-related quantities such as average model size and the proportion of truly
  active variables among the selected covariates.

In addition, the repository contains code for the real-data analysis based on the
METABRIC breast cancer cohort.

## Setup (run once)

This project needs CRAN + Bioconductor packages. Install R (>= 4.2 recommended) and Rtools (Windows).
Run the setup script once:

`source("scripts/00_setup_packages.R")`

Data download (automatic, cached)

Clinical METABRIC data is pulled from cBioPortal and cached locally:

cache folder: data_raw/cbioportal/

main cached file: `metabric_clinical.rds`

To refresh the cache, set force = TRUE in the data loader call.
## Reproducibility

All simulation results reported in the manuscript and Supplementary
Material can be reproduced using the scripts provided in this repository.
Random seeds are set to ensure reproducibility.

The code is written in R and relies on standard packages for survival analysis,
penalized regression, parallel computation, and selective inference, including
survival, glmnet, selectiveInference, and related dependencies.
All package versions used in the experiments are documented via `sessionInfo()`.


Depending on the number of Monte Carlo repetitions, reproducing the full simulation
study may require several hours of computation.

The code is written in R and relies on standard survival analysis and selective inference
packages.
