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

## Reproducibility

All simulation results, tables, and figures reported in the manuscript and Supplementary
Material can be reproduced using the scripts provided in this repository.
Random seeds are set to ensure reproducibility.

Details on how to run the simulations and regenerate the results are documented in the
corresponding scripts.
