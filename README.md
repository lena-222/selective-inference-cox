# selective-inference-cox

This repository contains the code accompanying the manuscript

**“Statistical inference after variable selection in Cox models: A simulation study”**.

It provides reproducible implementations of simulation designs and evaluation pipelines
for post-selection inference in Cox proportional hazards models after Lasso-type variable
selection, as well as a real-data analysis based on the METABRIC breast cancer cohort.

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

## Repository structure (high level)

- `scripts/`  
  Setup scripts, data download helpers, and analysis entry points.
- `R/` (or project root `*.R`)  
  Core functions for simulation, fitting, and selective inference.
- `results/`  
  Output written by the run scripts (CSV summaries, RDS objects, meta information).
- `data_raw/`  
  Cached raw data (not tracked or optionally ignored via `.gitignore`).

## Setup (run once)

This project uses CRAN packages (and optionally Bioconductor, depending on your pipeline).
Install R (>= 4.2 recommended) and (on Windows) Rtools. Run the setup script once:

```r
source("scripts/00_setup_packages.R")
```

All package versions used for experiments are documented via sessionInfo() written to the
corresponding *_meta.txt files in the results folders.

## METABRIC data access (example)
The real-data analysis uses clinical METABRIC data distributed via cBioPortal:

Data source:
https://www.cbioportal.org/study/summary?id=brca_metabric

Because API calls can occasionally be slow or rate-limited, the repository supports
automatic download with local caching, and (optionally) a manual fallback.

### Option 1: automatic download via cBioPortal API

The loader downloads the clinical data once and caches it locally.

cache folder: `data_raw/cbioportal/`

main cached file: `data_raw/cbioportal/metabric_clinical.rds`

Run (from project root):

```{r}
source("scripts/load_metabric_data.R")
metabric_df <- load_metabric_data()             # uses cache if present
metabric_df <- load_metabric_data(force = TRUE) # refresh cache
```

### Option 2: manual download + local import

If the API download fails repeatedly, you can download the clinical data manually:

Go to the METABRIC study page on cBioPortal:
https://www.cbioportal.org/study/summary?id=brca_metabric

Use the Download function to export the clinical table (CSV/TSV).

Place the downloaded file into `data_raw/cbioportal/` (e.g., `metabric_clinical_manual.csv`).

Load it in R and create the same cached RDS:

```{r}
# Example manual import (adjust filename/sep depending on your export)
dir.create("data_raw/cbioportal", recursive = TRUE, showWarnings = FALSE)

clin <- read.csv("data_raw/cbioportal/metabric_clinical_manual.csv",
                 stringsAsFactors = TRUE, check.names = FALSE)

saveRDS(clin, "data_raw/cbioportal/metabric_clinical.rds")
```

After that, the project will use the cached metabric_clinical.rds as usual.

# How a simulation run works (end-to-end)

A typical run script follows the same high-level steps:

1. Parse command line arguments (e.g., n, p, n_sim, censoring, tuning rule).

2. Generate or sample a dataset

  - generic simulation designs: generate (time, status, X) from the chosen baseline hazard and covariate structure

  - METABRIC designs: simulate from a fitted Weibull-PH model or sample from the real cohort

3. Run the analysis pipeline for each Monte Carlo repetition

  - variable selection (Lasso / tuning strategy)

  - post-selection inference methods (split / PSI / debiased, etc.)

  - predictive metrics (IBS / C-index) and selection metrics

4. Aggregate metrics across repetitions.

5. Save outputs

  - `*_outs.rds`: full per-repetition results

  - `*_perf_beta.csv`: aggregated metrics (coverage/power/type I error etc., when truth is available)

  - `*_meta.txt`: run configuration + sessionInfo()

## Example: running a generic simulation

```
Rscript run_analyze_res2.R \
  --beta_type=realistic --n=200 --p=20 --n_sim=200 \
  --dist=weibull --rho=0.1 --target_censoring=0.3 \
  --lambda_choice=cv --alpha=0.10
```

## Example: running a METABRIC-based simulation

Depending on your setup, you may run one of:

  - `sim_design=metabric_weibph` (METABRIC-based generator, no “truth”, used within the simuulation study)

  - `sim_design=metabric_weibph_truth` (METABRIC generator with a defined beta target, used for "oracle" simulation)

  - `sim_design=metabric_real` (subsampling/bootstrapping from the real cohort, not used within the simulation study)

```
Rscript run_analyze_meta2.R \
  --sim_design=metabric_weibph_truth --n=200 --n_sim=200 \
  --target_censoring=0.1 --lambda_choice=fix --alpha=0.10
```

# Reproducibility

All simulation results reported in the manuscript and Supplementary Material can be
reproduced using the scripts provided in this repository. Random seeds are set to ensure
reproducibility.

Depending on the number of Monte Carlo repetitions and scenarios, reproducing the full
simulation study may require several days of computation.
