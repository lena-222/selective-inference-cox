# Statistical inference after variable selection in Cox models: a neutral simulation study

This repository contains the R code accompanying the published simulation study:

> Schemet L, Friedrich-Welz S. **Statistical inference after variable selection in Cox models: a neutral simulation study.** *BMC Medical Research Methodology*. 2026;26:143. https://doi.org/10.1186/s12874-026-02887-0

The study compares inference procedures after variable selection in Cox regression, including sample splitting, exact post-selection inference, and debiased-Lasso-based inference. The simulation code covers different coefficient patterns, censoring settings, correlation structures, tuning rules, and both toy and METABRIC-informed scenarios.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── docs/
│   └── PORTABILITY_CHANGES.md
└── code/
    ├── analyze_all_meta.R
    ├── analyze_toy_sim.R
    ├── run_analyze_meta.R
    ├── run_analyze_toy.R
    ├── example.R
    ├── config/
    ├── helpers/
    ├── scripts/
    ├── data_raw/
    └── results/
```

## Main entry points

- `code/run_analyze_toy.R`: generic/toy simulation settings.
- `code/run_analyze_meta.R`: METABRIC-informed simulation settings.
- `code/example.R`: applied METABRIC example.
- `code/scripts/fixed_lambda_sims.R`: fixed-lambda simulation runs.
- `code/scripts/selective_metrics_by_folder.R`: aggregate selective-inference metrics from saved runs.

The original paper code has been kept as close as possible to the archived version. Only machine-specific paths and repository documentation were adjusted; see `docs/PORTABILITY_CHANGES.md`.

## R packages

From the `code/` directory, install the main dependencies with:

```r
source("scripts/00_setup_packages.R")
```

The scripts also contain local checks for required packages. Package availability can change over time, so reproducing the exact historical software environment may require versions close to those used for the publication.

## Quick start: toy simulation

Run commands from the `code/` directory. For example:

```bash
Rscript run_analyze_toy.R --n=250 --p=10 --n_sim=10 --lambda_choice=fix
```

Most simulation settings can be overridden with `--key=value` command-line arguments. Defaults are documented directly in `run_analyze_toy.R` and `run_analyze_meta.R`.

## METABRIC data

The repository does **not** include the METABRIC source data. The applied example and METABRIC-informed simulation use the publicly available METABRIC clinical data from cBioPortal.

Place the downloaded clinical CSV at:

```text
code/data_raw/cbioportal/Breast Cancer METABRIC.csv
```

or set an environment variable pointing to the file before running R:

```bash
# Linux/macOS
export METABRIC_CSV="/path/to/Breast Cancer METABRIC.csv"

# Windows PowerShell
$env:METABRIC_CSV="C:\path\to\Breast Cancer METABRIC.csv"
```

The expected variables are those used in `code/helpers/funs.metabric_weibph.R` and `code/example.R`. The source study is the METABRIC cohort available through cBioPortal (`brca_metabric`).

Raw data and generated output are excluded from Git by default.

## Output

Simulation output is written as CSV/RDS/text files by the run scripts. Generated files should be stored locally under `code/results/` or another designated output directory and are not intended to be version-controlled.

## Reproducibility notes

- Random-number-generator settings are explicitly defined in the run scripts.
- Fixed-lambda values used by the simulation are stored in `code/config/lambda_values_used_in_plots_n100000.csv`.
- The code distinguishes method-specific inferential targets; interpretation of coverage and related metrics should follow the published article.
- Some package versions and external data interfaces may have changed since the original analyses in 2026.

## Citation

If you use this code, please cite the accompanying article. Citation metadata are also provided in `CITATION.cff`.

## License

No software license is assigned in this repository at present. Unless a license is added by the authors, reuse is governed by applicable copyright law and the terms of the published article/data sources.

## License

The source code in this repository is licensed under the MIT License.
See [LICENSE](LICENSE) for details.

This license applies to the source code and repository-specific configuration
files only. It does not grant rights to third-party datasets, external
software, or the published article.
