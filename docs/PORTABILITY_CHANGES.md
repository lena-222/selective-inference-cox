# Portability changes relative to the archived code bundle

The scientific analysis code is preserved as closely as possible. The following repository-level changes were made when preparing the GitLab version:

1. Removed the machine-specific absolute Windows path `C:/forschung/submission simstudy paper/code/data_raw/cbioportal/metabric_weibph_state.rds`.
2. The METABRIC-informed scripts now read the clinical CSV from either the `METABRIC_CSV` environment variable or the repository-relative path `code/data_raw/cbioportal/Breast Cancer METABRIC.csv`.
3. Added explicit checks with informative errors when the METABRIC CSV is unavailable.
4. Expanded `scripts/00_setup_packages.R` so that it lists packages used across the supplied scripts rather than only the applied example.
5. Corrected the repository-relative lookup of `lambda_values_used_in_plots_n100000.csv` to the existing `config/` directory.
6. The archived bundle contained `data_raw/cbioportal/metabric_weibph_state.rds`. This serialized object is intentionally omitted from the GitLab-ready repository because fitted R objects can embed row-level source data. The METABRIC-informed model is instead reconstructed from the separately downloaded public clinical CSV.
7. Added repository documentation, citation metadata, `.gitignore`, and placeholder directories for data/results.

No simulation estimands, model definitions, inference procedures, tuning rules, or numerical algorithms were intentionally changed.
