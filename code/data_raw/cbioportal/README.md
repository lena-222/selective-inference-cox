# METABRIC data

The METABRIC clinical source data are not distributed with this repository.

Download the publicly available METABRIC clinical data from cBioPortal (`brca_metabric`) and either save the clinical CSV as `Breast Cancer METABRIC.csv` in this directory or set the `METABRIC_CSV` environment variable to the downloaded file.

The scripts expect the clinical variables referenced in `../../helpers/funs.metabric_weibph.R` and `../../example.R`.

Do not commit local data files to the repository.
