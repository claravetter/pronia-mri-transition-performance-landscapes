# PRONIA performance-landscape analysis

Research software for PRONIA context/landscape analyses. One driver consumes one scientific configuration, computes from authorized cohort records and frozen upstream model outputs, then renders that run's statistical exports. No saved statistical results or synthetic dataset is used as a substitute for computation.

The repository contains source, tests, configuration, dependency locks and documentation. Participant data, fitted models, empirical outputs, manuscripts, anatomical exports and atlas spreadsheets are supplied separately and are not distributed here. Original MRI preprocessing/model training and historical BrainAGE correction are upstream dependencies.

## Install

Supported environment: R 4.5.2 and Python 3.14.4. R and any required system compilers/libraries must be installed first. Run inside this repository:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock
Rscript scripts/restore_environment.R
export R_LIBS_USER="$PWD/.Rlibrary"
```

`renv.lock` specifies exact R package versions. The restore script installs packages into `.Rlibrary` and verifies exact versions, including recommended packages supplied by R itself. It does not read participant data or run analyses. `requirements.lock` pins the Python packages, including transitive dependencies. These are software environment locks, not statistical result snapshots.

## Run

A complete small-resample check generates its own synthetic input and writes to a new external directory:

```sh
.venv/bin/python python/run_pipeline.py synthetic --output ../local_runs/synthetic_run --cores 2
Rscript tests/verify_computed_run.R ../local_runs/synthetic_run
```

Validate input fields, cohort counts and PRONIA clinical-group labels without inference. Complete-case population counts are checked during a full run:

```sh
.venv/bin/python python/run_pipeline.py check \
  --input /authorized/all_cohorts_data.csv \
  --correction /authorized/pronia_duplicate_resolution.json
```

For authorized local empirical analysis:

```sh
.venv/bin/python python/run_pipeline.py run \
  --input /authorized/all_cohorts_data.csv \
  --correction /authorized/pronia_duplicate_resolution.json \
  --output /authorized/results/pronia_NEW \
  --config config/study.yml --cores 4 --authorize-private-analysis
```

See [input contract](data/README.md) for required columns and duplicate identity verification. Omit `--correction` if the file is already uniquely reconciled. Standardization is computed from the supplied unstandardized variables within each complete-case or training population.

Statistical figures are generated automatically. Anatomical Figure3 and Supplementary Figure3 additionally require the original, authorized fixed assets:

```sh
.venv/bin/python python/run_pipeline.py render \
  --output /authorized/results/pronia_NEW \
  --anatomy-assets /authorized/anatomy_assets
```

The external directory must contain `source_panels/`, `yeo/` and `aal3/` as described in [anatomy inputs](docs/ANATOMY_INPUTS.md). On an initial run without assets, the manifest records that anatomical figures were not requested. A later render with assets updates their hashes and provenance. A statistical-only rerender preserves previously generated anatomy only if its recorded output hashes still match, retaining its separate generating-source record in `anatomy_provenance`. Statistical analyses do not depend on those assets. Optional `--component-map-dir` records availability of the original NIfTI maps without changing fixed expression inputs.

Plot an already computed run without refitting:

```sh
.venv/bin/python python/run_pipeline.py render --output /authorized/results/pronia_NEW
```

Full settings include 5000 permutations, 2000 local AUC/pooled-calibration/LOCO bootstraps and separately defined 500-replicate calibration assessments. Synthetic settings are smaller and explicitly labelled. Production refuses existing output directories and paths inside the checkout, including symlink aliases. `--help` and imports do not read participant data or perform inference.

## Tests

```sh
Rscript tests/test_statistics.R
Rscript tests/test_production.R
.venv/bin/python -m unittest discover -s tests -p 'test_*.py' -v
```

Tests generate bounded examples in memory or temporary files. They verify production functions, not copied study results. The optional `tests/compare_results.py` compares a completed run with a separately supplied reference, reporting missing/schema/numeric differences. References are never production inputs.

Inspect `RUN_STATUS.txt`, `run_manifest.json`, `tables/resampling_validity_registry.csv` and `tables/method_review_required.csv` after execution. See [method boundaries](docs/METHOD_CONTRACT.md), [output map](docs/SOURCE_TO_MANUSCRIPT.md) and [tests](VALIDATION.md). In empirical runs, `private/` contains protected records, row identities and fitted objects. Keep the entire output directory outside this code repository and apply cohort disclosure rules before sharing any output.

## Licence

MIT; see [LICENSE](LICENSE). The licence permits academic and commercial software reuse with the required notices. It grants no rights to separately held participant data, fitted models, atlas assets or manuscript materials. Dependencies retain their own licences.

## Citation and access

Cite Clara Vetter, *PRONIA performance-landscape analysis*, software version 0.3.1, with the exact commit used and this repository URL: https://github.com/claravetter/pronia-mri-transition-performance-landscapes. The companion [Prediction Landscape Explorer](https://github.com/claravetter/prediction-landscape-explorer) is a standalone descriptive app. A final paper/preprint identifier is not specified here; use the authors' confirmed study citation when available.

See [data requests and field definitions](data/README.md). Neither the software licence nor this documentation grants access to cohort records or upstream assets.
