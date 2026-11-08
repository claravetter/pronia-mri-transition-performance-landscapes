# PRONIA performance-landscape analysis

Research software for PRONIA context/landscape analyses. One driver consumes one scientific configuration, computes from authorized cohort records and frozen upstream model outputs, then renders that run's statistical exports. No saved statistical results or synthetic dataset is used as a substitute for computation.

The repository contains source, tests, configuration and dependency locks. Participant data, fitted models, empirical outputs, manuscripts, anatomical exports and atlas spreadsheets are supplied separately and are not distributed here. Original MRI preprocessing/model training and historical BrainAGE correction are upstream dependencies.

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

The external directory must contain `source_panels/`, `yeo/` and `aal3/`. On an initial run without assets, the manifest records that anatomical figures were not requested. A later render with assets updates their hashes and provenance. A statistical-only rerender preserves previously generated anatomy only if its recorded output hashes still match, retaining its separate generating-source record in `anatomy_provenance`. Statistical analyses do not depend on those assets. Optional `--component-map-dir` records availability of the original NIfTI maps without changing fixed expression inputs.

Plot an already computed run without refitting:

```sh
.venv/bin/python python/run_pipeline.py render --output /authorized/results/pronia_NEW
```

Full settings include 5000 permutations, 2000 local AUC/pooled-calibration/LOCO bootstraps and separately defined 500-replicate calibration assessments. Synthetic settings are smaller and explicitly labelled. Production refuses existing output directories and paths inside the checkout, including symlink aliases. `--help` and imports do not read participant data or perform inference.

## External-validation analyses

The scripts in `R/analyses/` compute the remaining external-validation results of the paper. They reuse the helpers in `R/lib/`, so populations, windows and resampling streams match the main pipeline, and they write aggregate tables and figures only. Run them from the repository root with an output directory outside the checkout. If the input contains the documented duplicate key, set `PRONIA_CORRECTION_JSON` to the same record that `--correction` takes. `main_figures.R` also needs the R package `patchwork`.

```sh
export PRONIA_CORRECTION_JSON=/authorized/pronia_duplicate_resolution.json
INPUT=/authorized/all_cohorts_data.csv; OUT=/authorized/results/external_validation
Rscript R/analyses/external_auc_synthesis_and_calibration.R  $INPUT $OUT .
Rscript R/analyses/reference_shift_added_value_net_benefit.R $INPUT $OUT . 8
Rscript R/analyses/casemix_standardization.R                 $INPUT $OUT .
Rscript R/analyses/component_separation.R                    $INPUT $OUT . 8
Rscript R/analyses/scanner_site_check.R                      $INPUT $OUT . [SITE_MAP.csv]
Rscript R/analyses/landscape_power_simulation.R              $INPUT $OUT . 9 1000 199
Rscript R/analyses/negative_symptom_moderation_by_cohort.R   $INPUT $OUT .
Rscript R/analyses/napls3_subgroups_and_image_quality.R      $INPUT /authorized/napls3_covariates.csv $OUT .
Rscript R/analyses/pronia_profile_premorbid_onset_timing.R   $INPUT /authorized/pronia_pas_t0.csv /authorized/pronia_visits.csv $OUT .
Rscript R/analyses/main_figures.R                            $INPUT $OUT .
```

| Script | Content |
|---|---|
| `external_auc_synthesis_and_calibration.R` | Random-effects synthesis of external AUCs (Hartung-Knapp CI, prediction interval); cohorts used vs not used in the original report; calibration by cohort after intercept updating |
| `reference_shift_added_value_net_benefit.R` | PRONIA CHR-P reference; score distributions by cohort and outcome and what a per-cohort shift can explain; landscape robustness to population composition; symptom regression with and without the MRI score per held-out cohort; decision curves |
| `casemix_standardization.R` | Membership regressions and case-mix-standardized AUC of PRONIA CHR-P participants reweighted to each external cohort |
| `component_separation.R` | Separation of participants with and without transition on the MRI score and its 4 components, by cohort, with bootstrap tests of component-specific loss |
| `scanner_site_check.R` | Within-site AUC and site offsets of nontransition scores in multi-site cohorts |
| `landscape_power_simulation.R` | Operating characteristics of the sliding-window test at the study's sample sizes |
| `negative_symptom_moderation_by_cohort.R` | MRI-score AUC below vs above the cohort median SIPS negative score and the score × SIPS-N interaction per cohort |
| `napls3_subgroups_and_image_quality.R` | NAPLS-3 discrimination by CHR-P subtype, site, image quality and global gray matter volume |
| `pronia_profile_premorbid_onset_timing.R` | PRONIA CHR-P clinical profile by outcome, premorbid adjustment, and discrimination by whether the onset was documented at a study visit; needs the PRONIA participant number in `PatID` and the PRONIA PAS and visit exports |
| `main_figures.R` | Figures 1 and 2 and the supplementary figures built from the outputs above; runs last |

## Tests

```sh
Rscript tests/test_statistics.R
Rscript tests/test_production.R
.venv/bin/python -m unittest discover -s tests -p 'test_*.py' -v
```

Tests generate bounded examples in memory or temporary files. They verify production functions, not copied study results. The optional `tests/compare_results.py` compares a completed run with a separately supplied reference, reporting missing/schema/numeric differences. References are never production inputs.

Inspect `RUN_STATUS.txt`, `run_manifest.json`, `tables/resampling_validity_registry.csv` and `tables/method_review_required.csv` after execution. In empirical runs, `private/` contains protected records, row identities and fitted objects. Keep the entire output directory outside this code repository and apply cohort disclosure rules before sharing any output.

## Licence

MIT; see [LICENSE](LICENSE). The licence permits academic and commercial software reuse with the required notices. It grants no rights to separately held participant data, fitted models, atlas assets or manuscript materials. Dependencies retain their own licences.

## Citation and access

Cite Clara Vetter, *PRONIA performance-landscape analysis*, software version 0.3.1, with the exact commit used and this repository URL: https://github.com/claravetter/pronia-mri-transition-performance-landscapes. The companion [Prediction Landscape Explorer](https://github.com/claravetter/prediction-landscape-explorer) is a standalone descriptive app. A final paper/preprint identifier is not specified here; use the authors' confirmed study citation when available.

See [data requests and field definitions](data/README.md). Neither the software licence nor this documentation grants access to cohort records or upstream assets.
