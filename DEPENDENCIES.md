# Dependency restoration

`renv.lock` records exact R package versions; `scripts/restore_environment.R` restores packages into a project-local `.Rlibrary` and verifies exact versions, including recommended packages supplied by the R installation. Set `R_LIBS_USER` to that directory for subsequent commands. Use R 4.5.2 for the recorded environment.

`requirements.lock` pins the direct and transitive Python packages; `requirements.txt` lists the direct plotting/export requirements. Use Python 3.14.4 for the supported environment. Install into a separate virtual environment.

These locks do not contain result data, participant information, package binaries or fonts. Packages retain their original licences. External participant data, frozen model outputs and optional fixed anatomical assets are not software dependencies and are not bundled.

Parallel resampling uses Unix fork workers (tested on macOS). On Windows, use `--cores 1`. The operating system and numerical libraries can introduce small floating-point differences. Output comparisons use an explicit numerical tolerance (default 1e-10).

For strict repeat comparisons with OpenBLAS, set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` before starting the pipeline. `--cores` controls the independent R resampling workers. Other numerical libraries or thread settings may yield small differences in sensitive GAM diagnostics.
