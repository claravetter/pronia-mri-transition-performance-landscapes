# Per-stage contract

| Stage | Required inputs | Produced outputs and boundary |
|---|---|---|
| Upstream scoring | Authorized acquisition/QC/harmonization; frozen classifier, BrainAGE, component maps/expression and historical probability objects | Imported fixed values. Raw MRI retraining, historical correction and site-transfer recovery are not performed here. Confirmed fields and unresolved upstream metadata are specified in the input contract. |
| Canonicalization | Explicit CSV path, accepted cohort/outcome coding, optional identity-specific correction JSON | One reconciled dataset; strict fixed-field/schema/count checks; private row identities and input reconciliation. ZInEP outcome is psychosis transition during its own follow-up. |
| Populations/transforms | Canonical dataset, required-variable sets | Per-analysis N/events and cohort composition, complete-case scaling, matched shared-variable transforms, training-only LOCO transforms. |
| Descriptive/validation | Frozen scores/classes/probabilities; cohorts and recorded sex | Baselines/missingness; primary and sex-stratified AUC/DeLong CI, BACC/confusion metrics/LR/Brier; within-cohort/per-cohort permutations; random-effects external cohort meta-AUC. Pooled-participant AUC is distinct. |
| Landscapes | Six moderators, six pairs; shared fixed-range kernel | Every interval/cell, midpoint, participant mean, support count and local uncertainty; Dmax nulls; six-test BH families; unrestricted/blocked and population/support/width sensitivities. |
| Moderation | Symptom and functioning complete populations | Nine GLMs including matched comparison; coefficients, within-model interaction BH, six nested comparisons, diagnostics/simple slopes; canonical model/scaling/row IDs and NAPLS-3 plot grid. |
| GAM | Same canonical symptom population and preserved formulas | Four bounded specifications, scaling/row-order comparisons, concurvity, indexed k diagnostics and Bernoulli simulations; pooled-support segmented grids. No supported ribbon crosses a hole. |
| Calibration | Original probabilities and explicit pooled/cohort/held-out populations | Original, apparent intercept, apparent slope, source and held-out assessments. Pooled apparent bootstrap refits each update. Fixed-update cohort/held-out assessments stay conditional. Updating alpha/beta are separate from assessment intercept/slope. |
| Conditional LOCO | Frozen MRI scores and complete clinical covariates | Training-only means/SDs and second-stage fits; paired held-out predictions, pooled AUC/Brier differences, cohort-stratified conditional bootstrap. No upstream MRI retraining. |
| Correlation/components | Fixed expression columns Comp1–4 (COMP01/02/04/06), common masks | Original Pearson/Spearman families; subject-level links; all 288 context/summary associations; twelve BH families of 24; explicit 48-cell selected-statistic key. |
| Registries | Every resampling family and model | Attempted/valid/invalid/extreme counts, exact seeds/stream states, uncertainty/correction definitions and method-review rows. Undefined statistics and invalid resampling draws remain explicit. |
| Plots | Only the new run exports and original fixed atlas/map assets | Main F1–4, S1–4 and S8; nine analytical/atlas figures total. Static reporting illustrations are maintained separately. |

## Schema and sensitive products

See `../data/README.md` for all required fields. Numeric 0/1 labels are explicit; sex follows 1=male, 2=female. Clinical missingness is preserved. Fixed MRI score/class/probability/component inputs must be valid. Cohort alias NAPLES normalizes to NAPLS-3. A duplicate may differ only in the clinical field and values specified by its external correction record; any other field disagreement or unknown identity is rejected.

Private products include participant IDs/canonical rows, model frames/objects and individual predictions. A private output folder is not a data-release decision. Aggregate exports, small cells and fixed atlas assets also require applicable disclosure/ownership review before publication. This repository contains no participant records, fitted models or analysis outputs.

## RNG and validity

A stable SHA256 key of base seed and analysis ID initializes L'Ecuyer-CMRG; actual nextRNGStream states index replicates. Hash-state collisions are checked. Worker assignment does not choose seeds. Serial, parallel and reversed scheduling were tested. R normal/sample generators are recorded as Inversion/Rejection.

A finite-draw conditional permutation p is (extreme+1)/(valid+1). Invalid attempts remain visible and trigger method review because conditioning can change a null distribution. No p is emitted for an undefined observation/no valid draws. Fisher uses fixed-margin random tables and the stats::fisher.test probability-order tie tolerance. mgcv k checks preserve the diagnostic proportion without relabelling it as an add-one primary test. Bootstrap intervals report finite draws per metric; nonconverged/boundary calibration fits are invalid.
