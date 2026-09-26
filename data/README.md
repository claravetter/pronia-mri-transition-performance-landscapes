# Authorized participant input

Pass an authorized CSV explicitly with `--input`. Do not commit it.

Required fields:
`Cases`, `cohort`, `Studygroup`, `EXP_LABEL`, `PRED_LABEL`, `Mean_Score`, `Probs_platt`, `age`, `sex`, `COGDIS_score`, `SIPS_Positiv_Gesamt`, `SIPS_Negativ_Gesamt`, `Psychosoz_aequiv`, `BrainAGE_corr`, `PredictedAge_SBC_corr`, `BrainAGE_SBC_corr`, and `Comp1`–`Comp4`.

`Studygroup` must be `CHR` or `ROP` for every PRONIA record, with at least one CHR record present. Unknown or missing PRONIA labels are rejected; external cohorts may have missing group labels. The CHR sensitivity includes PRONIA CHR records and all external records. Sex codes are 1=male and 2=female. Missing clinical covariates remain missing and are handled by the existing analysis-specific complete-case rules.

Outcome 1 denotes cohort-defined transition in every cohort, including ZInEP. The input contains supplied MRI predictions and probabilities, BrainAGE values and fixed NeuroMiner component scores. Their upstream generation is not part of this pipeline.

Accept either the resolved unique participant file or an export with one authorized two-row clinical-field conflict. When duplicate cohort/participant keys are present, supply an external JSON record through `--correction` with exactly these fields:

| Field | Required value |
|---|---|
| `key_sha256` | Lowercase SHA-256 of `Cases` alone as a UTF-8 string, without a trailing newline; the cohort is matched separately |
| `cohort` | Cohort label after alias normalization |
| `field` | The clinical field whose two recorded values conflict |
| `conflicting_values` | Array of exactly two distinct, finite numeric values expected in that field |
| `retained_value` | One of those values, identifying the row to retain |

Correctable fields are `age`, `sex`, `COGDIS_score`, `SIPS_Positiv_Gesamt`, `SIPS_Negativ_Gesamt` and `Psychosoz_aequiv`. Every other field must agree across the two rows, including outcome, fixed predictions, BrainAGE and component scores. The importer rejects an unknown identity, mismatched cohort/values, additional duplicates or other field disagreements. It never averages values or selects by row position. Re-importing the unique resolved file is idempotent and does not require a correction record. Keep the completed correction JSON with the controlled inputs outside the repository.

Required canonical counts: PRONIA 329/26; NAPLS-3 458/51; ZInEP 157/23; FePsy 37/16; MUC-FRUE 33/16 (participants/transitions).

Participant identifiers are confined to the chosen output's `private/` directory. Fitted objects and individual predictions must remain controlled. Tables and figures contain aggregate results; disclosure approval is still required before publication.

Native component maps are optional for recording local map availability (`--component-map-dir`) or re-rendering anatomy with MRIcroGL. They are not required for analysis of the supplied component scores.

Imported derived z-score columns are removed before reconciliation; the workflow recomputes scaling within each analysis/training population. Only the explicit derived-column names in `canonicalize_data` are discarded. Every other raw/fixed field must agree except the field named in the authorized correction.

## Field dictionary

The importer accepts CSV missing values `NA` and empty fields. Do not replace missing clinical values with zero. Do not rename cohort endpoints to impose a common follow-up horizon. One record represents a participant within a cohort; `Cases` is combined with cohort to check identity.

| Fields | Meaning and encoding | Unit / provenance boundary |
|---|---|---|
| `Cases` | Nonempty cohort-local participant identifier | String; confidential; not an analysis predictor |
| `cohort` | PRONIA, NAPLS-3, ZInEP, FePsy or MUC-FRUE; `NAPLES` is normalized to NAPLS-3 | Study membership |
| `Studygroup` | Supplied clinical group; the CHR sensitivity retains PRONIA rows only when `Studygroup == "CHR"`, together with all external rows | PRONIA: `CHR` or `ROP`; provider must document harmonization |
| `EXP_LABEL` | 0=no recorded transition, 1=transition during the cohort-specific observation period | Binary endpoint; not a common-horizon survival outcome |
| `PRED_LABEL` | Inherited fixed classifier decision, 0/1 | Do not derive a replacement from score sign or a new probability threshold |
| `Mean_Score` | Frozen continuous MRI classifier score; higher values indicate greater transition risk | Model-specific score units, not probability |
| `Probs_platt` | Supplied original Platt-scaled probability in [0,1] | See [probability provenance](../docs/PROVENANCE.md); historical export settings require the provider's records |
| `age` | Chronological age at the supplied assessment | Years |
| `sex` | Recorded 1=male, 2=female | Categorical; missingness preserved |
| `COGDIS_score` | Supplied cognitive-disturbance summary | Instrument score; exact item aggregation/range must be confirmed by the input provider |
| `SIPS_Positiv_Gesamt` | Supplied SIPS positive-symptom total | Instrument score; preserve harmonized source coding |
| `SIPS_Negativ_Gesamt` | Supplied SIPS negative-symptom total | Instrument score; preserve harmonized source coding |
| `Psychosoz_aequiv` | Harmonized psychosocial-functioning measure | Exact cross-instrument conversion and scale anchors require provider confirmation |
| `BrainAGE_corr` | Supplied alternative corrected brain-age gap, used for identity/diagnostic comparisons | Years; correction is upstream |
| `PredictedAge_SBC_corr` | Supplied corrected MRI-predicted age | Years; correction/site transfer is upstream |
| `BrainAGE_SBC_corr` | Primary corrected predicted-age minus chronological-age gap | Years; identity with predicted age minus age is checked |
| `Comp1`–`Comp4` | Fixed subject-level component expression | Model-specific dot-product units; respectively COMP01, COMP02, COMP04, COMP06 |

Fixed outcome/class/score/probability/component values are validated before inference; required covariates determine each analysis's complete-case population. Imported derived z-score columns are not authoritative input fields. The [method contract](../docs/METHOD_CONTRACT.md) specifies complete-case and held-out scaling.

## Upstream models and cohort endpoints

The supplied study methods describe a PRONIA-derived linear SVM and frozen predictions evaluated without refitting the MRI model. The BrainAGE model used NeuroMiner 1.4, 1,105 healthy controls aged 12–65, five outer and five inner folds, and corrected MRI-predicted age minus chronological age. These are study-method descriptions, not verification of the exact historical model binary. Providers must supply model/export identifiers and document historical corrections/site transfer for upstream reproduction.

The supplied Supplementary Table 10 defines psychosis transition under each cohort's protocol. PRONIA, NAPLS-3, ZInEP, FePsy and MUC-FRUE retain their original endpoints and observation periods; the workflow does not infer follow-up from the binary label or estimate common-time risk. Exact cohort eligibility, follow-up and treatment-context metadata should accompany an authorized input agreement or the final study supplement, not be inferred from column names.

## Data requests

Requests for FePsy data should be directed to the Center for Gender Research and Early Detection, University of Basel Psychiatric Hospital, Basel, Switzerland.

For MUC-FRUE, direct requests to Nikolaos Koutsouleris. For the other cohorts and upstream model/anatomical assets, contact the relevant study custodians through the study authors. No response timeframe, automatic access entitlement, onward-sharing permission or approved derived-data scope is asserted here; those conditions must be agreed with the respective custodians.
