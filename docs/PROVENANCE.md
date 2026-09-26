# Scientific and software provenance

This workflow analyses fixed PRONIA classifier predictions, harmonized cohort records, BrainAGE values and NeuroMiner component expressions. It evaluates the complete context/landscape comparison families and sensitivity analyses. Statistical outputs are computed from the supplied input for each run; no saved reference outputs are production inputs.

Clara Vetter maintains this implementation. Statistical methods and upstream study/model contributions retain their original authorship. Scientific interpretation remains the authors' responsibility. Dependencies retain their own credits and licences.

Original MRI preprocessing/training, historical BrainAGE correction and site transfer are upstream dependencies, not implemented stages. The [input contract](../data/README.md) distinguishes verified information from unresolved upstream metadata. [Anatomical assets](ANATOMY_INPUTS.md) must be supplied separately with permission.

## Original probability generation

The inspected [NeuroMiner implementation](https://github.com/neurominer-git/NeuroMiner/tree/29401f700905387f49191823ac2aaf1f62ebde4c) at commit `29401f700905387f49191823ac2aaf1f62ebde4c` uses class-weighted Platt scaling ([`util/probtransform.m`](https://github.com/neurominer-git/NeuroMiner/blob/29401f700905387f49191823ac2aaf1f62ebde4c/util/probtransform.m), called by `gui/NM_Results_Viewer.m`). A binomial-logit model is fitted to decision scores, with positive-class weights equal to negative/positive counts and negative-class weights of one. Source probabilities use leave-one-out fits; external probabilities use the full source fit. The export column is `Probs_platt`. The inspected branch applies no subsequent prevalence-prior correction.

This verifies the implementation, not the exact historical export settings: model objects, selected GUI options and optional score-offset settings must be linked to the input export by its provider. The downstream workflow assesses the supplied probabilities unchanged as original probabilities; its explicit apparent and held-out updates are separate analyses. It does not reconstruct historical probabilities or interpret class-balanced probabilities as automatically calibrated absolute risks.

Run manifests retain statistical input/configuration/software hashes. A later render adds its own software and asset provenance without replacing the original statistical provenance. `render_provenance` identifies the latest rendering operation; `anatomy_provenance` separately binds the anatomy-generating source, runtime and asset hashes to the retained anatomy output hashes. Statistics-only rendering preserves this anatomy record. If an older manifest cannot establish that source, its status is explicitly `NOT_RECORDED`; supplying the original authorized assets for a new anatomy render establishes a new generating record. Run directories, participant records, fitted objects, figures and atlas records are outside the code distribution.

Statistics-only rendering requires a nonempty `anatomy_output_hashes` record whose files still match. Missing or mismatched output hashes stop rendering; supply authorized assets to regenerate. `NOT_RECORDED` applies only when output hashes can be verified but the generating source cannot be established.
