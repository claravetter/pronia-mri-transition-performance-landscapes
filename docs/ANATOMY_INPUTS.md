# Optional external anatomical inputs

Pass `--anatomy-assets DIR` to request the two anatomical/atlas figure exports. No assets are shipped with this code.

- `DIR/source_panels/comp_1_source.png` through `comp_4_source.png`: original signed MRIcroGL exports for COMP01, COMP02, COMP04, COMP06, using the orientation, slices and signed thresholds specified below.
- `DIR/yeo/Component_compXX_CVRabsGT3_AND_SBCGT1p3/YeoJ_thresh13_noConCorr.mat` and `.xlsx`, for XX01,02,04,06.
- `DIR/aal3/Component_compXX_CVRabsGT3_AND_SBCGT1p3/AAL3_thresh13_noConCorr.mat` and `.xlsx`, for the same components.

The exporters verify the MATLAB/spreadsheet correspondence and voxel denominators. Yeo radar ordering follows the 27-region union and spreadsheet order. Zero-voxel AAL3 parcels remain undefined. Figure3 remains signed panels a–d plus radar e, with the component-performance heatmap separate in Supplementary Figure4.

`python/render_anatomy_mricrogl.py` optionally re-renders original signed NIfTI maps using MRIcroGL and the required template supplied separately. It is not raw MRI training or a substitute for the original model maps. Standard statistical analysis/figures work without anatomical assets; their omission is recorded explicitly.

## File schemas and interpretation

Each MATLAB file must load with `scipy.io.loadmat(..., simplify_cells=True)` and contain `vROI`, an array of entries with `ref` records. Required fields are `name` (region label), `nvoxROI` (atlas-region voxel count), `nvoxK` (overlapping component voxels), and `percROI` (100 × nvoxK / nvoxROI). Labels must be unique after trimming whitespace. Zero-voxel AAL3 parcels must have zero overlap and undefined occupancy; they are not imputed as zero occupancy.

Each corresponding XLSX workbook must have a first sheet with columns `AnatomicalRegion`, `K_ROI_P1[%]`, `ROIvox`, and `K_P1[vox]`. The exporter uses non-`NONE` regions with positive occupancy, matches their names to MATLAB records, and checks counts and percentages (absolute tolerance 1e-10). Workbook display ordering uses descending occupancy and the component-wise union. The Yeo display contains 27 regions, which are atlas region labels rather than 27 distinct canonical networks.

The study methods specify the Yeo 17-network parcellation and AAL3, absolute back-projected weights >3 and log-scaled sign stability >1.3. The four signed PNG panels use the spm152 template, radiological axial orientation, z=-50,-35,-20,-5,10,25,40,55,70 mm, with positive red/negative blue and limits ±15. The Figure 3 radar scale is 0–45% ROI occupancy. Asset providers must identify the exact atlas/template build, voxel grid and generating model/export revision; filenames alone do not establish those versions. Supplied files are hashed when consumed.

To add anatomy to a completed statistical run, invoke `python/run_pipeline.py render --output EXTERNAL_RUN --anatomy-assets AUTHORIZED_ASSETS`. The manifest records input asset hashes, generated output hashes and a separate `anatomy_provenance` object containing the generating source/runtime and those same hashes. Later statistics-only renders verify the retained files and preserve this object. The original statistical source/input hashes also remain unchanged. For older manifests the anatomy generator is recovered only if the previous render record proves it generated the retained figures; otherwise it is explicitly `NOT_RECORDED`, never attributed to the latest statistics-only render. Regenerate with authorized assets to record a new generator. If rendering fails, `render_status` is `FAILED`; partial products are not a completed figure set.

Statistics-only rendering requires a nonempty `anatomy_output_hashes` record whose files still match. Missing or mismatched output hashes stop rendering; supply authorized assets to regenerate. `NOT_RECORDED` applies only when output hashes can be verified but the generating source cannot be established.

## MRIcroGL launch

Set `PRONIA_COMPONENT_MAP_DIR`, `PRONIA_PANEL_OUTPUT_DIR` (a new directory) and `PRONIA_MRICROGL_TEMPLATE`, then launch a fresh MRIcroGL process:

```sh
.venv/bin/python python/render_anatomy_mricrogl.py \
  --mricrogl /path/to/MRIcroGL --reset-preferences
```

The launcher passes `-r` before `-s`, selecting MRIcroGL's radiological default (left on right) at startup. This resets MRIcroGL's saved GUI preferences; use a dedicated rendering installation/account if those settings must be preserved. Calling `gl.resetdefaults()` inside an existing session does not set the orientation. The embedded script therefore requires the launcher. This follows MRIcroGL's [startup implementation](https://github.com/rordenlab/MRIcroGL/blob/master/mainunit.pas) and [radiological default](https://github.com/rordenlab/MRIcroGL/blob/master/prefs.pas).

The launcher verifies four nonempty panel files and the completion marker. It closes the process it started after rendering if the GUI remains open, and stops incomplete renders after `--timeout` seconds (default 300).
