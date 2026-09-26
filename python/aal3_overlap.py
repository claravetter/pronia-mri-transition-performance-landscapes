#!/usr/bin/env python3
"""Compute AAL3 regional occupancy and render Supplementary Figure 3 from atlas records."""
import argparse
import hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    a = parser.parse_args()
    full, selected, provenance = [], [], []
    fig, axes = plt.subplots(2, 2, figsize=(12, 14), layout='constrained')
    for i, (tag, ax, color) in enumerate(zip(['comp01','comp02','comp04','comp06'], axes.flat,
                                           ['#0f4c5c','#e36414','#6a994e','#7b2cbf']), 1):
        folder = a.source_dir / f'Component_{tag}_CVRabsGT3_AND_SBCGT1p3'
        xp, mp = folder/'AAL3_thresh13_noConCorr.xlsx', folder/'AAL3_thresh13_noConCorr.mat'
        for p in (xp, mp):
            provenance.append(dict(source=str(p.relative_to(a.source_dir)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        frame = pd.read_excel(xp)
        frame = frame.loc[(frame.AnatomicalRegion != 'NONE') & (frame['K_ROI_P1[%]'] > 0)].copy()
        records = [q['ref'] for q in loadmat(mp, simplify_cells=True)['vROI']]
        byname = {q['name'].strip(): q for q in records}
        assert len(byname) == len(records)
        shown = set(frame.AnatomicalRegion.str.strip())
        for _, row in frame.iterrows():
            q = byname[row.AnatomicalRegion.strip()]
            assert np.isclose(row['K_ROI_P1[%]'], q['percROI'], rtol=0, atol=1e-10)
            assert row.ROIvox == q['nvoxROI'] and row['K_P1[vox]'] == q['nvoxK']
        for q in records:
            if q['nvoxROI'] == 0:
                assert q['nvoxK'] == 0 and np.isnan(q['percROI'])  # absent atlas parcels remain undefined
            else:
                assert np.isclose(q['percROI'], 100*q['nvoxK']/q['nvoxROI'], rtol=0, atol=1e-10)
            full.append(dict(component=f'Comp {i}', saved_identifier=tag.upper(), region=q['name'].strip(),
                             ROI_voxels=q['nvoxROI'], overlap_voxels=q['nvoxK'], occupancy_percent=q['percROI'],
                             displayed=q['name'].strip() in shown, source_field='vROI.ref.percROI'))
        frame = frame.sort_values('K_ROI_P1[%]', ascending=False)
        selected.append(frame.assign(component=f'Comp {i}', saved_identifier=tag.upper(), display_order=np.arange(len(frame))+1))
        ax.barh(frame.AnatomicalRegion, frame['K_ROI_P1[%]'], color=color)
        ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel('Parcel occupancy (%)')
        ax.set_title(f'{chr(96+i)}  Comp {i}', loc='left', fontweight='bold')
        ax.tick_params(axis='y', labelsize=8)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('AAL3 summaries of fixed suprathreshold absolute component maps', fontsize=13)
    data = a.output_dir/'figure_data/aal3'; data.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(full).to_csv(data/'aal3_recorded_full_overlap.csv', index=False)
    pd.concat(selected).to_csv(data/'aal3_display.csv', index=False)
    pd.DataFrame(provenance).to_csv(data/'aal3_source_manifest.csv', index=False)
    figs = a.output_dir/'figures'; figs.mkdir(exist_ok=True)
    for ext in ['pdf','png']:
        fig.savefig(figs/f'supplementary_figure3_aal3.{ext}', dpi=200)
    plt.close(fig)
    print(f'PASS: AAL3 {len(full)} recorded values; {sum(len(x) for x in selected)} displayed parcels; original voxel denominators')

if __name__ == '__main__':
    main()
