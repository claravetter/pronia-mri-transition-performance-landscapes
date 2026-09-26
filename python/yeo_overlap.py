#!/usr/bin/env python3
"""Export recorded vROI occupancy values in the associated spreadsheet region order."""
import argparse
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    tables, full, manifest = [], [], []
    for i, tag in enumerate(['comp01','comp02','comp04','comp06'], 1):
        folder = a.source_dir / f'Component_{tag}_CVRabsGT3_AND_SBCGT1p3'
        xp, mp = folder/'YeoJ_thresh13_noConCorr.xlsx', folder/'YeoJ_thresh13_noConCorr.mat'
        for src in [xp,mp]:
            if not src.is_file(): raise FileNotFoundError(f'Required Yeo source file missing: {src}')
            manifest.append(dict(source_file=str(src.relative_to(a.source_dir)),sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
        t = pd.read_excel(xp)
        t = t[(t['AnatomicalRegion'] != 'NONE') & (t['K_ROI_P1[%]'] > 0)].sort_values('K_ROI_P1[%]',ascending=False)
        tables.append(t[['AnatomicalRegion','K_ROI_P1[%]']].rename(columns={'K_ROI_P1[%]':f'Comp {i}'}))
        records = [z['ref'] for z in loadmat(mp,simplify_cells=True)['vROI']]
        byname = {r['name'].strip():r for r in records}
        assert len(byname) == len(records)
        for _, r in t.iterrows():
            q = byname[r.AnatomicalRegion.strip()]
            assert np.isclose(q['percROI'],r['K_ROI_P1[%]'],rtol=0,atol=1e-10)
            assert q['nvoxK'] == r['K_P1[vox]'] and q['nvoxROI'] == r.ROIvox
        for r in records:
            assert np.isclose(r['percROI'],100*r['nvoxK']/r['nvoxROI'],rtol=0,atol=1e-10)
            full.append(dict(manuscript_label=f'Comp {i}',saved_identifier=tag.upper(),yeo_network=r['name'].strip(),
                network_roi_voxels=r['nvoxROI'],component_overlap_voxels=r['nvoxK'],network_occupancy_percent=r['percROI'],
                reported_in_thresholded_xlsx=r['name'].strip() in {s.strip() for s in t.AnatomicalRegion},
                source_field='vROI.ref.percROI',source_file=str(mp.relative_to(a.source_dir))))
    # Original outer merge and dominant-component ordering. NaN entries need no zero-fill to order maxima.
    merged=tables[0]
    for t in tables[1:]: merged=merged.merge(t,on='AnatomicalRegion',how='outer')
    m=merged.set_index('AnatomicalRegion')
    cols=list(m.columns); dominant=m.idxmax(axis=1); maximum=m.max(axis=1)
    order=sorted(m.index,key=lambda r:(cols.index(dominant[r]),-float(maximum[r])))
    frame=pd.DataFrame(full); order=[r.strip() for r in order]
    display=frame[frame.yeo_network.isin(order)].copy()
    display['network_order']=display.yeo_network.map({r:i+1 for i,r in enumerate(order)})
    display=display.sort_values(['network_order','manuscript_label'])
    assert len(display)==4*len(order) and not display.network_occupancy_percent.isna().any()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    frame.to_csv(a.output_dir/'yeo_recorded_full_overlap.csv',index=False)
    display.to_csv(a.output_dir/'yeo_radar_display.csv',index=False)
    pd.DataFrame(manifest).to_csv(a.output_dir/'yeo_radar_source_manifest.csv',index=False)
    print(f'PASS: {len(frame)} recorded values; {len(order)} original displayed networks; {len(display)} radar values; no imputation')

if __name__=='__main__': main()
