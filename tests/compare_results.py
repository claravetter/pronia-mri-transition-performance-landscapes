#!/usr/bin/env python3
"""Compare a NEW complete analysis to an explicit external reference; never generate results."""
import argparse
import csv
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--actual',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--tolerance',type=float,default=1e-10)
    a=p.parse_args()
    contract=Path(__file__).resolve().parents[1]/'config/output_contract.csv'
    with contract.open(newline='') as f:
        missing=[r['path'] for r in csv.DictReader(f) if not (a.actual/r['path']).is_file()]
    if missing:
        raise SystemExit('Actual run is incomplete: '+', '.join(missing))
    status=a.actual/'RUN_STATUS.txt'
    if not status.is_file() or 'status: COMPUTED' not in status.read_text():
        raise SystemExit('Actual run lacks a COMPUTED marker; refusing an equivalence claim.')
    if 'WITH_METHOD_REVIEW' in status.read_text():
        raise SystemExit('Resolve recorded method issues before requesting an equivalence claim.')
    rows=[]
    metadata={'runtime_rng_manifest.csv','component_runtime_rng_manifest.csv','component_identifier_mapping.csv',
              'component_import_matching_audit.csv','correction_registry.csv'}
    failures=[]
    for folder in ('tables','figure_data'):
        for path in sorted((a.actual/folder).glob('*.csv')):
            if path.name in metadata:
                rows.append(dict(file=str(path.relative_to(a.actual)),status='METADATA',detail='Run-specific or consolidated registry; inspect separately'))
                continue
            candidates=[a.reference/folder/path.name]
            matches=[x for x in candidates if x.is_file()]
            if len(matches)!=1:
                rows.append(dict(file=str(path.relative_to(a.actual)),status='NO_UNIQUE_REFERENCE',detail='Not numerically verified'))
                failures.append(path.name)
                continue
            x,y=pd.read_csv(path),pd.read_csv(matches[0])
            status='PASS';detail='All columns and rows match within declared tolerance'
            if len(x)!=len(y) or set(x.columns)!=set(y.columns):
                status='FAILED';detail='Row count or column set differs'
            else:
                for col in x:
                    if pd.api.types.is_numeric_dtype(x[col]) and pd.api.types.is_numeric_dtype(y[col]):
                        equal=np.allclose(x[col],y[col],rtol=a.tolerance,atol=a.tolerance,equal_nan=True)
                    else:
                        equal=x[col].fillna('<NA>').astype(str).equals(y[col].fillna('<NA>').astype(str))
                    if not equal:
                        status='FAILED';detail='First differing column: '+col;break
            rows.append(dict(file=str(path.relative_to(a.actual)),status=status,detail=detail))
            if status!='PASS':failures.append(path.name)
    if not rows:raise ValueError('No analysis tables or figure data found')
    a.report.parent.mkdir(parents=True,exist_ok=True)
    with a.report.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['file','status','detail']);w.writeheader();w.writerows(rows)
    if failures:raise SystemExit('Comparison incomplete or failed; inspect report. A full numerical-equivalence claim is not supported.')
    print('PASS: available comparable output files agree. Metadata rows require separate inspection.')

if __name__=='__main__':main()
