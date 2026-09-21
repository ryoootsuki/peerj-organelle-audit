#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import yaml


def check_file(path: Path):
    return path.exists() and path.stat().st_size>0

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--strict',action='store_true'); args=ap.parse_args()
    checks=[]
    for sid,expected in [('silene',2679),('chicken',141)]:
        comp=Path('results')/sid/'provenance'/'analysis_completion_summary.tsv'
        if check_file(comp):
            d=pd.read_csv(comp,sep='\t'); obs=int(d.iloc[0]['observed_primary_mapq_runs']); checks.append((f'{sid}_full_audit_count',obs==expected,f'observed={obs}, expected={expected}'))
        else: checks.append((f'{sid}_full_audit_count',False,'completion summary missing'))
        for label,path in [
            ('reference_provenance',Path('work')/sid/'references'/'reports'/'reference_provenance.tsv'),
            ('mask_summary',Path('work')/sid/'references'/'reports'/'mask.main.summary.tsv'),
            ('protocol_metadata',Path('results')/sid/'tables'/'library_protocol_metadata.tsv'),
            ('downstream_variant_deltas',Path('results')/sid/'tables'/'downstream_variant_deltas.tsv'),
            ('reference_mask_sensitivity',Path('results')/sid/'tables'/'reference_mask_sensitivity.tsv'),
            ('reference_mask_variant_deltas',Path('results')/sid/'tables'/'reference_mask_variant_deltas.tsv'),
        ]:
            checks.append((f'{sid}_{label}',check_file(path),str(path)))
    checks.append(('workflow_figure',check_file(Path('results/combined/figures/Figure_1_workflow_overview.png')),'results/combined/figures/Figure_1_workflow_overview.png'))
    checks.append(('combined_species_summary',check_file(Path('results/combined/species_validation_summary.tsv')),'results/combined/species_validation_summary.tsv'))
    out=Path('results/combined/report_materials/final_output_validation.tsv'); out.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame([{'check':n,'status':'PASS' if ok else 'FAIL','detail':detail} for n,ok,detail in checks]).to_csv(out,sep='\t',index=False)
    failed=[n for n,ok,_ in checks if not ok]
    print(f'[INFO] output validation: {len(checks)-len(failed)}/{len(checks)} checks passed')
    if failed: print('[WARN] failed checks: '+', '.join(failed))
    if args.strict and failed: raise SystemExit(1)
if __name__=='__main__': main()
