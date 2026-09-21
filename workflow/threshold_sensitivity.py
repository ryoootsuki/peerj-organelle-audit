#!/usr/bin/env python3
from __future__ import annotations
import argparse, itertools
import pandas as pd
from common import load_config, results_dir


def cat_one(x, low, high):
    return 'Retain' if x < low else ('Sensitivity analysis' if x < high else 'Exclude')


def cat_framework(strict, multi, low, high, ambiguity):
    if strict >= high:
        return 'Exclude'
    if strict < low and multi < ambiguity:
        return 'Retain'
    return 'Sensitivity analysis'


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); a=ap.parse_args()
    cfg=load_config(a.config); tdir=results_dir(cfg)/'tables'
    df=pd.read_csv(tdir/'audit_primary_mapq.tsv',sep='\t')
    for c in ['frac_strict_total','frac_multi_total','frac_inclusive_total']:
        df[c]=pd.to_numeric(df[c],errors='coerce')
    grid=cfg.get('triage_sensitivity',{})
    rows=[]
    for low,high,amb in itertools.product(grid.get('low_fraction_grid',[0.05]), grid.get('high_fraction_grid',[0.20]), grid.get('ambiguity_fraction_grid',[0.01])):
        if float(low)>=float(high):
            continue
        one=df.frac_inclusive_total.map(lambda x:cat_one(x,float(low),float(high)))
        fw=df.apply(lambda x:cat_framework(x.frac_strict_total,x.frac_multi_total,float(low),float(high),float(amb)),axis=1)
        rows.append({'low_fraction':low,'high_fraction':high,'ambiguity_fraction':amb,'n_runs':len(df),
                     'retain':int((fw=='Retain').sum()),'sensitivity_analysis':int((fw=='Sensitivity analysis').sum()),'exclude':int((fw=='Exclude').sum()),
                     'reclassified_vs_one_threshold':int((fw!=one).sum()),'reclassified_fraction':float((fw!=one).mean())})
    pd.DataFrame(rows).to_csv(tdir/'triage_threshold_sensitivity.tsv',sep='\t',index=False)
    print(f"[OK] threshold sensitivity written for {cfg['species']['id']}")

if __name__=='__main__':
    main()
