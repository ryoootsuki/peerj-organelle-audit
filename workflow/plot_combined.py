#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    out=Path('results/combined/figures'); out.mkdir(parents=True,exist_ok=True)
    frames=[]
    for sid in ['silene','chicken']:
        p=Path('results')/sid/'tables'/'audit_primary_mapq.tsv'
        if p.exists():
            d=pd.read_csv(p,sep='\t'); d['species']=sid; frames.append(d)
    if len(frames)<2:
        print('[WARN] combined species figure skipped because both species are not complete'); return
    d=pd.concat(frames,ignore_index=True)
    fig,ax=plt.subplots(figsize=(7,4.8)); order=['silene','chicken']; labels=[r'$\it{Silene\ latifolia}$',r'$\it{Gallus\ gallus}$']; ax.boxplot([100*d.loc[d.species==s,'frac_strict_total'].dropna() for s in order],tick_labels=labels,showfliers=False); ax.set_ylabel('Strict organelle-compatible fraction (%)'); ax.set_title('Cross-taxon transferability check at primary MAPQ'); fig.tight_layout(); fig.savefig(out/'Figure_10_silene_chicken_comparison.png',dpi=300,bbox_inches='tight'); fig.savefig(out/'Figure_10_silene_chicken_comparison.pdf',bbox_inches='tight'); plt.close(fig)
    summary=d.groupby('species').agg(n_runs=('Run','nunique'),median_strict=('frac_strict_total','median'),p95_strict=('frac_strict_total',lambda x:x.quantile(.95)),median_multi=('frac_multi_total','median'),max_inclusive=('frac_inclusive_total','max')).reset_index(); summary.to_csv('results/combined/species_validation_summary.tsv',sep='\t',index=False)
    print('[OK] combined Silene/chicken outputs written')
if __name__=='__main__': main()
