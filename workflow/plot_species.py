#!/usr/bin/env python3
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common import load_config, results_dir

ORDER=['Retain','Sensitivity analysis','Exclude']

def sci_math(name: str) -> str:
    return r'$\it{' + r'\ '.join(str(name).split()) + '}$'


def save(fig,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.tight_layout()
    fig.savefig(path,dpi=300,bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight')
    plt.close(fig)


def horizontal_two_condition(d, value, ylabel, title, path, smaller_is_worse=False):
    if d.empty: return
    agg=d.groupby('Run')[value].agg(lambda x: x.min() if smaller_is_worse else x.abs().max())
    runs=agg.sort_values(ascending=smaller_is_worse).index.tolist()
    y=np.arange(len(runs),dtype=float)
    fig_h=max(5.0,0.28*len(runs)+1.8)
    fig,ax=plt.subplots(figsize=(9.5,fig_h))
    conditions=list(dict.fromkeys(d['condition'].astype(str)))
    offsets=np.linspace(-0.16,0.16,max(len(conditions),1))
    markers=['o','x','s','^']
    for i,cond in enumerate(conditions):
        sub=d.loc[d.condition.astype(str).eq(cond)].set_index('Run').reindex(runs)
        ax.scatter(sub[value],y+offsets[i],label=cond,marker=markers[i%len(markers)],s=28)
    if not smaller_is_worse: ax.axvline(0,linewidth=.8)
    ax.set_yticks(y,runs,fontsize=7)
    ax.set_xlabel(ylabel)
    ax.set_ylabel('Run')
    ax.set_title(title)
    ax.legend(fontsize=8)
    save(fig,path)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); args=ap.parse_args()
    cfg=load_config(args.config); sid=cfg['species']['id']; sci=cfg['species']['scientific_name']; sci_t=sci_math(sci)
    t=results_dir(cfg)/'tables'; f=results_dir(cfg)/'figures'; f.mkdir(parents=True,exist_ok=True)
    audit=pd.read_csv(t/'audit_mapq_sweep.tsv',sep='\t'); primary=pd.read_csv(t/'one_threshold_vs_framework.tsv',sep='\t')
    nuclear_available=bool(cfg['references'].get('nuclear_assembly_accession')); qs=sorted(audit.MAPQ.unique())
    metric='frac_strict_total' if nuclear_available else 'frac_inclusive_total'
    label='Strict organelle-compatible fraction (%)' if nuclear_available else 'Inclusive organelle-compatible fraction (%)'
    stem='strict' if nuclear_available else 'inclusive'

    # MAPQ summary figure: medians + IQR + p95; full numerical values are also written to mapq_distribution_summary.tsv.
    summ=[]
    for q in qs:
        s=100*audit.loc[audit.MAPQ==q,metric].dropna()
        summ.append((q,s.median(),s.quantile(.25),s.quantile(.75),s.quantile(.95)))
    a=np.array(summ,float)
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    ax.errorbar(a[:,0],a[:,1],yerr=np.vstack([a[:,1]-a[:,2],a[:,3]-a[:,1]]),fmt='o-',capsize=3,label='median and IQR')
    ax.plot(a[:,0],a[:,4],'x--',label='95th percentile')
    ax.set(xlabel='Minimum MAPQ',ylabel=label,title=f'{sci_t}: organelle-compatible signal across MAPQ')
    ax.legend(); save(fig,f/f'Figure_2_{stem}_by_MAPQ.png')

    if nuclear_available:
        fig,ax=plt.subplots(figsize=(7.5,4.8))
        ax.boxplot([100*audit.loc[audit.MAPQ==q,'frac_multi_total'].dropna() for q in qs],tick_labels=[str(q) for q in qs],showfliers=False)
        ax.set(xlabel='Minimum MAPQ',ylabel='Multi-reference fraction (%)',title=f'{sci_t}: nuclear–organelle ambiguity across MAPQ')
        save(fig,f/'Figure_3_multi_by_MAPQ.png')

    if 'Model' in primary and primary.Model.nunique()>1:
        counts=primary.Model.value_counts(); models=counts[counts>=5].index.tolist()
        if models:
            fig,ax=plt.subplots(figsize=(9,5))
            ax.boxplot([100*primary.loc[primary.Model==m,'frac_strict_total'].dropna() for m in models],tick_labels=[f'{m}\n(n={counts[m]})' for m in models],showfliers=False)
            ax.tick_params(axis='x',rotation=25); ax.set(ylabel='Strict fraction (%)',title=f'{sci_t}: strict fraction by sequencing model')
            save(fig,f/'Figure_4_strict_by_model.png')

    if nuclear_available:
        cross=pd.crosstab(pd.Categorical(primary.decision_one_threshold,categories=ORDER,ordered=True),pd.Categorical(primary.decision_framework,categories=ORDER,ordered=True)).reindex(index=ORDER,columns=ORDER,fill_value=0)
        fig,ax=plt.subplots(figsize=(6.8,5)); im=ax.imshow(cross.values,cmap='cividis',aspect='auto')
        for i in range(3):
            for j in range(3): ax.text(j,i,str(int(cross.iloc[i,j])),ha='center',va='center',color='white' if cross.iloc[i,j]>cross.values.max()/2 else 'black')
        ax.set_xticks(range(3),ORDER,rotation=20,ha='right'); ax.set_yticks(range(3),ORDER); ax.set_xlabel('Framework category'); ax.set_ylabel('One-threshold category'); ax.set_title(f'{sci_t}: category changes after separating ambiguity'); fig.colorbar(im,ax=ax,label='Run count')
        save(fig,f/'Figure_5_triage_heatmap_cividis.png')

    if (t/'downstream_alignment_depth_deltas.tsv').exists():
        d=pd.read_csv(t/'downstream_alignment_depth_deltas.tsv',sep='\t')
        horizontal_two_condition(d,'delta_callable_positions','Change in callable positions vs baseline',f'{sci_t}: downstream depth sensitivity',f/'Figure_6_downstream_callable_delta.png')
    if (t/'downstream_variant_deltas.tsv').exists():
        d=pd.read_csv(t/'downstream_variant_deltas.tsv',sep='\t')
        horizontal_two_condition(d,'delta_variants','Change in called variants vs baseline',f'{sci_t}: variant-call sensitivity',f/'Figure_7_variant_delta.png')
        horizontal_two_condition(d,'jaccard_vs_baseline','Variant-site Jaccard vs baseline',f'{sci_t}: overlap of variant calls after filtering',f/'Figure_8_variant_jaccard.png',smaller_is_worse=True)
    if (t/'reference_mask_sensitivity.tsv').exists():
        r=pd.read_csv(t/'reference_mask_sensitivity.tsv',sep='\t'); piv=r.pivot(index='Run',columns='nuclear_reference',values='multi_qnames_primary')
        fig,ax=plt.subplots(figsize=(8.5,5)); piv.plot(kind='bar',ax=ax); ax.set_ylabel('Multi-reference qnames at primary MAPQ'); ax.set_title(f'{sci_t}: sensitivity to nuclear reference masking'); ax.tick_params(axis='x',rotation=90,labelsize=7); save(fig,f/'Figure_9_reference_mask_sensitivity.png')
        if 'n_variants' in r.columns:
            vp=r.pivot(index='Run',columns='nuclear_reference',values='n_variants'); fig,ax=plt.subplots(figsize=(8.5,5)); vp.plot(kind='bar',ax=ax); ax.set_ylabel('Called variants'); ax.set_title(f'{sci_t}: variant sensitivity to nuclear reference masking'); ax.tick_params(axis='x',rotation=90,labelsize=7); save(fig,f/'Figure_9b_reference_mask_variant_sensitivity.png')
    print(f'[OK] figures written for {sid}: {f}')

if __name__=='__main__': main()
