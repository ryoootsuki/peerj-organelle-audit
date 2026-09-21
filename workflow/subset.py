#!/usr/bin/env python3
from __future__ import annotations
import argparse
from collections import Counter
import numpy as np
import pandas as pd
from common import load_config, results_dir, workdir


def diverse_top(df: pd.DataFrame, n: int, sort_cols, ascending, group: str, reason: str, max_per_project: int = 2) -> pd.DataFrame:
    ranked=df.sort_values(sort_cols,ascending=ascending).copy()
    chosen=[]; counts=Counter()
    for _,row in ranked.iterrows():
        project=str(row.get('BioProject','not_reported'))
        if counts[project] >= max_per_project: continue
        chosen.append(row); counts[project]+=1
        if len(chosen)>=n: break
    if len(chosen)<n:
        already={x.Run for x in chosen}
        for _,row in ranked.loc[~ranked.Run.isin(already)].iterrows():
            chosen.append(row)
            if len(chosen)>=n: break
    out=pd.DataFrame(chosen)
    if len(out): out['selection_group']=group; out['selection_reason']=reason; out['selection_rank']=range(1,len(out)+1)
    return out


def matched_low_controls(df: pd.DataFrame, anchors: pd.DataFrame, n: int, excluded: set[str]) -> pd.DataFrame:
    pool=df.loc[(df.frac_strict_total<=df.frac_strict_total.quantile(.25)) & (df.frac_multi_total<=df.frac_multi_total.quantile(.25)) & ~df.Run.isin(excluded)].copy()
    pool['log_spots']=np.log10(pool.spots.clip(lower=1)); selected=[]
    for _,a in anchors.head(n).iterrows():
        cand=pool.loc[~pool.Run.isin([x.Run for x in selected])].copy()
        for col in ['BioProject','LibraryLayout','LibraryStrategy']:
            if col in cand and col in a.index:
                same=cand.loc[cand[col].astype(str)==str(a[col])]
                if len(same): cand=same
        if cand.empty: continue
        cand['distance']=abs(cand.log_spots-np.log10(max(float(a.spots),1)))
        selected.append(cand.sort_values(['distance','frac_strict_total','frac_multi_total','Run']).iloc[0])
    if len(selected)<n:
        used={x.Run for x in selected}; extra=pool.loc[~pool.Run.isin(used)].sort_values(['frac_strict_total','frac_multi_total','Run']).head(n-len(selected))
        selected.extend([r for _,r in extra.iterrows()])
    out=pd.DataFrame(selected)
    if len(out): out['selection_group']='low_control'; out['selection_reason']='Low strict/ambiguity controls matched where possible by BioProject, layout, strategy, and log10(spots).'; out['selection_rank']=range(1,len(out)+1)
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); args=ap.parse_args(); cfg=load_config(args.config); sid=cfg['species']['id']
    tdir=results_dir(cfg)/'tables'
    if not cfg['references'].get('nuclear_assembly_accession'):
        raise SystemExit('[ERROR] downstream filtering/variant analysis requires a nuclear reference; organelle-only mode supports inclusive archive screening only')
    df=pd.read_csv(tdir/'one_threshold_vs_framework.tsv',sep='\t')
    n=int(cfg['subset']['silene_runs_per_group'] if sid=='silene' else cfg['subset']['chicken_runs_per_group'])
    groups=[]
    if sid=='silene':
        hs=diverse_top(df,n,['frac_strict_total','frac_multi_total','Run'],[False,False,True],'high_strict','Highest strict organelle-compatible fractions at the prespecified primary MAPQ, with project diversity where possible.')
        hm=diverse_top(df,n,['frac_multi_total','frac_strict_total','Run'],[False,False,True],'high_multi','Highest nuclear-organelle multi-reference fractions at the prespecified primary MAPQ, with project diversity where possible.')
        re=df.loc[df.reclassified].copy(); re['framework_difference']=(re.frac_inclusive_total-re.frac_strict_total).abs()
        rr=diverse_top(re,n,['framework_difference','frac_multi_total','Run'],[False,False,True],'reclassified','Runs whose triage category changed after strict signal and ambiguity were reported separately.')
        anchors=pd.concat([hs,hm],ignore_index=True).drop_duplicates('Run')
        low=matched_low_controls(df,anchors,n,set(anchors.Run)|set(rr.Run if len(rr) else []))
        prj=df.loc[df.BioProject.astype(str).eq('PRJNA285775')].copy()
        prj_group=diverse_top(prj,n,['frac_strict_total','frac_multi_total','Run'],[False,False,True],'prjna285775_sbfI_rad_contrast','SbfI RAD-seq contrast runs from PRJNA285775. This replaces the invalid WGS label after publication-based metadata curation.',max_per_project=n)
        groups=[hs,hm,rr,low,prj_group]
    else:
        hi=diverse_top(df,n,['frac_inclusive_total','frac_multi_total','Run'],[False,False,True],'high_organelle_compatible','Highest mitochondrial-compatible fractions within PRJNA573756; “high” is relative to this low-signal dataset.',max_per_project=n)
        hm=diverse_top(df,n,['frac_multi_total','frac_inclusive_total','Run'],[False,False,True],'high_multi','Highest nuclear-mitochondrial ambiguity fractions within PRJNA573756; overlap with the compatible group is retained transparently.',max_per_project=n)
        low=matched_low_controls(df,pd.concat([hi,hm]).drop_duplicates('Run'),n,set(hi.Run)|set(hm.Run))
        groups=[hi,hm,low]
    sel=pd.concat([x for x in groups if x is not None and len(x)],ignore_index=True)
    if sel.empty: raise SystemExit('[ERROR] no downstream runs could be selected')
    keep=['Run','selection_group','selection_rank','selection_reason','spots','LibraryStrategy','LibraryLayout','Platform','Model','BioProject','frac_strict_total','frac_inclusive_total','frac_multi_total','decision_one_threshold','decision_framework','reclassified']
    keep=[c for c in keep if c in sel.columns]
    sel[keep].to_csv(tdir/'downstream_subset.tsv',sep='\t',index=False)
    membership=(sel.groupby('Run').agg(selection_groups=('selection_group',lambda x:'|'.join(sorted(set(x)))),n_group_memberships=('selection_group','nunique')).reset_index())
    membership.to_csv(tdir/'downstream_unique_runs.tsv',sep='\t',index=False)
    out=workdir(cfg)/'manifest'/'downstream_runs.txt'; out.write_text('\n'.join(sorted(membership.Run))+'\n',encoding='utf-8')
    print(f'[OK] downstream subset selected for {sid}: memberships={len(sel)}, unique runs={len(membership)}')
if __name__=='__main__': main()
