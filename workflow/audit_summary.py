#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from common import load_config, results_dir, workdir


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)==0 or len(b)==0: return float('nan')
    # Rank-based identity: delta = 2U/(n1*n2)-1.
    u=mannwhitneyu(a,b,alternative='two-sided').statistic
    return float(2*u/(len(a)*len(b))-1)


def bootstrap_median_difference(a,b,n,seed):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)==0 or len(b)==0: return (np.nan,np.nan,np.nan)
    rng=np.random.default_rng(seed); vals=np.empty(n)
    for i in range(n): vals[i]=np.median(rng.choice(a,len(a),replace=True))-np.median(rng.choice(b,len(b),replace=True))
    return float(np.median(a)-np.median(b)),float(np.quantile(vals,.025)),float(np.quantile(vals,.975))


def grouped(df,col):
    rows=[]
    for key,sub in df.groupby(col,dropna=False):
        rows.append({
            col:key,'n_runs':len(sub),'median_strict':sub.frac_strict_total.median(),'iqr_strict':sub.frac_strict_total.quantile(.75)-sub.frac_strict_total.quantile(.25),
            'p95_strict':sub.frac_strict_total.quantile(.95),'max_strict':sub.frac_strict_total.max(),
            'median_inclusive':sub.frac_inclusive_total.median(),'median_multi':sub.frac_multi_total.median(),'median_spots':sub.spots.median()
        })
    return pd.DataFrame(rows).sort_values(['n_runs',col],ascending=[False,True]) if rows else pd.DataFrame()


def fastp_row(run_dir: Path, run: str):
    p=run_dir/'fastp'/f'{run}.json'
    if not p.exists(): return {'Run':run,'fastp_status':'missing'}
    try:
        x=json.loads(p.read_text()); s=x.get('summary',{}); before=s.get('before_filtering',{}); after=s.get('after_filtering',{})
        return {'Run':run,'fastp_status':'parsed','reads_before':before.get('total_reads'),'bases_before':before.get('total_bases'),'q20_rate_before':before.get('q20_rate'),'q30_rate_before':before.get('q30_rate'),'reads_after':after.get('total_reads'),'bases_after':after.get('total_bases'),'q20_rate_after':after.get('q20_rate'),'q30_rate_after':after.get('q30_rate')}
    except Exception as exc: return {'Run':run,'fastp_status':f'parse_error:{exc}'}


def category_one(x,low,high): return 'Retain' if x<low else ('Sensitivity analysis' if x<high else 'Exclude')
def category_framework(strict,multi,low,high,multi_cut):
    if strict>=high: return 'Exclude'
    if strict<low and multi<multi_cut: return 'Retain'
    return 'Sensitivity analysis'


def compare_two_groups(df,group_col,g1,g2,metric,cfg,label):
    a=df.loc[df[group_col]==g1,metric].dropna().to_numpy(); b=df.loc[df[group_col]==g2,metric].dropna().to_numpy()
    if len(a)==0 or len(b)==0: return None
    u=mannwhitneyu(a,b,alternative='two-sided'); diff,lo,hi=bootstrap_median_difference(a,b,int(cfg['statistics']['bootstrap_replicates']),int(cfg['workflow']['random_seed']))
    return {'comparison':label,'group_column':group_col,'group_1':g1,'group_2':g2,'metric':metric,'n_1':len(a),'n_2':len(b),'median_1':np.median(a),'median_2':np.median(b),'median_difference_1_minus_2':diff,'bootstrap_95ci_low':lo,'bootstrap_95ci_high':hi,'mann_whitney_u':u.statistic,'p_value_unadjusted':u.pvalue,'cliffs_delta':cliffs_delta(a,b),'interpretation_note':'Descriptive association only; project, protocol, tissue, depth, and platform are potential confounders.'}


def benjamini_hochberg(values):
    import numpy as np
    p=np.asarray(values,dtype=float); n=len(p)
    if n==0: return []
    order=np.argsort(p); ranked=p[order]; adj=np.empty(n,float); running=1.0
    for i in range(n-1,-1,-1):
        running=min(running, ranked[i]*n/(i+1)); adj[order[i]]=min(running,1.0)
    return adj.tolist()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); args=ap.parse_args(); cfg=load_config(args.config); sid=cfg['species']['id']
    rdir=results_dir(cfg); tdir=rdir/'tables'; pdir=rdir/'provenance'; nuclear_available=bool(cfg['references'].get('nuclear_assembly_accession')); tdir.mkdir(parents=True,exist_ok=True); pdir.mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(workdir(cfg)/'manifest/manifest.all.tsv',sep='\t')
    audit_frames=[]; statuses=[]; fastp=[]
    for _,m in manifest.iterrows():
        run=m.Run; rd=workdir(cfg)/'runs'/run; sp=rd/'audit.status.json'
        if sp.exists():
            try: statuses.append(json.loads(sp.read_text()))
            except Exception as exc: statuses.append({'Run':run,'status':'status_parse_error','error':str(exc)})
        else:
            ev = str(m.evaluable).strip().lower() in {'true','1','yes'}
            statuses.append({'Run':run,'status':'not_started' if ev else 'excluded_before_download','error':m.exclusion_reason})
        apath=rd/'audit_rows.tsv'
        if apath.exists() and apath.stat().st_size: audit_frames.append(pd.read_csv(apath,sep='\t'))
        fastp.append(fastp_row(rd,run))
    status=pd.DataFrame(statuses); status.to_csv(pdir/'run_completion_status.tsv',sep='\t',index=False)
    pd.DataFrame(fastp).merge(manifest,on='Run',how='left').to_csv(tdir/'fastp_run_metrics.tsv',sep='\t',index=False)
    if not audit_frames: raise SystemExit('[ERROR] no completed audit tables found')
    audit=pd.concat(audit_frames,ignore_index=True)
    metadata_cols=[c for c in manifest.columns if c not in audit.columns or c=='Run']
    audit=audit.merge(manifest[metadata_cols],on='Run',how='left',validate='many_to_one')
    audit.to_csv(tdir/'audit_mapq_sweep.tsv',sep='\t',index=False)
    primary=audit.loc[audit.MAPQ==int(cfg['workflow']['primary_mapq'])].copy()
    primary.to_csv(tdir/'audit_primary_mapq.tsv',sep='\t',index=False)

    mapq_rows=[]
    for q,sub in audit.groupby('MAPQ'):
        for metric in ['frac_strict_total','frac_inclusive_total','frac_multi_total']:
            s=sub[metric].dropna()
            if len(s)==0: continue
            mapq_rows.append({'MAPQ':q,'metric':metric,'n_runs':len(s),'median':s.median(),'q1':s.quantile(.25),'q3':s.quantile(.75),'iqr':s.quantile(.75)-s.quantile(.25),'p95':s.quantile(.95),'maximum':s.max()})
    pd.DataFrame(mapq_rows).to_csv(tdir/'mapq_distribution_summary.tsv',sep='\t',index=False)

    for col in ['BioProject','Model','LibraryLayout','LibraryStrategy','Platform']:
        if col in primary.columns: grouped(primary,col).to_csv(tdir/f'grouped_by_{col}.tsv',sep='\t',index=False)

    tri=cfg['triage']; low=float(tri['low_fraction']); high=float(tri['high_fraction']); mc=float(tri['ambiguity_fraction'])
    primary['decision_one_threshold']=primary.frac_inclusive_total.map(lambda x:category_one(x,low,high))
    if nuclear_available:
        primary['decision_framework']=primary.apply(lambda x:category_framework(x.frac_strict_total,x.frac_multi_total,low,high,mc),axis=1)
        primary['reclassified']=primary.decision_one_threshold!=primary.decision_framework
    else:
        primary['decision_framework']='Not available: nuclear reference absent'
        primary['reclassified']=False
    primary.to_csv(tdir/'one_threshold_vs_framework.tsv',sep='\t',index=False)
    primary.to_csv(tdir/'audit_primary_mapq.tsv',sep='\t',index=False)

    tests=[]
    comparison_metric='frac_strict_total' if nuclear_available else 'frac_inclusive_total'
    if 'LibraryLayout' in primary:
        x=compare_two_groups(primary,'LibraryLayout','PAIRED','SINGLE',comparison_metric,cfg,'paired_vs_single');
        if x: tests.append(x)
    if 'LibraryStrategy' in primary:
        strategies=set(primary.LibraryStrategy.astype(str))
        if 'WGS' in strategies and 'RAD-Seq' in strategies:
            x=compare_two_groups(primary,'LibraryStrategy','WGS','RAD-Seq',comparison_metric,cfg,'wgs_vs_radseq');
            if x: tests.append(x)
            if nuclear_available:
                x=compare_two_groups(primary,'LibraryStrategy','WGS','RAD-Seq','frac_multi_total',cfg,'wgs_vs_radseq_ambiguity');
                if x: tests.append(x)
    comparison_columns=['comparison','group_column','group_1','group_2','metric','n_1','n_2','median_1','median_2','median_difference_1_minus_2','bootstrap_95ci_low','bootstrap_95ci_high','mann_whitney_u','p_value_unadjusted','cliffs_delta','interpretation_note','p_value_bh']
    test_df=pd.DataFrame(tests)
    if len(test_df):
        test_df['p_value_bh']=benjamini_hochberg(test_df['p_value_unadjusted'].tolist())
        test_df=test_df.reindex(columns=comparison_columns)
    else:
        test_df=pd.DataFrame(columns=comparison_columns)
    test_df.to_csv(tdir/'technical_group_comparisons.tsv',sep='\t',index=False)

    expected=int(cfg['species']['expected_evaluable_runs']); observed=primary.Run.nunique(); failed=max(expected-observed,0)
    completion=pd.DataFrame([{'species':sid,'expected_total_runs':int(cfg['species']['expected_total_runs']),'expected_evaluable_runs':expected,'observed_primary_mapq_runs':observed,'missing_or_failed_evaluable_runs':failed,'full_audit_complete':observed==expected,'primary_mapq':int(cfg['workflow']['primary_mapq'])}])
    completion.to_csv(pdir/'analysis_completion_summary.tsv',sep='\t',index=False)
    print(f'[OK] aggregate audit complete for {sid}: primary runs={observed}/{expected}')

if __name__=='__main__': main()
