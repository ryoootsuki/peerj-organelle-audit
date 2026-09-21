#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd
from common import load_config, results_dir, workdir


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); args=ap.parse_args(); cfg=load_config(args.config); sid=cfg['species']['id']
    tdir=results_dir(cfg)/'tables'; subset=pd.read_csv(tdir/'downstream_subset.tsv',sep='\t'); rows=[]; vars=[]; refs=[]; status=[]
    for run in subset.Run.drop_duplicates():
        rd=workdir(cfg)/'downstream'/run
        if (rd/'status.json').exists(): status.append(json.loads((rd/'status.json').read_text()))
        if (rd/'downstream_metrics.tsv').exists(): rows.append(pd.read_csv(rd/'downstream_metrics.tsv',sep='\t'))
        if (rd/'variant_metrics.tsv').exists(): vars.append(pd.read_csv(rd/'variant_metrics.tsv',sep='\t'))
        if (rd/'reference_sensitivity.tsv').exists(): refs.append(pd.read_csv(rd/'reference_sensitivity.tsv',sep='\t'))
    pd.DataFrame(status).to_csv(results_dir(cfg)/'provenance'/'downstream_completion_status.tsv',sep='\t',index=False)
    if not rows or not vars or not refs: raise SystemExit('[ERROR] incomplete downstream outputs; metrics, variant, or reference sensitivity tables are absent')
    m=pd.concat(rows,ignore_index=True); v=pd.concat(vars,ignore_index=True); r=pd.concat(refs,ignore_index=True)
    group=(subset.groupby('Run').agg(selection_group=('selection_group',lambda x:'|'.join(sorted(set(x)))),selection_reason=('selection_reason',lambda x:' | '.join(dict.fromkeys(map(str,x))))).reset_index())
    m=m.merge(group,on='Run',how='left'); v=v.merge(group,on='Run',how='left'); r=r.merge(group,on='Run',how='left')
    m.to_csv(tdir/'downstream_alignment_depth_metrics.tsv',sep='\t',index=False); v.to_csv(tdir/'downstream_variant_metrics.tsv',sep='\t',index=False); r.to_csv(tdir/'reference_mask_sensitivity.tsv',sep='\t',index=False)

    base_m=m.loc[m.condition=='baseline'].set_index('Run')
    md=[]
    for _,x in m.loc[m.condition!='baseline'].iterrows():
        b=base_m.loc[x.Run]
        md.append({'Run':x.Run,'selection_group':x.selection_group,'condition':x.condition,'removed_primary_qnames':x.removed_primary_qnames,'delta_callable_positions':x.callable_positions_ge_min_depth-b.callable_positions_ge_min_depth,'delta_covered_positions':x.covered_positions-b.covered_positions,'delta_mean_depth_covered':x.mean_depth_covered-b.mean_depth_covered})
    pd.DataFrame(md).to_csv(tdir/'downstream_alignment_depth_deltas.tsv',sep='\t',index=False)
    base_v=v.loc[v.condition=='baseline'].set_index('Run')
    vd=[]
    for _,x in v.loc[v.condition!='baseline'].iterrows():
        b=base_v.loc[x.Run]
        vd.append({'Run':x.Run,'selection_group':x.selection_group,'condition':x.condition,'delta_variants':x.n_variants-b.n_variants,'delta_snps':x.n_snps-b.n_snps,'delta_indels':x.n_indels-b.n_indels,'jaccard_vs_baseline':x.jaccard_vs_baseline,'sites_lost_vs_baseline':x.sites_lost_vs_baseline,'sites_gained_vs_baseline':x.sites_gained_vs_baseline})
    pd.DataFrame(vd).to_csv(tdir/'downstream_variant_deltas.tsv',sep='\t',index=False)


    ref_base=r.loc[r.nuclear_reference=='nuclear_masked_main'].set_index('Run')
    ref_delta=[]
    for _,x in r.loc[r.nuclear_reference!='nuclear_masked_main'].iterrows():
        b=ref_base.loc[x.Run]
        ref_delta.append({'Run':x.Run,'selection_group':x.selection_group,'nuclear_reference':x.nuclear_reference,'delta_nuclear_qnames_primary':x.nuclear_qnames_primary-b.nuclear_qnames_primary,'delta_multi_qnames_primary':x.multi_qnames_primary-b.multi_qnames_primary,'delta_callable_positions':x.callable_positions_ge_min_depth-b.callable_positions_ge_min_depth,'delta_variants':x.n_variants-b.n_variants,'variant_jaccard_vs_masked_main':x.variant_jaccard_vs_masked_main,'variant_sites_lost_vs_masked_main':x.variant_sites_lost_vs_masked_main,'variant_sites_gained_vs_masked_main':x.variant_sites_gained_vs_masked_main})
    pd.DataFrame(ref_delta).to_csv(tdir/'reference_mask_variant_deltas.tsv',sep='\t',index=False)

    ref_piv=r.pivot(index=['Run','selection_group'],columns='nuclear_reference',values=['nuclear_qnames_primary','multi_qnames_primary','strict_organelle_qnames_primary']).reset_index()
    ref_piv.columns=['__'.join([str(y) for y in x if str(y)]) if isinstance(x,tuple) else x for x in ref_piv.columns]
    ref_piv.to_csv(tdir/'reference_mask_sensitivity_wide.tsv',sep='\t',index=False)
    print(f'[OK] downstream outputs aggregated for {sid}: runs={subset.Run.nunique()}')
if __name__=='__main__': main()
