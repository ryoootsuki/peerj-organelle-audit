#!/usr/bin/env python3
"""Check key retained numerical summaries against compact result evidence."""
from pathlib import Path
import math
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
E=ROOT/'examples/results_summary'

def close(a,b,tol):
    assert abs(float(a)-float(b)) <= tol, (a,b,tol)

# Silene MAPQ 20 claims.
m=pd.read_csv(E/'silene_mapq_distribution_summary.tsv',sep='\t')
def row(metric):
    q=m[(m.MAPQ==20)&(m.metric==metric)]
    assert len(q)==1
    return q.iloc[0]
s=row('frac_strict_total'); i=row('frac_inclusive_total'); u=row('frac_multi_total')
assert int(s.n_runs)==2679
close(100*s['median'],2.381,0.001)
close(100*i['median'],2.894,0.001)
close(100*u['median'],0.458,0.001)
close(100*s['p95'],21.25,0.01)
close(100*s['maximum'],32.15,0.01)

# Primary triage setting.
g=pd.read_csv(E/'silene_triage_threshold_sensitivity.tsv',sep='\t')
q=g[(g.low_fraction==0.05)&(g.high_fraction==0.20)&(g.ambiguity_fraction==0.01)]
assert len(q)==1
r=q.iloc[0]
assert int(r.n_runs)==2679 and int(r.reclassified_vs_one_threshold)==208
close(100*r.reclassified_fraction,7.76,0.01)
assert (int(r.retain),int(r.sensitivity_analysis),int(r.exclude))==(2002,519,158)

# Paired/single descriptive association in the current manuscript.
t=pd.read_csv(E/'silene_technical_group_comparisons_curated.tsv',sep='\t')
p=t[t.comparison.eq('paired_vs_single')].iloc[0]
close(100*p.median_1,21.10,0.01)
close(100*p.median_2,2.206,0.001)
close(p.cliffs_delta,0.743,0.001)

# Downstream Silene claims.
d=pd.read_csv(E/'silene_downstream_variant_deltas.tsv',sep='\t')
assert d.Run.nunique()==26 and len(d)==52
assert int(d.delta_variants.astype(float).abs().max())==4644
close(d.jaccard_vs_baseline.astype(float).min(),0.989992,0.000001)
assert not d.astype(str).apply(lambda c:c.str.contains('wgs_contrast',case=False,na=False)).any().any()

# Chicken claims.
cm=pd.read_csv(E/'chicken_mapq_distribution_summary.tsv',sep='\t')
cr=cm[(cm.MAPQ==20)&(cm.metric=='frac_inclusive_total')].iloc[0]
close(100*cr.maximum,0.00664,0.00001)
cd=pd.read_csv(E/'chicken_downstream_variant_deltas.tsv',sep='\t')
assert int(cd.delta_variants.astype(float).abs().max())==1

print('[OK] manuscript numerical concordance checks passed')
