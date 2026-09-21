#!/usr/bin/env python3
from pathlib import Path
import re, sys
ROOT=Path(__file__).resolve().parents[1]
errors=[]
for p in ROOT.rglob('*'):
    if p.resolve() == Path(__file__).resolve(): continue
    if not p.is_file() or p.suffix.lower() in {'.png','.pdf','.zip','.gz'}: continue
    try: s=p.read_text(encoding='utf-8',errors='ignore')
    except Exception: continue
    for pat in [r'/media/ryo/',r'/mnt/20TB',r'AUTHOR ACTION REQUIRED',r'\[PAGE/LINE\]',r'PENDING_OR_CONFIRM',r'Current Draft',r'Remaining author actions',r'release blocker',r'reviewer-only',r'editor-only',r'manual confirmation from supplementary methods required',r'\breviewer\b']:
        if re.search(pat,s,flags=re.IGNORECASE): errors.append(f'{p.relative_to(ROOT)} contains disallowed pattern {pat}')
# required scientific assets
for rel in ['metadata/bioproject_protocol_metadata_curated.tsv','examples/results_summary/silene_triage_threshold_sensitivity.tsv','examples/results_summary/silene_mapq_distribution_summary.tsv','figures/manuscript/silene/Figure_5_triage_heatmap_cividis.png']:
    if not (ROOT/rel).is_file(): errors.append(f'missing {rel}')
if errors:
    print('[FAIL] public safety audit')
    for e in errors: print(' -',e)
    raise SystemExit(1)
print('[OK] public safety audit passed')
