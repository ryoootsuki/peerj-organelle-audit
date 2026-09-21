#!/usr/bin/env python3
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path


def merge_intervals(intervals):
    intervals=sorted(intervals)
    out=[]
    for s,e in intervals:
        if not out or s>out[-1][1]: out.append([s,e])
        else: out[-1][1]=max(out[-1][1],e)
    return [(s,e) for s,e in out]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--paf',required=True)
    ap.add_argument('--whole-contigs')
    ap.add_argument('--mask-bed')
    ap.add_argument('--whole-min-coverage',type=float,default=0.80)
    ap.add_argument('--whole-min-identity',type=float,default=0.95)
    ap.add_argument('--mask-min-length',type=int,default=200)
    ap.add_argument('--mask-min-identity',type=float,default=0.90)
    args=ap.parse_args()
    intervals=defaultdict(list); lengths={}; identity_num=defaultdict(float); aln_bp=defaultdict(int)
    mask_rows=[]
    with open(args.paf,encoding='utf-8') as fh:
        for line in fh:
            if not line.strip(): continue
            f=line.rstrip('\n').split('\t')
            if len(f)<12: continue
            target=f[5]; tlen=int(f[6]); start=int(f[7]); end=int(f[8]); matches=int(f[9]); alen=int(f[10])
            if alen<=0: continue
            ident=matches/alen
            lengths[target]=tlen
            intervals[target].append((min(start,end),max(start,end)))
            identity_num[target]+=matches; aln_bp[target]+=alen
            if alen>=args.mask_min_length and ident>=args.mask_min_identity:
                mask_rows.append((target,min(start,end),max(start,end),ident,alen,f[0]))
    if args.whole_contigs:
        p=Path(args.whole_contigs); p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('w',encoding='utf-8') as out:
            out.write('contig\tcontig_length\tmerged_aligned_bp\tcoverage_fraction\tweighted_identity\tremove\n')
            for target in sorted(lengths):
                merged=merge_intervals(intervals[target]); cov=sum(e-s for s,e in merged); frac=cov/lengths[target] if lengths[target] else 0
                ident=identity_num[target]/aln_bp[target] if aln_bp[target] else 0
                remove=frac>=args.whole_min_coverage and ident>=args.whole_min_identity
                out.write(f'{target}\t{lengths[target]}\t{cov}\t{frac:.8f}\t{ident:.8f}\t{str(remove).lower()}\n')
    if args.mask_bed:
        p=Path(args.mask_bed); p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('w',encoding='utf-8') as out:
            for target,s,e,ident,alen,query in sorted(mask_rows,key=lambda x:(x[0],x[1],x[2])):
                out.write(f'{target}\t{s}\t{e}\t{query};identity={ident:.6f};aln_len={alen}\n')
if __name__=='__main__': main()
