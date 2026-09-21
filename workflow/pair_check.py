#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip
from itertools import zip_longest
from common import normalize_qname

def opener(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p,'r',encoding='utf-8')
def names(path):
    with opener(path) as f:
        i=0
        while True:
            h=f.readline()
            if not h: break
            s=f.readline(); plus=f.readline(); q=f.readline()
            if not (s and plus and q): raise SystemExit(f'[ERROR] truncated FASTQ: {path}')
            if not h.startswith('@'): raise SystemExit(f'[ERROR] invalid FASTQ header: {path}: {h[:80]}')
            yield normalize_qname(h[1:]); i+=1

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--r1',required=True); ap.add_argument('--r2',required=True); args=ap.parse_args()
    n=0
    for a,b in zip_longest(names(args.r1),names(args.r2)):
        if a is None or b is None: raise SystemExit('[ERROR] paired FASTQ record counts differ')
        if a!=b: raise SystemExit(f'[ERROR] paired FASTQ names differ at pair {n+1}: {a} != {b}')
        n+=1
    if n == 0: raise SystemExit('[ERROR] synchronized paired FASTQ contains zero pairs')
    print(f'[OK] synchronized paired FASTQ: {n} pairs')
if __name__=='__main__': main()
