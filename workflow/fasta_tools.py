#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
from common import iter_fasta, write_fasta_record


def cmd_header_ids(args):
    explicit={x.split('.')[0] for x in args.accessions}
    terms=[x.lower() for x in args.keywords]
    with open(args.output,'w',encoding='utf-8') as out:
        for header,_ in iter_fasta(args.fasta):
            sid=header.split()[0]
            h=header.lower()
            if sid.split('.')[0] in explicit or any(t in h for t in terms):
                out.write(sid+'\n')


def cmd_double(args):
    with open(args.output,'w',encoding='utf-8') as out:
        for header,seq in iter_fasta(args.fasta):
            write_fasta_record(out,header.split()[0]+'_double',seq+seq)


def cmd_mask_summary(args):
    total=0; intervals=0
    if Path(args.bed).exists():
        with open(args.bed,encoding='utf-8') as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'): continue
                f=line.split('\t'); total+=int(f[2])-int(f[1]); intervals+=1
    Path(args.output).write_text(f'mask_name\tn_intervals\tmasked_bp\n{args.name}\t{intervals}\t{total}\n',encoding='utf-8')


def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='cmd',required=True)
    p=sub.add_parser('header-ids'); p.add_argument('--fasta',required=True); p.add_argument('--output',required=True); p.add_argument('--accessions',nargs='*',default=[]); p.add_argument('--keywords',nargs='*',default=['mitochondrion','chloroplast','plastid']); p.set_defaults(func=cmd_header_ids)
    p=sub.add_parser('double'); p.add_argument('--fasta',required=True); p.add_argument('--output',required=True); p.set_defaults(func=cmd_double)
    p=sub.add_parser('mask-summary'); p.add_argument('--bed',required=True); p.add_argument('--name',required=True); p.add_argument('--output',required=True); p.set_defaults(func=cmd_mask_summary)
    args=ap.parse_args(); args.func(args)
if __name__=='__main__': main()
