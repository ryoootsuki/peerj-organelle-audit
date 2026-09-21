#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

from common import normalize_qname

EXCLUDE_FLAGS = 0x4 | 0x100 | 0x800  # unmapped, secondary, supplementary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--thresholds", required=True, nargs="+", type=int)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--compression-level", type=int, default=1)
    args = ap.parse_args()

    max_mapq: dict[str, int] = {}
    n_records = 0
    n_primary_mapped_records = 0
    malformed = 0

    for line in sys.stdin:
        if not line or line.startswith("@"):
            continue
        n_records += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            malformed += 1
            continue
        try:
            flag = int(f[1])
            mapq = int(f[4])
        except ValueError:
            malformed += 1
            continue
        if flag & EXCLUDE_FLAGS:
            continue
        n_primary_mapped_records += 1
        qname = normalize_qname(f[0])
        old = max_mapq.get(qname)
        if old is None or mapq > old:
            max_mapq[qname] = mapq

    items = sorted(max_mapq.items())
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    counts: dict[int, int] = {}
    for q in sorted(set(args.thresholds)):
        path = outdir / f"{args.tag}.q{q}.qnames.txt.gz"
        count = 0
        with gzip.open(path, "wt", compresslevel=max(1, min(9, args.compression_level))) as out:
            for name, mq in items:
                if mq >= q:
                    out.write(name + "\n")
                    count += 1
        counts[q] = count

    summary = {
        "tag": args.tag,
        "sam_records_seen": n_records,
        "primary_mapped_records": n_primary_mapped_records,
        "unique_primary_mapped_qnames": len(max_mapq),
        "malformed_records": malformed,
        "threshold_counts": counts,
        "algorithm": "single_pass_max_mapq_per_normalized_qname",
        "excluded_flags": EXCLUDE_FLAGS,
    }
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[OK] {args.tag}: {len(max_mapq)} unique primary mapped qnames", file=sys.stderr)


if __name__ == "__main__":
    main()
