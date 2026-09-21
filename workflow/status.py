#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from common import load_config, species_id, workdir


def summarize(cfg_path: str) -> None:
    cfg = load_config(cfg_path)
    sid = species_id(cfg)
    wd = workdir(cfg)
    for stage, base, pattern in [
        ("audit", wd / "runs", "audit.status.json"),
        ("downstream", wd / "downstream", "status.json"),
    ]:
        c = Counter()
        stage_c = Counter()
        for p in base.glob(f"*/{pattern}") if base.exists() else []:
            try:
                d = json.loads(p.read_text())
                c[str(d.get("status", "unknown"))] += 1
                stage_c[str(d.get("stage", "unknown"))] += 1
            except Exception:
                c["unreadable"] += 1
        print(f"[{sid}] {stage}: " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())) if c else f"[{sid}] {stage}: no status files")
        if stage_c:
            print("  stages: " + ", ".join(f"{k}={v}" for k, v in sorted(stage_c.items())))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="*", default=["config/species/silene.yaml", "config/species/chicken.yaml"])
    args = ap.parse_args()
    for c in args.configs:
        summarize(c)


if __name__ == "__main__":
    main()
