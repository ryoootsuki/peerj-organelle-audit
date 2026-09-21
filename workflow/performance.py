#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd


def parse_time(x: str | None):
    if not x:
        return None
    try:
        return datetime.fromisoformat(x.replace("Z", "+00:00"))
    except Exception:
        return None


def collect(species: str, stage: str, base: Path, filename: str) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    starts = []
    ends = []
    for p in base.glob(f"*/{filename}") if base.exists() else []:
        try:
            d = json.loads(p.read_text())
        except Exception:
            continue
        if d.get("status") != "ready":
            continue
        wall = float(d.get("wall_seconds", 0) or 0)
        rows.append({
            "species": species,
            "analysis_stage": stage,
            "Run": d.get("Run", p.parent.name),
            "wall_seconds": wall,
            "layout": d.get("layout", ""),
            "spots": d.get("spots", ""),
        })
        st = parse_time(d.get("started_utc"))
        en = parse_time(d.get("completed_utc"))
        if st:
            starts.append(st)
        if en:
            ends.append(en)
    if rows:
        s = pd.Series([r["wall_seconds"] for r in rows], dtype=float)
        elapsed = (max(ends) - min(starts)).total_seconds() if starts and ends else 0.0
        summary = {
            "species": species,
            "analysis_stage": stage,
            "n_completed": len(rows),
            "median_run_wall_seconds": s.median(),
            "p95_run_wall_seconds": s.quantile(0.95),
            "maximum_run_wall_seconds": s.max(),
            "sum_run_wall_hours": s.sum() / 3600,
            "observed_stage_elapsed_hours": elapsed / 3600 if elapsed > 0 else 0.0,
            "observed_runs_per_hour": len(rows) / (elapsed / 3600) if elapsed > 0 else 0.0,
        }
    else:
        summary = {
            "species": species, "analysis_stage": stage, "n_completed": 0,
            "median_run_wall_seconds": 0, "p95_run_wall_seconds": 0,
            "maximum_run_wall_seconds": 0, "sum_run_wall_hours": 0,
            "observed_stage_elapsed_hours": 0, "observed_runs_per_hour": 0,
        }
    return rows, summary


def main() -> None:
    all_rows: list[dict] = []
    summaries: list[dict] = []
    for sid in ["silene", "chicken"]:
        r, s = collect(sid, "archive_audit", Path("work") / sid / "runs", "audit.status.json")
        all_rows.extend(r)
        summaries.append(s)
        r, s = collect(sid, "downstream", Path("work") / sid / "downstream", "status.json")
        all_rows.extend(r)
        summaries.append(s)
    out = Path("results/combined/provenance")
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(out / "run_wall_times.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(out / "performance_summary.tsv", sep="\t", index=False)
    print(f"[OK] performance summary written: {out / 'performance_summary.tsv'}")


if __name__ == "__main__":
    main()
