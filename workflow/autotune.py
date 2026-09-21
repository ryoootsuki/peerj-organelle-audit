#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
from pathlib import Path

import yaml


def physical_cores() -> int:
    try:
        p = subprocess.run(["lscpu", "-p=CORE,SOCKET"], check=True, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        pairs = {line for line in p.stdout.splitlines() if line and not line.startswith("#")}
        if pairs:
            return len(pairs)
    except Exception:
        pass
    return max(1, (os.cpu_count() or 1) // 2)


def memory_gb() -> float:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) / 1024 / 1024
    except Exception:
        pass
    return 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="config/runtime.auto.yaml")
    ap.add_argument("--report", default="results/combined/provenance/hardware_profile.tsv")
    args = ap.parse_args()

    logical = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    physical = physical_cores()
    mem = memory_gb()

    # Conservative throughput-oriented defaults. Scientific output is invariant
    # to these settings; they only control concurrency and temporary compression.
    audit_workers = max(1, min(8, physical // 8, logical // 12))
    audit_threads = max(4, min(10, logical // max(1, audit_workers) - 4))
    downstream_workers = max(1, min(4, physical // 16, logical // 24))
    downstream_threads = max(8, min(16, logical // max(1, downstream_workers) - 8))
    downstream_sort_threads = max(2, min(6, logical // max(1, downstream_workers) - downstream_threads))
    prefetch_slots = max(1, min(4, audit_workers))
    fasterq_slots = max(1, min(4, audit_workers))
    fasterq_threads = max(2, min(6, logical // max(1, fasterq_slots * 4)))
    reference_threads = max(4, min(48, physical - 8 if physical > 16 else physical))
    repair_mem = 16 if mem >= 256 else 8 if mem >= 128 else 4
    sort_mem = "3G" if mem >= 256 else "2G" if mem >= 128 else "1G"

    cfg = {
        "resources": {
            "profile_name": "autotuned_runtime",
            "detected_physical_cores": physical,
            "detected_logical_cpus": logical,
            "detected_memory_gb": round(mem, 1),
            "audit_workers": audit_workers,
            "audit_threads_per_run": audit_threads,
            "prefetch_slots": prefetch_slots,
            "fasterq_slots": fasterq_slots,
            "fasterq_threads": fasterq_threads,
            "downstream_workers": downstream_workers,
            "downstream_threads_per_run": downstream_threads,
            "downstream_sort_threads": downstream_sort_threads,
            "downstream_sort_mem": sort_mem,
            "bcftools_threads": min(4, downstream_sort_threads),
            "reference_threads": reference_threads,
            "repair_memory_gb": repair_mem,
            "workers": audit_workers,
            "threads_per_run": audit_threads,
            "sort_threads": downstream_sort_threads,
            "sort_mem": sort_mem,
            "download_threads": fasterq_threads,
        }
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    report = Path(args.report)
    report.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        ("platform", platform.platform()),
        ("physical_cores", physical),
        ("logical_cpus_available", logical),
        ("memory_gb", f"{mem:.1f}"),
        ("audit_workers", audit_workers),
        ("audit_threads_per_run", audit_threads),
        ("prefetch_slots", prefetch_slots),
        ("fasterq_slots", fasterq_slots),
        ("fasterq_threads", fasterq_threads),
        ("downstream_workers", downstream_workers),
        ("downstream_threads_per_run", downstream_threads),
        ("downstream_sort_threads", downstream_sort_threads),
        ("reference_threads", reference_threads),
        ("runtime_config", str(out)),
    ]
    report.write_text("item\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in rows) + "\n", encoding="utf-8")
    print(json.dumps(dict(rows), ensure_ascii=False, indent=2))
    print(f"[OK] runtime resource override written: {out}")


if __name__ == "__main__":
    main()
