#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

from common import command_exists, ensure_dirs, load_config, resource_value, sha256

CORE_REQUIRED = [
    "python", "fastp", "repair.sh", "seqkit", "minimap2", "bedtools",
    "bwa-mem2", "samtools", "bcftools"
]

NETWORK_READ_TOOLS = ["prefetch", "fasterq-dump"]
NETWORK_REFERENCE_TOOLS = ["datasets", "efetch", "unzip"]


VERSION_COMMANDS = {
    "python": ["python", "--version"],
    "prefetch": ["prefetch", "--version"],
    "fasterq-dump": ["fasterq-dump", "--version"],
    "fastp": ["fastp", "--version"],
    "seqkit": ["seqkit", "version"],
    "datasets": ["datasets", "version"],
    "minimap2": ["minimap2", "--version"],
    "bedtools": ["bedtools", "--version"],
    "bwa-mem2": ["bwa-mem2", "version"],
    "samtools": ["samtools", "--version"],
    "bcftools": ["bcftools", "--version"],
}


def version(cmd):
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
        return p.stdout.strip().splitlines()[0] if p.stdout.strip() else "available"
    except Exception as exc:
        return f"available; version query failed: {exc}"



def conda_package_version(package: str) -> str:
    """Return a Conda package version from the active prefix without invoking network tools."""
    meta = Path(sys.prefix) / "conda-meta"
    if meta.is_dir():
        hits = sorted(meta.glob(f"{package}-*.json"))
        for path in hits:
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                if str(data.get("name", "")).lower() == package.lower() and data.get("version"):
                    return str(data["version"])
            except Exception:
                continue
    return "available; exact package version not resolved from active Conda prefix"

def memory_gb() -> float:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) / 1024 / 1024
    except Exception:
        return 0.0
    return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    required = list(CORE_REQUIRED)
    if bool(cfg.get("network", {}).get("allow_read_download", False)):
        required.extend(NETWORK_READ_TOOLS)
    if bool(cfg.get("network", {}).get("allow_reference_download", False)):
        required.extend(NETWORK_REFERENCE_TOOLS)
    required = list(dict.fromkeys(required))
    missing = [x for x in required if not command_exists(x)]
    if missing:
        raise SystemExit("[ERROR] missing required commands for the configured source policy: " + ", ".join(missing))

    r = cfg.get("resources", {})
    logical = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    mem = memory_gb()
    audit_workers = int(resource_value(cfg, "audit_workers", "workers", 1))
    audit_threads = int(resource_value(cfg, "audit_threads_per_run", "threads_per_run", 1))
    downstream_workers = int(resource_value(cfg, "downstream_workers", None, 1))
    downstream_threads = int(resource_value(cfg, "downstream_threads_per_run", "threads_per_run", 1))
    downstream_sort_threads = int(resource_value(cfg, "downstream_sort_threads", "sort_threads", 1))

    warnings: list[str] = []
    if audit_workers * audit_threads > logical:
        warnings.append(f"audit_workers × audit_threads_per_run = {audit_workers * audit_threads} exceeds {logical} available logical CPUs")
    if downstream_workers * (downstream_threads + downstream_sort_threads) > logical:
        warnings.append(
            "downstream_workers × (mapping + sort threads) = "
            f"{downstream_workers * (downstream_threads + downstream_sort_threads)} exceeds {logical} available logical CPUs"
        )
    expected_mem = float(r.get("expected_memory_gb", 0) or 0)
    if expected_mem and mem and mem < expected_mem * 0.85:
        warnings.append(f"detected memory {mem:.1f} GB is lower than profile expectation {expected_mem:.0f} GB")

    out = Path("results") / cfg["species"]["id"] / "provenance" / "software_versions.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        ("workflow_version", cfg["workflow"]["version"]),
        ("platform", platform.platform()),
        ("resource_profile_name", r.get("profile_name", "unspecified")),
        ("available_logical_cpus", logical),
        ("detected_memory_gb", f"{mem:.1f}"),
        ("audit_workers", audit_workers),
        ("audit_threads_per_run", audit_threads),
        ("downstream_workers", downstream_workers),
        ("downstream_threads_per_run", downstream_threads),
        ("downstream_sort_threads", downstream_sort_threads),
        ("prefetch_slots", r.get("prefetch_slots", "")),
        ("fasterq_slots", r.get("fasterq_slots", "")),
        ("fasterq_threads", r.get("fasterq_threads", "")),
        ("heavy_run_slots", r.get("heavy_run_slots", "")),
        ("reference_threads", r.get("reference_threads", "")),
        ("repair_memory_gb", r.get("repair_memory_gb", "")),
        ("tmp_root", r.get("tmp_root", "tmp")),
        ("cpu_affinity", r.get("cpu_affinity", False)),
        ("taskset_available", shutil.which("taskset") is not None),
    ]
    rows.extend([
        ("read_source_policy", cfg.get("local_inputs", {}).get("read_source_policy", "prefer_local")),
        ("reference_source_policy", cfg.get("local_inputs", {}).get("reference_source_policy", "prefer_local")),
        ("allow_read_download", cfg.get("network", {}).get("allow_read_download", False)),
        ("allow_reference_download", cfg.get("network", {}).get("allow_reference_download", False)),
    ])
    for tool in required:
        if tool == "repair.sh":
            rows.append(("bbmap", conda_package_version("bbmap")))
            rows.append(("repair.sh_path", shutil.which("repair.sh") or "not_found"))
        else:
            rows.append((tool, version(VERSION_COMMANDS.get(tool, [tool, "--version"]))))
    rows.extend([
        ("common_config_sha256", sha256(cfg["_common_config"])),
        ("resource_profile_sha256", sha256(cfg["_resource_profile"]) if cfg.get("_resource_profile") and Path(cfg["_resource_profile"]).exists() else "not_used"),
        ("runtime_config_sha256", sha256(cfg["_runtime_config"]) if cfg.get("_runtime_config") and Path(cfg["_runtime_config"]).exists() else "not_used"),
        ("species_config_sha256", sha256(cfg["_species_config"])),
        ("exact_manifest_sha256", sha256(cfg["species"]["exact_manifest"])),
    ])
    if warnings:
        rows.extend(("resource_warning", w) for w in warnings)
    out.write_text("item\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in rows) + "\n", encoding="utf-8")
    for w in warnings:
        print(f"[WARN] {w}")
    print(f"[OK] environment, hardware profile, and provenance checks passed for {cfg['species']['id']}")


if __name__ == "__main__":
    main()
