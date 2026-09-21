#!/usr/bin/env python3
"""Single entry point for the organelle-audit workflow."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from common import load_config, species_id, workflow_env

ROOT = Path(__file__).resolve().parents[1]
SPECIES_CONFIGS = {
    "silene": "config/species/silene.yaml",
    "chicken": "config/species/chicken.yaml",
}


def configs(target: str) -> list[str]:
    return list(SPECIES_CONFIGS.values()) if target == "both" else [SPECIES_CONFIGS[target]]


def call(script: str, *args: str) -> None:
    cmd = [sys.executable, f"workflow/{script}", *args]
    print("[STEP] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=workflow_env())


def timed(label: str, fn) -> None:
    start = time.monotonic()
    print(f"\n=== {label} ===", flush=True)
    fn()
    print(f"[DONE] {label} ({time.monotonic() - start:.1f} s)", flush=True)


def prepare(cfg_path: str, mode: str, force_index: bool = False) -> None:
    call("doctor.py", "--config", cfg_path)
    call("manifest.py", "--config", cfg_path, "--mode", mode)
    idx_args = ["--config", cfg_path]
    if force_index:
        idx_args.append("--force")
    call("local_inputs.py", *idx_args)
    call("references.py", "--config", cfg_path)


def audit(cfg_path: str) -> None:
    call("scheduler.py", "--config", cfg_path, "--stage", "audit")
    call("audit_summary.py", "--config", cfg_path)


def downstream(cfg_path: str) -> None:
    call("subset.py", "--config", cfg_path)
    call("scheduler.py", "--config", cfg_path, "--stage", "downstream")
    call("downstream_summary.py", "--config", cfg_path)


def report_all(cfg_paths: list[str]) -> None:
    for cfg in cfg_paths:
        call("threshold_sensitivity.py", "--config", cfg)
        call("plot_species.py", "--config", cfg)
    call("plot_combined.py")
    call("plot_workflow.py")
    call("report_materials.py")
    call("validate_outputs.py")
    call("performance.py")


def show_plan(cfg_paths: list[str], mode: str) -> None:
    for path in cfg_paths:
        cfg = load_config(path)
        sid = species_id(cfg)
        print(f"\n[{sid}]")
        print(f"  scientific_name: {cfg['species']['scientific_name']}")
        print(f"  reference_mode: {cfg['species']['mode']}")
        print(f"  pinned_manifest: {cfg['species']['exact_manifest']}")
        print(f"  expected_runs: {cfg['species'].get('expected_evaluable_runs')}")
        print(f"  mode: {mode}")
        print("  read_roots:")
        for p in cfg.get("local_inputs", {}).get("read_roots", []):
            print(f"    - {p}")
        print("  references:")
        refs = cfg.get("local_inputs", {}).get("references", {})
        print(f"    nuclear: {refs.get('nuclear_fasta', 'not configured')}")
        for name, spec in refs.get("organelles", {}).items():
            print(f"    {name}: {spec.get('fasta') or 'extract ' + str(spec.get('extract_accession_from_nuclear'))}")
        r = cfg.get("resources", {})
        print(f"  audit_parallelism: {r.get('audit_workers')} workers x {r.get('audit_threads_per_run')} threads")
        print(f"  downstream_parallelism: {r.get('downstream_workers')} workers x {r.get('downstream_threads_per_run')} threads")
        print(f"  temp_root: {r.get('tmp_root')}")
        print("  network_downloads: disabled")


def clean(level: str, target: str) -> None:
    selected = ["silene", "chicken"] if target == "both" else [target]
    paths: list[Path] = []
    if level in {"tmp", "intermediate", "all"}:
        paths.append(Path("tmp"))
    if level in {"intermediate", "all"}:
        paths.extend(Path("work") / s / p for s in selected for p in ["runs", "downstream"])
        paths.append(Path("work/.resource_locks"))
    if level == "all":
        paths.extend(Path("work") / s for s in selected)
        paths.extend(Path("results") / s for s in selected)
        paths.extend(Path("logs") / s for s in selected)
        paths.append(Path("results/combined"))
    for p in sorted(set(paths), key=lambda x: len(str(x)), reverse=True):
        if p.exists():
            print(f"[CLEAN] {p}")
            shutil.rmtree(p)


def main() -> None:
    os.chdir(ROOT)
    os.environ.setdefault("PYTHONPATH", str(ROOT / "workflow"))
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=[
        "plan", "autotune", "doctor", "local-check", "references", "pilot", "full",
        "audit", "downstream", "report", "status", "readiness", "performance", "clean",
    ])
    ap.add_argument("target", nargs="?", default="both", choices=["both", "silene", "chicken"])
    ap.add_argument("--clean-level", default="tmp", choices=["tmp", "intermediate", "all"])
    args = ap.parse_args()
    cfgs = configs(args.target)

    if args.command == "plan":
        show_plan(cfgs, "full")
        return
    if args.command == "autotune":
        call("autotune.py")
        return
    if args.command == "status":
        call("status.py", *cfgs)
        return
    if args.command == "performance":
        call("performance.py")
        return
    if args.command == "readiness":
        call("validate_outputs.py", "--strict")
        return
    if args.command == "clean":
        clean(args.clean_level, args.target)
        return

    if args.command == "doctor":
        for cfg in cfgs:
            timed(f"doctor {cfg}", lambda c=cfg: call("doctor.py", "--config", c))
        return
    if args.command == "local-check":
        for cfg in cfgs:
            def check_local(c=cfg):
                call("doctor.py", "--config", c)
                call("manifest.py", "--config", c, "--mode", "full")
                call("local_inputs.py", "--config", c, "--force")
            timed(f"local input check {cfg}", check_local)
        return
    if args.command == "references":
        for cfg in cfgs:
            timed(f"reference preparation {cfg}", lambda c=cfg: prepare(c, "pilot"))
        return
    if args.command == "audit":
        for cfg in cfgs:
            timed(f"audit {cfg}", lambda c=cfg: (prepare(c, "full"), audit(c)))
        return
    if args.command == "downstream":
        for cfg in cfgs:
            timed(f"downstream {cfg}", lambda c=cfg: downstream(c))
        report_all(cfgs)
        return
    if args.command == "report":
        report_all(cfgs)
        return

    mode = args.command
    for cfg in cfgs:
        timed(f"prepare {cfg}", lambda c=cfg: prepare(c, mode))
        timed(f"audit {cfg}", lambda c=cfg: audit(c))
        timed(f"downstream {cfg}", lambda c=cfg: downstream(c))
    timed("combined reports", lambda: report_all(cfgs))


if __name__ == "__main__":
    main()
