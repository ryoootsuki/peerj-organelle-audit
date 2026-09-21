#!/usr/bin/env python3
"""Dynamic, resource-aware run scheduler.

One scheduler serves both archive audit and selected downstream analysis. Runs are
queued largest-first using SRA spots as a stable cost proxy. Each worker owns a
fixed CPU-affinity slot and pulls the next pending run when it becomes free,
which avoids the long-tail idle time of static partitions.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import queue
import subprocess
import sys
import threading
from pathlib import Path

import pandas as pd

from common import load_config, resource_value, taskset_prefix, workdir, workflow_env

STAGES = {
    "audit": {
        "worker": "workflow/audit_worker.py",
        "done_glob": "runs/*/.audit.done",
        "run_file": None,
        "workers_key": "audit_workers",
        "failure_file": "audit_failed_runs.txt",
    },
    "downstream": {
        "worker": "workflow/downstream_worker.py",
        "done_glob": "downstream/*/.done",
        "run_file": "manifest/downstream_runs.txt",
        "workers_key": "downstream_workers",
        "failure_file": "downstream_failed_runs.txt",
    },
}


def pending_table(cfg: dict, stage: str, force: bool) -> pd.DataFrame:
    wd = workdir(cfg)
    manifest = pd.read_csv(wd / "manifest/manifest.evaluable.tsv", sep="\t", dtype={"Run": str})
    spec = STAGES[stage]
    if stage == "audit":
        active = pd.read_csv(wd / "manifest/manifest.active.tsv", sep="\t", dtype={"Run": str})
        wanted = set(active["Run"].astype(str))
    else:
        p = wd / str(spec["run_file"])
        if not p.exists():
            raise SystemExit(f"[ERROR] selected-run list does not exist: {p}; run subset selection first")
        wanted = {x.strip() for x in p.read_text().splitlines() if x.strip()}
    out = manifest.loc[manifest["Run"].astype(str).isin(wanted)].copy()
    if not force:
        done = {p.parent.name for p in wd.glob(str(spec["done_glob"]))}
        out = out.loc[~out["Run"].astype(str).isin(done)].copy()
    out["spots"] = pd.to_numeric(out.get("spots", 1), errors="coerce").fillna(1)
    return out.sort_values(["spots", "Run"], ascending=[False, True]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--stage", required=True, choices=sorted(STAGES))
    ap.add_argument("--workers", type=int)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    spec = STAGES[args.stage]
    table = pending_table(cfg, args.stage, args.force)
    if table.empty:
        print(f"[OK] all {args.stage} runs are already complete")
        return

    n_workers = args.workers or int(resource_value(cfg, str(spec["workers_key"]), None, 1))
    n_workers = max(1, min(n_workers, len(table)))
    jobs: queue.PriorityQueue[tuple[float, str]] = queue.PriorityQueue()
    for row in table.itertuples(index=False):
        jobs.put((-float(getattr(row, "spots", 1) or 1), str(row.Run)))

    print(
        f"[INFO] {args.stage} scheduler: policy=dynamic-largest-first, "
        f"workers={n_workers}, pending_runs={len(table)}, weight=spots",
        flush=True,
    )
    failures: list[str] = []
    lock = threading.Lock()

    def worker_loop(slot: int) -> None:
        prefix = taskset_prefix(cfg, slot, n_workers)
        while True:
            try:
                _, run = jobs.get_nowait()
            except queue.Empty:
                return
            cmd = prefix + [
                sys.executable,
                str(spec["worker"]),
                "--config", args.config,
                "--run", run,
            ]
            if args.force:
                cmd.append("--force")
            try:
                proc = subprocess.run(cmd, env=workflow_env())
                rc = proc.returncode
            except Exception as exc:
                rc = 99
                print(f"[ERROR] scheduler could not launch stage={args.stage} run={run}: {exc}", flush=True)
            print(f"[STATUS] stage={args.stage} slot={slot} run={run} exit={rc}", flush=True)
            if rc:
                with lock:
                    failures.append(run)
            jobs.task_done()

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [pool.submit(worker_loop, slot) for slot in range(n_workers)]
        for fut in concurrent.futures.as_completed(futures):
            fut.result()

    fail_path = workdir(cfg) / "manifest" / str(spec["failure_file"])
    if failures:
        fail_path.write_text("\n".join(sorted(set(failures))) + "\n", encoding="utf-8")
        msg = f"{len(set(failures))} {args.stage} run(s) failed; list={fail_path}"
        if args.stage == "downstream" or cfg.get("workflow", {}).get("missing_read_policy") == "fail":
            raise SystemExit("[ERROR] " + msg)
        print("[WARN] " + msg)
    else:
        fail_path.unlink(missing_ok=True)
        print(f"[OK] {args.stage} completed for {len(table)} pending run(s)")


if __name__ == "__main__":
    main()
