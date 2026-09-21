#!/usr/bin/env python3
from __future__ import annotations

import sys
import argparse
import json
import shlex
import shutil
import subprocess
import time
import traceback
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from common import (
    atomic_json,
    ensure_dirs,
    file_slot,
    load_config,
    read_qname_set,
    resolve_local_read_inputs,
    local_source_policy,
    resource_value,
    require_free_space,
    scratch_dir,
    run_cmd,
    run_cmd_retry,
    species_id,
    workdir,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def shell(cmd: str, log: Path) -> None:
    run_cmd(["bash", "-o", "pipefail", "-c", cmd], log=log)


def find_sra(sra_root: Path, run: str) -> Path:
    candidates = [sra_root / run / f"{run}.sra", sra_root / f"{run}.sra"]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p
    found = list(sra_root.rglob(f"{run}.sra"))
    if found:
        return found[0]
    raise FileNotFoundError(f"SRA file not found after prefetch for {run}")


def qname_path(run_dir: Path, tag: str, q: int) -> Path:
    return run_dir / "sets" / f"{tag}.q{q}.qnames.txt.gz"


def classify(run: str, spots: int, compartments: list[str], mapqs: list[int],
             run_dir: Path, out_path: Path, nuclear_available: bool) -> None:
    rows = []
    for q in mapqs:
        nuclear = read_qname_set(qname_path(run_dir, "nuclear", q)) if nuclear_available else set()
        comp_sets = {c: read_qname_set(qname_path(run_dir, c, q)) for c in compartments}
        organelle = set().union(*comp_sets.values()) if comp_sets else set()
        strict = organelle - nuclear if nuclear_available else set()
        multi = organelle & nuclear if nuclear_available else set()
        row = {
            "Run": run,
            "MAPQ": q,
            "spots": spots,
            "reference_mode": "nuclear_plus_organelle" if nuclear_available else "organelle_only",
            "nuclear_qnames": len(nuclear) if nuclear_available else float("nan"),
            "organelle_inclusive_qnames": len(organelle),
            "organelle_strict_qnames": len(strict) if nuclear_available else float("nan"),
            "multi_reference_qnames": len(multi) if nuclear_available else float("nan"),
            "frac_inclusive_total": len(organelle) / spots if spots else float("nan"),
            "frac_strict_total": len(strict) / spots if (spots and nuclear_available) else float("nan"),
            "frac_multi_total": len(multi) / spots if (spots and nuclear_available) else float("nan"),
        }
        for c, s in comp_sets.items():
            row[f"{c}_qnames"] = len(s)
            row[f"frac_{c}_total"] = len(s) / spots if spots else float("nan")
        rows.append(row)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def run_streaming_mapping(ref: Path, fastqs: list[str], tag: str, run: str,
                          threads: int, qvals: list[int], run_dir: Path,
                          cfg: dict, log: Path) -> None:
    rg = f"@RG\tID:{run}.{tag}\tSM:{run}\tPL:ILLUMINA"
    bwa = ["bwa-mem2", "mem", "-t", str(threads), "-R", rg, str(ref), *fastqs]
    parser = [
        sys.executable, "workflow/sam_qnames.py",
        "--output-dir", str(run_dir / "sets"),
        "--tag", tag,
        "--thresholds", *map(str, qvals),
        "--summary", str(run_dir / "sets" / f"{tag}.mapping_summary.json"),
        "--compression-level", str(int(resource_value(cfg, "audit_qname_compression_level", None, 1))),
    ]
    shell(shlex.join(bwa) + " | " + shlex.join(parser), log)


def run_legacy_bam_mapping(ref: Path, fastqs: list[str], tag: str, run: str,
                           threads: int, sort_threads: int, sort_mem: str,
                           qvals: list[int], run_dir: Path, tmp_dir: Path, log: Path) -> None:
    bam = run_dir / "bam" / f"{tag}.bam"
    rg = f"@RG\tID:{run}.{tag}\tSM:{run}\tPL:ILLUMINA"
    bwa = shlex.join(["bwa-mem2", "mem", "-t", str(threads), "-R", rg, str(ref), *fastqs])
    sort = shlex.join([
        "samtools", "sort", "-@", str(sort_threads), "-m", sort_mem,
        "-l", "1", "-T", str(tmp_dir / tag), "-o", str(bam), "-",
    ])
    shell(bwa + " | " + sort, log)
    run_cmd(["samtools", "index", "-@", str(sort_threads), str(bam)], log=log)
    sam = shlex.join(["samtools", "view", "-@", str(sort_threads), "-F", "2308", str(bam)])
    parser = shlex.join([
        sys.executable, "workflow/sam_qnames.py", "--output-dir", str(run_dir / "sets"),
        "--tag", tag, "--thresholds", *map(str, qvals),
        "--summary", str(run_dir / "sets" / f"{tag}.mapping_summary.json"),
        "--compression-level", "1",
    ])
    shell(sam + " | " + parser, log)
    bam.unlink(missing_ok=True)
    Path(str(bam) + ".bai").unlink(missing_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    ensure_dirs(cfg)
    sid = species_id(cfg)
    manifest = pd.read_csv(workdir(cfg) / "manifest/manifest.active.tsv", sep="\t")
    sub = manifest.loc[manifest["Run"] == args.run]
    if sub.empty:
        raise SystemExit(f"[ERROR] run not in active pinned manifest: {args.run}")
    row = sub.iloc[0]
    layout = str(row["LibraryLayout"]).upper()
    strategy = str(row.get("LibraryStrategy", "")).upper()
    spots = int(row["spots"])
    run_dir = workdir(cfg) / "runs" / args.run
    tmp_dir = scratch_dir(cfg, "audit", args.run)
    require_free_space(tmp_dir, float(resource_value(cfg, "minimum_free_tmp_gb", None, 20)))
    log = run_dir / "audit.log"
    status_path = run_dir / "audit.status.json"
    done = run_dir / ".audit.done"
    if done.exists() and not args.force:
        print(f"[SKIP] audit complete: {args.run}")
        return
    for d in ["sra", "raw", "paired", "singletons", "trimmed", "bam", "sets", "tmp", "fastp"]:
        (run_dir / d).mkdir(parents=True, exist_ok=True)

    started = time.monotonic()
    status = {
        "Run": args.run, "species": sid, "layout": layout, "spots": spots,
        "stage": "start", "status": "running", "started_utc": utc_now(),
    }
    atomic_json(status_path, status)

    threads = int(resource_value(cfg, "audit_threads_per_run", "threads_per_run", 8))
    sort_threads = int(resource_value(cfg, "downstream_sort_threads", "sort_threads", 4))
    sort_mem = str(resource_value(cfg, "downstream_sort_mem", "sort_mem", "1G"))
    fasterq_threads = int(resource_value(cfg, "fasterq_threads", "download_threads", 4))
    prefetch_slots = int(resource_value(cfg, "prefetch_slots", None, 2))
    fasterq_slots = int(resource_value(cfg, "fasterq_slots", None, 2))
    repair_mem = int(resource_value(cfg, "repair_memory_gb", None, 8))
    qvals = list(map(int, cfg["workflow"]["mapq_thresholds"]))
    primary = int(cfg["workflow"]["primary_mapq"])
    compartments = list(cfg["references"]["organelles"].keys())
    nuclear_available = bool(cfg["references"].get("nuclear_assembly_accession"))

    heavy_strategies = {str(x).upper() for x in cfg.get("resources", {}).get("heavy_run_strategies", ["WGS"])}
    heavy_threshold = int(resource_value(cfg, "heavy_run_spots_threshold", None, 20000000))
    heavy_slots = int(resource_value(cfg, "heavy_run_slots", None, 2))
    is_heavy = strategy in heavy_strategies or spots >= heavy_threshold
    stack = ExitStack()
    if is_heavy:
        stack.enter_context(file_slot(cfg, "heavy_run", heavy_slots))
        status["heavy_run_throttled"] = True
        atomic_json(status_path, status)

    try:
        local_assignment = resolve_local_read_inputs(cfg, args.run, layout)
        read_policy = local_source_policy(cfg, "reads")
        source_r1: Path | None = None
        source_r2: Path | None = None
        source_single: Path | None = None
        if local_assignment is not None:
            status["stage"] = "local_fastq_reuse"
            status["read_source"] = "local_fastq"
            status["local_r1"] = local_assignment.get("r1", "")
            status["local_r2"] = local_assignment.get("r2", "")
            status["local_single"] = local_assignment.get("single", "")
            atomic_json(status_path, status)
            if layout == "PAIRED":
                source_r1 = Path(str(local_assignment["r1"]))
                source_r2 = Path(str(local_assignment["r2"]))
            else:
                source_single = Path(str(local_assignment["single"]))
            print(f"[INFO] reusing local FASTQ for {args.run}")
        else:
            allow_download = bool(cfg.get("network", {}).get("allow_read_download", True))
            if read_policy == "local_required" or not allow_download:
                raise FileNotFoundError(
                    f"validated local FASTQ assignment missing for {args.run}; "
                    f"run workflow/local_inputs.py --config {args.config} and inspect "
                    f"{workdir(cfg) / 'local_inputs' / 'unresolved_runs.tsv'}"
                )
            status["stage"] = "prefetch"
            status["read_source"] = "SRA_download"
            atomic_json(status_path, status)
            with file_slot(cfg, "prefetch", prefetch_slots):
                prefetch_cmd = [
                    "prefetch", args.run, "-O", str(run_dir / "sra"),
                    "--max-size", str(resource_value(cfg, "prefetch_max_size", None, "200G")),
                ]
                run_cmd_retry(prefetch_cmd, attempts=int(cfg["network"]["retries"]),
                              sleep_seconds=int(cfg["network"]["retry_sleep_seconds"]), log=log)

            sra = find_sra(run_dir / "sra", args.run)
            status["stage"] = "fasterq_dump"
            atomic_json(status_path, status)
            with file_slot(cfg, "fasterq", fasterq_slots):
                run_cmd_retry([
                    "fasterq-dump", str(sra), "--split-3", "-e", str(fasterq_threads),
                    "-t", str(tmp_dir), "-O", str(run_dir / "raw"), "--force",
                ], attempts=2, sleep_seconds=10, log=log)
            if layout == "PAIRED":
                source_r1 = run_dir / "raw" / f"{args.run}_1.fastq"
                source_r2 = run_dir / "raw" / f"{args.run}_2.fastq"
            else:
                candidates = [run_dir / "raw" / f"{args.run}.fastq", run_dir / "raw" / f"{args.run}_1.fastq"]
                source_single = next((p for p in candidates if p.exists() and p.stat().st_size > 0), None)

        status["stage"] = "repair_preprocess"
        atomic_json(status_path, status)
        pp = cfg["preprocess"]
        fastp_json = run_dir / "fastp" / f"{args.run}.json"
        fastp_html = run_dir / "fastp" / f"{args.run}.html"
        compression = str(int(pp.get("compression_level", 1)))
        if layout == "PAIRED":
            r1 = source_r1
            r2 = source_r2
            if r1 is None or r2 is None or not (r1.exists() and r2.exists() and r1.stat().st_size > 0 and r2.stat().st_size > 0):
                raise FileNotFoundError(f"paired FASTQ inputs missing: {r1}, {r2}")
            p1 = run_dir / "paired" / f"{args.run}_1.paired.fastq.gz"
            p2 = run_dir / "paired" / f"{args.run}_2.paired.fastq.gz"
            single = run_dir / "singletons" / f"{args.run}.singletons.fastq.gz"
            run_cmd([
                "repair.sh", f"-Xmx{repair_mem}g", f"in1={r1}", f"in2={r2}",
                f"out1={p1}", f"out2={p2}", f"outs={single}",
                "repair=t", "overwrite=t", "zl=1",
            ], log=log)
            t1 = run_dir / "trimmed" / f"{args.run}_1.trim.fastq.gz"
            t2 = run_dir / "trimmed" / f"{args.run}_2.trim.fastq.gz"
            cmd = [
                "fastp", "--thread", str(threads), "--in1", str(p1), "--in2", str(p2),
                "--out1", str(t1), "--out2", str(t2), "--compression", compression,
                "--qualified_quality_phred", str(pp["qualified_quality_phred"]),
                "--length_required", str(pp["length_required"]),
                "--unqualified_percent_limit", str(pp["unqualified_percent_limit"]),
                "--n_base_limit", str(pp["n_base_limit"]),
                "--cut_right_window_size", str(pp["cut_right_window_size"]),
                "--cut_right_mean_quality", str(pp["cut_right_mean_quality"]),
                "--json", str(fastp_json), "--html", str(fastp_html),
            ]
            if pp.get("detect_adapter_for_pe", True):
                cmd.append("--detect_adapter_for_pe")
            if pp.get("trim_poly_g", True):
                cmd.append("--trim_poly_g")
            if pp.get("cut_right", True):
                cmd.append("--cut_right")
            run_cmd(cmd, log=log)
            run_cmd([sys.executable, "workflow/pair_check.py", "--r1", str(t1), "--r2", str(t2)], log=log)
            fastqs = [str(t1), str(t2)]
        else:
            raw = source_single
            if raw is None or not raw.exists() or raw.stat().st_size <= 0:
                raise FileNotFoundError(f"single-end FASTQ input missing for {args.run}: {raw}")
            t1 = run_dir / "trimmed" / f"{args.run}.trim.fastq.gz"
            cmd = [
                "fastp", "--thread", str(threads), "--in1", str(raw), "--out1", str(t1),
                "--compression", compression,
                "--qualified_quality_phred", str(pp["qualified_quality_phred"]),
                "--length_required", str(pp["length_required"]),
                "--unqualified_percent_limit", str(pp["unqualified_percent_limit"]),
                "--n_base_limit", str(pp["n_base_limit"]),
                "--cut_right_window_size", str(pp["cut_right_window_size"]),
                "--cut_right_mean_quality", str(pp["cut_right_mean_quality"]),
                "--json", str(fastp_json), "--html", str(fastp_html),
            ]
            if pp.get("trim_poly_g", True):
                cmd.append("--trim_poly_g")
            if pp.get("cut_right", True):
                cmd.append("--cut_right")
            run_cmd(cmd, log=log)
            fastqs = [str(t1)]

        status["stage"] = "mapping_stream"
        atomic_json(status_path, status)
        prep = workdir(cfg) / "references" / "prepared"
        refs = {c: prep / f"{c}.fa" for c in compartments}
        if nuclear_available:
            refs["nuclear"] = prep / "nuclear.masked.main.fa"
        streaming = bool(cfg["workflow"].get("audit_streaming_qnames", True))
        for tag, ref in refs.items():
            if not ref.exists():
                raise FileNotFoundError(f"prepared reference missing: {ref}")
            if streaming:
                run_streaming_mapping(ref, fastqs, tag, args.run, threads, qvals, run_dir, cfg, log)
            else:
                run_legacy_bam_mapping(ref, fastqs, tag, args.run, threads, sort_threads,
                                       sort_mem, qvals, run_dir, tmp_dir, log)

        audit_path = run_dir / "audit_rows.tsv"
        classify(args.run, spots, compartments, qvals, run_dir, audit_path, nuclear_available)
        if cfg["workflow"].get("retain_primary_mapq_sets", True):
            for p in (run_dir / "sets").glob("*.q*.qnames.txt.gz"):
                if f".q{primary}." not in p.name:
                    p.unlink()
        else:
            shutil.rmtree(run_dir / "sets", ignore_errors=True)
        if cfg["workflow"].get("cleanup_raw_sra_after_run", True):
            shutil.rmtree(run_dir / "sra", ignore_errors=True)
        if cfg["workflow"].get("cleanup_raw_fastq_after_run", True):
            for d in ["raw", "paired", "singletons"]:
                shutil.rmtree(run_dir / d, ignore_errors=True)
        if cfg["workflow"].get("cleanup_trimmed_after_audit", True):
            shutil.rmtree(run_dir / "trimmed", ignore_errors=True)
        shutil.rmtree(tmp_dir, ignore_errors=True)

        status.update({
            "stage": "complete", "status": "ready", "audit_table": str(audit_path),
            "completed_utc": utc_now(), "wall_seconds": round(time.monotonic() - started, 3),
            "audit_streaming_qnames": streaming,
        })
        atomic_json(status_path, status)
        done.write_text("ok\n", encoding="utf-8")
        print(f"[OK] audit run complete: {args.run} wall={status['wall_seconds']}s")
    except Exception as exc:
        status.update({
            "status": "failed", "error": str(exc), "traceback": traceback.format_exc(),
            "failed_utc": utc_now(), "wall_seconds": round(time.monotonic() - started, 3),
        })
        atomic_json(status_path, status)
        print(f"[ERROR] audit run failed {args.run}: {exc}", flush=True)
        raise
    finally:
        stack.close()


if __name__ == "__main__":
    main()
