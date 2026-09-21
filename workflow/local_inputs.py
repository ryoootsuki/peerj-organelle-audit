#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from common import load_config, normalize_qname, species_id, workdir

RUN_SEARCH_RE = re.compile(r"(?:SRR|ERR|DRR)[0-9]{6,}", re.IGNORECASE)


@dataclass(frozen=True)
class Candidate:
    run: str
    role: str
    path: str
    root_index: int
    pattern_score: int
    pair_key: str
    size_bytes: int
    mtime_ns: int


def strip_fastq_extension(name: str, accepted: list[str]) -> tuple[str, str] | None:
    lower = name.lower()
    for ext in sorted(accepted, key=len, reverse=True):
        if lower.endswith(ext.lower()):
            return name[: -len(ext)], ext
    return None


def classify_stem(run: str, stem: str) -> tuple[str, int, str] | None:
    """Return role, score and pair key for common SRA FASTQ naming conventions."""
    if stem == run:
        return "single", 110, run

    exact_patterns = [
        (rf"^{re.escape(run)}_([12])$", 110),
        (rf"^{re.escape(run)}\.([12])$", 108),
        (rf"^{re.escape(run)}-([12])$", 106),
        (rf"^{re.escape(run)}_R([12])$", 104),
        (rf"^{re.escape(run)}\.R([12])$", 102),
        (rf"^{re.escape(run)}-R([12])$", 100),
    ]
    for pat, score in exact_patterns:
        m = re.match(pat, stem, flags=re.IGNORECASE)
        if m:
            mate = m.group(1)
            key = re.sub(r"([_.-])R?[12]$", r"\1<M>", stem, flags=re.IGNORECASE)
            return ("r1" if mate == "1" else "r2"), score, key

    # Tolerate a limited suffix after the mate token, but score below exact raw names.
    m = re.match(
        rf"^{re.escape(run)}([_.-])R?([12])([_.-].+)$",
        stem,
        flags=re.IGNORECASE,
    )
    if m:
        mate = m.group(2)
        key = f"{run}{m.group(1)}<M>{m.group(3)}"
        return ("r1" if mate == "1" else "r2"), 80, key
    return None


def iter_fastq_files(root: Path, recursive: bool, followlinks: bool) -> Iterable[Path]:
    if root.is_file():
        yield root
        return
    if not root.is_dir():
        return
    if recursive:
        for dirpath, _, filenames in os.walk(root, followlinks=followlinks):
            base = Path(dirpath)
            for fn in filenames:
                yield base / fn
    else:
        with os.scandir(root) as it:
            for ent in it:
                if ent.is_file(follow_symlinks=followlinks):
                    yield Path(ent.path)


def open_fastq(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.name.lower().endswith(".gz") else path.open("rt", encoding="utf-8", errors="replace")


def first_record(path: Path) -> tuple[str, str]:
    with open_fastq(path) as fh:
        h = fh.readline().rstrip("\r\n")
        seq = fh.readline().rstrip("\r\n")
        plus = fh.readline().rstrip("\r\n")
        qual = fh.readline().rstrip("\r\n")
    if not h.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
        raise ValueError(f"invalid first FASTQ record: {path}")
    return normalize_qname(h[1:]), seq


def choose_pair(r1s: list[Candidate], r2s: list[Candidate]) -> tuple[Candidate, Candidate, int]:
    pairs: list[tuple[tuple[int, ...], Candidate, Candidate]] = []
    for a in r1s:
        for b in r2s:
            same_key = int(a.pair_key == b.pair_key)
            same_parent = int(Path(a.path).parent == Path(b.path).parent)
            same_root = int(a.root_index == b.root_index)
            score = (
                same_key,
                same_parent,
                same_root,
                min(a.pattern_score, b.pattern_score),
                -max(a.root_index, b.root_index),
                -abs(len(a.path) - len(b.path)),
                -len(a.path) - len(b.path),
            )
            pairs.append((score, a, b))
    if not pairs:
        raise ValueError("no paired combination")
    pairs.sort(key=lambda x: (x[0], x[1].path, x[2].path), reverse=True)
    return pairs[0][1], pairs[0][2], len(pairs) - 1


def choose_single(xs: list[Candidate]) -> tuple[Candidate, int]:
    if not xs:
        raise ValueError("no single-end candidate")
    xs = sorted(xs, key=lambda x: (x.pattern_score, -x.root_index, -len(x.path), x.path), reverse=True)
    return xs[0], len(xs) - 1


def signature(cfg: dict, manifest: pd.DataFrame) -> str:
    payload = {
        "species": species_id(cfg),
        "local_inputs": cfg.get("local_inputs", {}),
        "common_local_inputs": cfg.get("local_inputs", {}),
        "runs": manifest[["Run", "LibraryLayout"]].astype(str).to_dict("records"),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def existing_resolution_valid(path: Path) -> bool:
    try:
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except Exception:
        return False
    if df.empty or (df["status"] != "resolved").any():
        return False
    for col in ["r1", "r2", "single"]:
        if col not in df.columns:
            continue
        for p in df[col]:
            if p and (not Path(p).is_file() or Path(p).stat().st_size <= 0):
                return False
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description="Index and validate already-downloaded FASTQ files once per active manifest.")
    ap.add_argument("--config", required=True)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--deep-gzip-test", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    sid = species_id(cfg)
    manifest_path = workdir(cfg) / "manifest" / "manifest.evaluable.tsv"
    if not manifest_path.exists():
        raise SystemExit(f"[ERROR] active manifest missing: {manifest_path}")
    manifest = pd.read_csv(manifest_path, sep="\t", dtype={"Run": str}).fillna("")
    indexed_runs = set(manifest["Run"].astype(str))

    outdir = workdir(cfg) / "local_inputs"
    outdir.mkdir(parents=True, exist_ok=True)
    resolution_path = outdir / "read_resolution.tsv"
    sig_path = outdir / "read_index_config.sha256"
    done_path = outdir / ".local_inputs.done"
    sig = signature(cfg, manifest)
    if (
        not args.force and done_path.exists() and sig_path.exists()
        and sig_path.read_text().strip() == sig and existing_resolution_valid(resolution_path)
    ):
        print(f"[SKIP] local FASTQ index is current for {sid}: {resolution_path}")
        return

    local_cfg = cfg.get("local_inputs", {})
    roots = [Path(str(x)).expanduser() for x in local_cfg.get("read_roots", [])]
    if not roots:
        raise SystemExit(f"[ERROR] no local_inputs.read_roots configured for {sid}")
    missing_roots = [str(x) for x in roots if not x.exists()]
    if missing_roots:
        raise SystemExit("[ERROR] configured local read roots do not exist:\n  " + "\n  ".join(missing_roots))

    accepted = [str(x).lower() for x in local_cfg.get("accepted_extensions", [".fastq.gz", ".fq.gz", ".fastq", ".fq"])]
    excludes = [str(x).lower() for x in local_cfg.get("exclude_name_tokens", [])]
    recursive = bool(local_cfg.get("recursive_scan", True))
    followlinks = bool(local_cfg.get("follow_symlinks", False))

    candidates: list[Candidate] = []
    scanned = accepted_files = 0
    for root_index, root in enumerate(roots):
        print(f"[INFO] scanning local FASTQ root {root_index + 1}/{len(roots)}: {root}")
        for p in iter_fastq_files(root, recursive, followlinks):
            scanned += 1
            parsed = strip_fastq_extension(p.name, accepted)
            if parsed is None:
                continue
            if any(tok in p.name.lower() for tok in excludes):
                continue
            accepted_files += 1
            stem, _ = parsed
            m = RUN_SEARCH_RE.search(stem)
            if not m:
                continue
            run = m.group(0).upper()
            if run not in indexed_runs:
                continue
            classified = classify_stem(run, stem)
            if classified is None:
                continue
            role, pattern_score, pair_key = classified
            try:
                st = p.stat()
            except OSError:
                continue
            if st.st_size <= 0:
                continue
            candidates.append(Candidate(run, role, str(p.resolve()), root_index, pattern_score, pair_key, st.st_size, st.st_mtime_ns))

    by_run: dict[str, list[Candidate]] = {}
    for c in candidates:
        by_run.setdefault(c.run, []).append(c)

    validate_first = bool(local_cfg.get("validate_fastq_first_record", True))
    deep_test = args.deep_gzip_test or bool(local_cfg.get("deep_gzip_test", False))
    rows: list[dict[str, object]] = []
    for _, mrow in manifest.iterrows():
        run = str(mrow["Run"])
        layout = str(mrow["LibraryLayout"]).upper()
        cs = by_run.get(run, [])
        rec: dict[str, object] = {
            "Run": run, "LibraryLayout": layout, "status": "missing", "reason": "no_matching_local_fastq",
            "source_type": "local_fastq", "r1": "", "r2": "", "single": "",
            "r1_size_bytes": 0, "r2_size_bytes": 0, "single_size_bytes": 0,
            "alternate_candidate_count": 0, "first_qname": "", "first_pair_qname_match": "",
        }
        try:
            if layout == "PAIRED":
                a, b, alternates = choose_pair([x for x in cs if x.role == "r1"], [x for x in cs if x.role == "r2"])
                rec.update({
                    "status": "resolved", "reason": "", "r1": a.path, "r2": b.path,
                    "r1_size_bytes": a.size_bytes, "r2_size_bytes": b.size_bytes,
                    "alternate_candidate_count": alternates,
                })
                if validate_first:
                    q1, _ = first_record(Path(a.path)); q2, _ = first_record(Path(b.path))
                    rec["first_qname"] = q1
                    rec["first_pair_qname_match"] = str(q1 == q2).lower()
                if deep_test:
                    for p in [a.path, b.path]:
                        if p.lower().endswith(".gz"):
                            subprocess.run(["gzip", "-t", p], check=True)
            elif layout == "SINGLE":
                usable = [x for x in cs if x.role == "single"] or [x for x in cs if x.role == "r1"]
                a, alternates = choose_single(usable)
                rec.update({
                    "status": "resolved", "reason": "", "single": a.path,
                    "single_size_bytes": a.size_bytes, "alternate_candidate_count": alternates,
                })
                if validate_first:
                    q1, _ = first_record(Path(a.path)); rec["first_qname"] = q1
                if deep_test and a.path.lower().endswith(".gz"):
                    subprocess.run(["gzip", "-t", a.path], check=True)
            else:
                raise ValueError(f"unsupported LibraryLayout={layout}")
        except Exception as exc:
            rec["status"] = "invalid"
            rec["reason"] = str(exc)
        rows.append(rec)

    resolution = pd.DataFrame(rows)
    resolution.to_csv(resolution_path, sep="\t", index=False)
    pd.DataFrame([asdict(x) for x in candidates]).to_csv(outdir / "read_candidates.tsv", sep="\t", index=False)

    status_counts = resolution["status"].value_counts().to_dict()
    summary = pd.DataFrame([{
        "species": sid,
        "indexed_evaluable_runs": len(manifest),
        "resolved_runs": int((resolution.status == "resolved").sum()),
        "missing_runs": int((resolution.status == "missing").sum()),
        "invalid_runs": int((resolution.status == "invalid").sum()),
        "candidate_fastq_files": len(candidates),
        "filesystem_entries_scanned": scanned,
        "accepted_extension_files_seen": accepted_files,
        "read_roots": ";".join(map(str, roots)),
        "source_policy": local_cfg.get("read_source_policy", "prefer_local"),
    }])
    summary.to_csv(outdir / "read_inventory_summary.tsv", sep="\t", index=False)

    unresolved = resolution.loc[resolution.status != "resolved", ["Run", "LibraryLayout", "status", "reason"]]
    if not unresolved.empty:
        unresolved.to_csv(outdir / "unresolved_runs.tsv", sep="\t", index=False)
    else:
        (outdir / "unresolved_runs.tsv").write_text("Run\tLibraryLayout\tstatus\treason\n", encoding="utf-8")

    policy = str(local_cfg.get("read_source_policy", "prefer_local"))
    strict = policy == "local_required"
    if strict and not unresolved.empty:
        done_path.unlink(missing_ok=True)
        print(f"[ERROR] local FASTQ resolution failed for {len(unresolved)} of {len(manifest)} evaluable {sid} runs")
        print(unresolved.head(30).to_string(index=False))
        raise SystemExit(1)

    sig_path.write_text(sig + "\n", encoding="utf-8")
    done_path.write_text("ok\n", encoding="utf-8")
    print(f"[OK] local FASTQ index for {sid}: resolved={status_counts.get('resolved', 0)}/{len(manifest)} table={resolution_path}")


if __name__ == "__main__":
    main()
