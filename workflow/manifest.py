#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from common import RUN_RE, ensure_dirs, load_config, results_dir, sha256, workdir

ALIASES = {
    "Run": ["Run", "run_accession", "run"],
    "spots": ["spots", "read_count", "spots_count"],
    "LibraryStrategy": ["LibraryStrategy", "library_strategy"],
    "LibraryLayout": ["LibraryLayout", "library_layout"],
    "Platform": ["Platform", "instrument_platform"],
    "Model": ["Model", "instrument_model"],
    "BioProject": ["BioProject", "study_accession", "bioproject"],
    "ScientificName": ["ScientificName", "scientific_name"],
    "Experiment": ["Experiment", "experiment_accession"],
    "BioSample": ["BioSample", "sample_accession"],
    "SampleName": ["SampleName", "sample_alias"],
    "Sex": ["Sex", "sex"],
    "LibrarySelection": ["LibrarySelection", "library_selection"],
    "LibrarySource": ["LibrarySource", "library_source"],
    "avgLength": ["avgLength", "average_length"],
    "bases": ["bases", "base_count"],
}


def pick_column(df: pd.DataFrame, choices: list[str]) -> str | None:
    for c in choices:
        if c in df.columns:
            return c
    return None


def _read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"[ERROR] metadata curation table is missing: {p}")
    return pd.read_csv(p, sep="\t", dtype=str).fillna("")


def apply_strategy_curation(out: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    out = out.copy()
    out["ArchiveLibraryStrategy"] = out["LibraryStrategy"].astype(str)
    out["LibraryStrategy"] = out["ArchiveLibraryStrategy"]
    out["StrategyCurationSource"] = ""
    out["StrategyCurationNote"] = ""
    path = cfg.get("metadata_curation", {}).get("strategy_overrides")
    if not path:
        return out
    cur = _read_table(path)
    required = {"BioProject", "CuratedLibraryStrategy"}
    missing = required - set(cur.columns)
    if missing:
        raise SystemExit(f"[ERROR] strategy override table lacks columns: {sorted(missing)}")
    if cur["BioProject"].duplicated().any():
        raise SystemExit("[ERROR] duplicate BioProject rows in strategy override table")
    cols = [c for c in ["BioProject", "CuratedLibraryStrategy", "curation_basis", "curation_note"] if c in cur.columns]
    out = out.merge(cur[cols], on="BioProject", how="left", validate="many_to_one")
    mask = out["CuratedLibraryStrategy"].fillna("").astype(str).str.strip().ne("")
    out.loc[mask, "LibraryStrategy"] = out.loc[mask, "CuratedLibraryStrategy"].astype(str)
    if "curation_basis" in out:
        out.loc[mask, "StrategyCurationSource"] = out.loc[mask, "curation_basis"].astype(str)
    if "curation_note" in out:
        out.loc[mask, "StrategyCurationNote"] = out.loc[mask, "curation_note"].astype(str)
    return out.drop(columns=[c for c in ["CuratedLibraryStrategy", "curation_basis", "curation_note"] if c in out.columns])


def write_protocol_metadata(out: pd.DataFrame, cfg: dict) -> None:
    cols = [c for c in [
        "Run", "BioProject", "Experiment", "BioSample", "SampleName",
        "ArchiveLibraryStrategy", "LibraryStrategy", "LibraryLayout", "LibrarySelection", "LibrarySource",
        "Platform", "Model", "ScientificName", "Sex", "StrategyCurationSource", "StrategyCurationNote"
    ] if c in out.columns]
    protocol = out[cols].copy()
    path = cfg.get("metadata_curation", {}).get("protocol_table")
    if path:
        project = _read_table(path)
        if "BioProject" not in project.columns or project["BioProject"].duplicated().any():
            raise SystemExit("[ERROR] protocol table must contain one unique row per BioProject")
        protocol = protocol.merge(project, on="BioProject", how="left", validate="many_to_one")
        protocol["protocol_metadata_source"] = "publication_curated_BioProject_table"
    else:
        for col in ["tissue", "dna_extraction", "restriction_enzyme_1", "restriction_enzyme_2", "fragment_selection", "pcr_cycles", "primary_publication", "publication_doi"]:
            protocol[col] = "not_reported_in_pinned_sra_manifest"
        protocol["manual_verification_required"] = True
        protocol["protocol_metadata_source"] = "pinned_archive_metadata_only"
    protocol["metadata_interpretation"] = "descriptive_metadata; missing protocol fields are not inferred"
    protocol.to_csv(results_dir(cfg) / "tables" / "library_protocol_metadata.tsv", sep="\t", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--mode", choices=["full", "pilot"], default="full")
    args = ap.parse_args()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    src = Path(cfg["species"]["exact_manifest"])
    if not src.is_file():
        raise SystemExit(f"[ERROR] pinned manifest is missing: {src}")
    expected_hash = str(cfg["species"].get("manifest_sha256", "")).strip()
    if expected_hash:
        observed_hash = sha256(src)
        if observed_hash != expected_hash:
            raise SystemExit(f"[ERROR] pinned manifest SHA-256 mismatch for {src}: observed={observed_hash}, expected={expected_hash}")
    fmt = cfg["species"].get("manifest_format", "tsv")
    raw = pd.read_csv(src, low_memory=False) if fmt == "sra_runinfo_csv" else pd.read_csv(src, sep="\t", low_memory=False)

    out = pd.DataFrame(index=raw.index)
    for canonical, aliases in ALIASES.items():
        col = pick_column(raw, aliases)
        out[canonical] = raw[col] if col else ""
    out["Run"] = out["Run"].astype(str).str.strip()
    invalid = out.loc[~out["Run"].map(lambda x: bool(RUN_RE.fullmatch(x))), "Run"].tolist()
    if invalid:
        raise SystemExit(f"[ERROR] invalid run accessions in pinned manifest: {invalid[:20]}")
    if out["Run"].duplicated().any():
        raise SystemExit(f"[ERROR] duplicate runs in pinned manifest: {out.loc[out['Run'].duplicated(False), 'Run'].tolist()[:20]}")

    out["spots"] = pd.to_numeric(out["spots"], errors="coerce").fillna(0).astype("int64")
    out["bases"] = pd.to_numeric(out["bases"], errors="coerce").fillna(0).astype("int64")
    out["avgLength"] = pd.to_numeric(out["avgLength"], errors="coerce").fillna(0)
    out["LibraryLayout"] = out["LibraryLayout"].astype(str).str.upper().replace({"NAN": "UNKNOWN", "": "UNKNOWN"})
    out["LibraryStrategy"] = out["LibraryStrategy"].astype(str).replace({"nan": "not_reported", "": "not_reported"})
    out["ArchiveLibraryLayout"] = out["LibraryLayout"]
    out = apply_strategy_curation(out, cfg)
    out["evaluable"] = out["spots"] > 0
    out["exclusion_reason"] = ""
    out.loc[~out["evaluable"], "exclusion_reason"] = "metadata_spots_zero"

    expected_total = int(cfg["species"]["expected_total_runs"])
    expected_eval = int(cfg["species"]["expected_evaluable_runs"])
    if len(out) != expected_total:
        raise SystemExit(f"[ERROR] pinned manifest run count {len(out)} != expected {expected_total}")
    if int(out["evaluable"].sum()) != expected_eval:
        raise SystemExit(f"[ERROR] evaluable run count {int(out['evaluable'].sum())} != expected {expected_eval}")
    for strategy, expected in cfg["species"].get("expected_library_strategy_counts", {}).items():
        observed = int((out["LibraryStrategy"] == strategy).sum())
        if observed != int(expected):
            raise SystemExit(f"[ERROR] curated {strategy} count {observed} != expected {expected}")

    manifest_dir = workdir(cfg) / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(manifest_dir / "manifest.all.tsv", sep="\t", index=False)
    evaluable = out.loc[out["evaluable"]].copy()
    evaluable.to_csv(manifest_dir / "manifest.evaluable.tsv", sep="\t", index=False)

    if args.mode == "pilot":
        n = int(cfg["workflow"]["pilot_runs_per_species"])
        sort_cols = [c for c in ["BioProject", "LibraryLayout", "LibraryStrategy", "Run"] if c in evaluable.columns]
        groups = [c for c in ["BioProject", "LibraryLayout", "LibraryStrategy"] if c in evaluable.columns]
        picks = [sub.iloc[0] for _, sub in evaluable.sort_values(sort_cols).groupby(groups, dropna=False, sort=True)] if groups else []
        pilot = pd.DataFrame(picks)
        remaining = evaluable.loc[~evaluable["Run"].isin(pilot.get("Run", []))].sort_values("Run")
        if len(pilot) < n:
            pilot = pd.concat([pilot, remaining.head(n - len(pilot))], ignore_index=True)
        active = pilot.head(n)
    else:
        active = evaluable
    active.to_csv(manifest_dir / "manifest.active.tsv", sep="\t", index=False)

    out.groupby(["LibraryStrategy", "LibraryLayout"], dropna=False).size().reset_index(name="n_runs").to_csv(
        results_dir(cfg) / "tables" / "manifest_summary.tsv", sep="\t", index=False
    )
    write_protocol_metadata(out, cfg)
    print(f"[OK] {cfg['species']['id']} pinned manifest validated: total={len(out)}, evaluable={len(evaluable)}, active={len(active)}")

if __name__ == "__main__":
    main()
