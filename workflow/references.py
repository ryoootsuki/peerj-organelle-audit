#!/usr/bin/env python3
from __future__ import annotations

import sys
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

from common import ensure_dirs, load_config, resource_value, run_cmd, sha256, species_id, workdir


def source_stat(path: Path) -> dict[str, object]:
    st = path.stat()
    return {"path": str(path.resolve()), "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns}


def config_signature(cfg: dict) -> str:
    local = cfg.get("local_inputs", {}).get("references", {})
    payload: dict[str, object] = {
        "species": cfg["species"]["id"],
        "references": cfg["references"],
        "reference_preparation": cfg["reference_preparation"],
        "reference_source_policy": cfg.get("local_inputs", {}).get("reference_source_policy", "prefer_local"),
        "local_reference_config": local,
        "sources": [],
    }
    paths: list[Path] = []
    nuc = local.get("nuclear_fasta")
    if nuc:
        paths.append(Path(str(nuc)).expanduser())
    for spec in local.get("organelles", {}).values():
        if isinstance(spec, dict) and spec.get("fasta"):
            paths.append(Path(str(spec["fasta"])).expanduser())
    for p in paths:
        payload["sources"].append(source_stat(p) if p.exists() else {"path": str(p), "missing": True})
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def retry(cmd: list[str], log: Path, attempts: int = 4) -> None:
    for i in range(1, attempts + 1):
        try:
            run_cmd(cmd, log=log)
            return
        except Exception:
            if i == attempts:
                raise
            import time
            time.sleep(i * 10)


def clean_fasta(src: Path, dst: Path, log: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    run_cmd(["seqkit", "seq", "--only-id", str(src), "-o", str(dst)], log=log)
    if not dst.exists() or dst.stat().st_size == 0:
        raise RuntimeError(f"seqkit produced an empty FASTA: {dst}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    ensure_dirs(cfg)
    sid = species_id(cfg)
    refs_cfg = cfg["references"]
    prep_cfg = cfg["reference_preparation"]
    local_cfg = cfg.get("local_inputs", {})
    local_refs = local_cfg.get("references", {})
    source_policy = str(local_cfg.get("reference_source_policy", "prefer_local"))
    allow_download = bool(cfg.get("network", {}).get("allow_reference_download", True))
    ref_threads = int(resource_value(cfg, "reference_threads", None, 8))

    refdir = workdir(cfg) / "references"
    source = refdir / "source"
    prepared = refdir / "prepared"
    report = refdir / "reports"
    log = Path("logs") / sid / "03_references.log"
    for p in [source, prepared, report, log.parent]:
        p.mkdir(parents=True, exist_ok=True)

    sig = config_signature(cfg)
    expected = prepared / "nuclear.masked.main.fa" if refs_cfg.get("nuclear_assembly_accession") else prepared / "organelles.combined.fa"
    if (
        not args.force
        and (report / ".references.done").exists()
        and (report / "reference_config.sha256").exists()
        and (report / "reference_config.sha256").read_text().strip() == sig
        and expected.exists() and expected.stat().st_size > 0
    ):
        print(f"[SKIP] prepared references are current for {sid}")
        return

    (report / ".references.done").unlink(missing_ok=True)
    provenance: list[dict[str, object]] = []

    nuclear_accession = str(refs_cfg.get("nuclear_assembly_accession") or "")
    nuclear_available = bool(nuclear_accession)
    nuclear_clean: Path | None = None
    local_nuclear = local_refs.get("nuclear_fasta")

    if nuclear_available:
        if local_nuclear:
            nuc_src = Path(str(local_nuclear)).expanduser()
            if not nuc_src.is_file() or nuc_src.stat().st_size <= 0:
                raise SystemExit(f"[ERROR] configured local nuclear FASTA is missing or empty: {nuc_src}")
            source_kind = "local_fasta"
        else:
            if source_policy == "local_required" or not allow_download:
                raise SystemExit(f"[ERROR] local nuclear FASTA is required for {sid} but local_inputs.references.nuclear_fasta is not usable")
            zip_path = source / f"{nuclear_accession}.zip"
            dataset_dir = source / f"{nuclear_accession}_dataset"
            if not zip_path.exists():
                retry(["datasets", "download", "genome", "accession", nuclear_accession, "--include", "genome,seq-report", "--filename", str(zip_path)], log)
            if not dataset_dir.exists():
                dataset_dir.mkdir(parents=True)
                run_cmd(["unzip", "-q", "-o", str(zip_path), "-d", str(dataset_dir)], log=log)
            hits = sorted([p for p in dataset_dir.rglob("*") if p.is_file() and (p.name.endswith("genomic.fna") or p.suffix in {".fna", ".gz"})])
            if not hits:
                raise RuntimeError(f"nuclear FASTA not found in downloaded dataset for {nuclear_accession}")
            nuc_src = hits[0]
            source_kind = "NCBI_datasets_download"
        nuclear_clean = source / "nuclear.clean.fa"
        clean_fasta(nuc_src, nuclear_clean, log)
        provenance.append({"component": "nuclear", "accession": nuclear_accession, "source_kind": source_kind, **source_stat(nuc_src)})

    organelle_fastas: dict[str, Path] = {}
    organelle_local_specs = local_refs.get("organelles", {})
    for name, accession_value in refs_cfg.get("organelles", {}).items():
        accession = str(accession_value)
        spec = organelle_local_specs.get(name, {}) if isinstance(organelle_local_specs, dict) else {}
        local_fasta = spec.get("fasta") if isinstance(spec, dict) else None
        extract_acc = spec.get("extract_accession_from_nuclear") if isinstance(spec, dict) else None
        raw = source / f"{name}.{accession}.source.fa"
        if local_fasta:
            org_src = Path(str(local_fasta)).expanduser()
            if not org_src.is_file() or org_src.stat().st_size <= 0:
                raise SystemExit(f"[ERROR] configured local {name} FASTA is missing or empty: {org_src}")
            source_kind = "local_fasta"
            clean_fasta(org_src, prepared / f"{name}.fa", log)
            provenance.append({"component": name, "accession": accession, "source_kind": source_kind, **source_stat(org_src)})
        elif extract_acc:
            if nuclear_clean is None:
                raise SystemExit(f"[ERROR] {name} extraction requested but no nuclear assembly is available")
            extract_acc = str(extract_acc)
            # nuclear.clean.fa was created with SeqKit --only-id, so an exact
            # identifier match is safer and more portable than a Perl-style
            # regular expression (Go/RE2 does not support all PCRE syntax).
            run_cmd(["seqkit", "grep", "-p", extract_acc, str(nuclear_clean), "-o", str(raw)], log=log)
            if not raw.exists() or raw.stat().st_size <= 0:
                raise RuntimeError(f"failed to extract {extract_acc} from local nuclear assembly {nuclear_clean}")
            clean_fasta(raw, prepared / f"{name}.fa", log)
            provenance.append({
                "component": name, "accession": accession,
                "source_kind": "extracted_from_local_nuclear_assembly",
                **source_stat(raw),
                "parent_source_path": str(Path(str(local_nuclear)).expanduser().resolve()) if local_nuclear else str(nuclear_clean.resolve()),
            })
        else:
            if source_policy == "local_required" or not allow_download:
                raise SystemExit(f"[ERROR] local source for {name} ({accession}) is required but no fasta/extraction rule is configured")
            retry(["bash", "-o", "pipefail", "-c", f"efetch -db nucleotide -id {accession} -format fasta > {raw}"], log)
            clean_fasta(raw, prepared / f"{name}.fa", log)
            provenance.append({"component": name, "accession": accession, "source_kind": "NCBI_efetch_download", **source_stat(raw)})
        organelle_fastas[name] = prepared / f"{name}.fa"

    if not organelle_fastas:
        raise SystemExit("[ERROR] at least one organelle reference is required")
    with (prepared / "organelles.combined.fa").open("wb") as out:
        for p in organelle_fastas.values():
            with p.open("rb") as src:
                shutil.copyfileobj(src, out)

    whole_cov = float(prep_cfg["whole_contig_min_coverage"])
    whole_id = float(prep_cfg["whole_contig_min_identity"])
    main_len = int(prep_cfg["main_mask_min_alignment_length"])
    main_id = float(prep_cfg["main_mask_min_identity"])
    strict_len = int(prep_cfg["strict_mask_min_alignment_length"])
    strict_id = float(prep_cfg["strict_mask_min_identity"])
    mm_preset = str(prep_cfg.get("minimap2_preset", "asm20"))
    mm_n = int(prep_cfg.get("minimap2_max_secondary", 10000))
    mm_p = float(prep_cfg.get("minimap2_secondary_score_ratio", 0.20))

    if nuclear_available and nuclear_clean is not None:
        explicit = [str(x) for x in refs_cfg.get("explicit_organelle_accessions_to_remove", [])]
        run_cmd([sys.executable, "workflow/fasta_tools.py", "header-ids", "--fasta", str(nuclear_clean), "--output", str(report / "header_organelle_contigs.ids"), "--accessions", *explicit], log=log)
        paf_raw = report / "organelles_vs_nuclear_raw.paf"
        with paf_raw.open("w", encoding="utf-8") as out, log.open("a", encoding="utf-8") as err:
            subprocess.run(["minimap2", "-x", mm_preset, "--secondary=yes", "-N", str(mm_n), "-p", str(mm_p), "-t", str(ref_threads), str(nuclear_clean), str(prepared / "organelles.combined.fa")], check=True, stdout=out, stderr=err)
        run_cmd([sys.executable, "workflow/paf_intervals.py", "--paf", str(paf_raw), "--whole-contigs", str(report / "whole_contig_alignment_report.tsv"), "--whole-min-coverage", str(whole_cov), "--whole-min-identity", str(whole_id)], log=log)

        aligned_ids: list[str] = []
        whole_lines = (report / "whole_contig_alignment_report.tsv").read_text().splitlines()
        for line in whole_lines[1:]:
            f = line.split("\t")
            if len(f) >= 6 and f[5].lower() == "true":
                aligned_ids.append(f[0])
        (report / "alignment_organelle_contigs.ids").write_text("\n".join(sorted(set(aligned_ids))) + ("\n" if aligned_ids else ""), encoding="utf-8")
        header_ids = set(x.strip() for x in (report / "header_organelle_contigs.ids").read_text().splitlines() if x.strip())
        removed = sorted(header_ids | set(aligned_ids))
        (report / "removed_organelle_contigs.ids").write_text("\n".join(removed) + ("\n" if removed else ""), encoding="utf-8")
        contig_filtered = prepared / "nuclear.contig_filtered.fa"
        if removed:
            run_cmd(["seqkit", "grep", "-v", "-f", str(report / "removed_organelle_contigs.ids"), str(nuclear_clean), "-o", str(contig_filtered)], log=log)
        else:
            shutil.copy2(nuclear_clean, contig_filtered)

        paf_filtered = report / "organelles_vs_contig_filtered.paf"
        with paf_filtered.open("w", encoding="utf-8") as out, log.open("a", encoding="utf-8") as err:
            subprocess.run(["minimap2", "-x", mm_preset, "--secondary=yes", "-N", str(mm_n), "-p", str(mm_p), "-t", str(ref_threads), str(contig_filtered), str(prepared / "organelles.combined.fa")], check=True, stdout=out, stderr=err)
        for mode, min_len, min_ident in [("main", main_len, main_id), ("strict", strict_len, strict_id)]:
            raw_bed = report / f"mask.{mode}.raw.bed"
            merged_bed = report / f"mask.{mode}.merged.bed"
            masked = prepared / f"nuclear.masked.{mode}.fa"
            run_cmd([sys.executable, "workflow/paf_intervals.py", "--paf", str(paf_filtered), "--mask-bed", str(raw_bed), "--mask-min-length", str(min_len), "--mask-min-identity", str(min_ident)], log=log)
            if raw_bed.exists() and raw_bed.stat().st_size > 0:
                cmd = f"cut -f1-3 {raw_bed!s} | sort -k1,1 -k2,2n | bedtools merge -i - > {merged_bed!s}"
                run_cmd(["bash", "-o", "pipefail", "-c", cmd], log=log)
                run_cmd(["bedtools", "maskfasta", "-fi", str(contig_filtered), "-bed", str(merged_bed), "-fo", str(masked)], log=log)
            else:
                merged_bed.write_text("", encoding="utf-8")
                shutil.copy2(contig_filtered, masked)
            run_cmd([sys.executable, "workflow/fasta_tools.py", "mask-summary", "--bed", str(merged_bed), "--name", mode, "--output", str(report / f"mask.{mode}.summary.tsv")], log=log)
    else:
        for name in ["header_organelle_contigs.ids", "alignment_organelle_contigs.ids", "removed_organelle_contigs.ids", "mask.main.merged.bed", "mask.strict.merged.bed"]:
            (report / name).write_text("", encoding="utf-8")

    index_fastas = list(organelle_fastas.values())
    if nuclear_available:
        index_fastas += [prepared / "nuclear.contig_filtered.fa", prepared / "nuclear.masked.main.fa", prepared / "nuclear.masked.strict.fa"]
    for fa in index_fastas:
        if not (fa.exists() and fa.stat().st_size > 0):
            raise RuntimeError(f"prepared reference missing: {fa}")
        if not Path(str(fa) + ".fai").exists():
            run_cmd(["samtools", "faidx", str(fa)], log=log)
        if not Path(str(fa) + ".bwt.2bit.64").exists():
            run_cmd(["bwa-mem2", "index", str(fa)], log=log)

    stats_inputs = list(organelle_fastas.values())
    if nuclear_available and nuclear_clean is not None:
        stats_inputs = [nuclear_clean, prepared / "nuclear.contig_filtered.fa", prepared / "nuclear.masked.main.fa", prepared / "nuclear.masked.strict.fa", *organelle_fastas.values()]
    with (report / "reference_seqkit_stats.tsv").open("w", encoding="utf-8") as out:
        subprocess.run(["seqkit", "stats", "-T", *map(str, stats_inputs)], check=True, stdout=out)
    with (report / "reference_sha256.txt").open("w", encoding="utf-8") as out:
        for fa in stats_inputs:
            out.write(f"{sha256(fa)}  {fa}\n")
    for rec in provenance:
        p = Path(str(rec["path"]))
        rec["sha256"] = sha256(p) if p.exists() else "not_available"
    cols = ["component", "accession", "source_kind", "path", "parent_source_path", "size_bytes", "mtime_ns", "sha256"]
    with (report / "reference_source_provenance.tsv").open("w", encoding="utf-8") as out:
        out.write("\t".join(cols) + "\n")
        for rec in provenance:
            out.write("\t".join(str(rec.get(c, "")) for c in cols) + "\n")
    with (report / "reference_provenance.tsv").open("w", encoding="utf-8") as out:
        out.write("species\tnuclear_available\tnuclear_assembly\torganelle_references\treference_source_policy\twhole_contig_min_coverage\twhole_contig_min_identity\tmain_mask_min_length\tmain_mask_min_identity\tstrict_mask_min_length\tstrict_mask_min_identity\tminimap2_preset\tminimap2_max_secondary\tminimap2_secondary_score_ratio\treference_threads\n")
        orgs = ",".join(f"{k}={v}" for k, v in refs_cfg.get("organelles", {}).items())
        out.write(f"{sid}\t{str(nuclear_available).lower()}\t{nuclear_accession or 'not_available'}\t{orgs}\t{source_policy}\t{whole_cov}\t{whole_id}\t{main_len}\t{main_id}\t{strict_len}\t{strict_id}\t{mm_preset}\t{mm_n}\t{mm_p}\t{ref_threads}\n")
    (report / "reference_config.sha256").write_text(sig + "\n", encoding="utf-8")
    (report / ".references.done").write_text("ok\n", encoding="utf-8")
    print(f"[OK] local references prepared for {sid}: {prepared}")


if __name__ == "__main__":
    main()
