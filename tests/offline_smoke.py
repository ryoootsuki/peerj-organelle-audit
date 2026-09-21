#!/usr/bin/env python3
"""Fast, network-free tests for manifests, local input indexing, PAF parsing, and SAM qname extraction."""
from __future__ import annotations

import gzip
import os
import py_compile
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
WF = ROOT / "workflow"
ENV = os.environ.copy()
ENV["PYTHONPATH"] = str(WF)


def run(cmd, *, cwd=ROOT, input_text=None):
    return subprocess.run(cmd, cwd=cwd, env=ENV, text=True, input=input_text, check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def write_fastq(path: Path, name: str):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        fh.write(f"@{name}\nACGTACGT\n+\nFFFFFFFF\n")


def main():
    for path in sorted(WF.glob("*.py")):
        py_compile.compile(str(path), doraise=True)
    print(f"[OK] Python syntax: {len(list(WF.glob('*.py')))} modules")

    sil = pd.read_csv(ROOT / "manifests/silene_SraRunInfo_retrieved_2025-12-25.csv", low_memory=False)
    chick = pd.read_csv(ROOT / "manifests/chicken_PRJNA573756_exact_141_runs.tsv", sep="\t")
    assert len(sil) == 2686
    assert (pd.to_numeric(sil["spots"], errors="coerce").fillna(0) > 0).sum() == 2679
    assert len(chick) == 141
    override = pd.read_csv(ROOT / "metadata/bioproject_strategy_overrides.tsv", sep="\t", dtype=str).fillna("")
    r = override.loc[override.BioProject.eq("PRJNA285775")]
    assert len(r) == 1 and r.iloc[0].CuratedLibraryStrategy == "RAD-Seq"
    print("[OK] pinned manifests: Silene=2686/2679 evaluable, chicken=141; PRJNA285775 curated as RAD-Seq")

    with tempfile.TemporaryDirectory(prefix="peerj_v2_test_") as td_s:
        td = Path(td_s)
        (td / "config").mkdir(); (td / "manifests").mkdir(); (td / "local_reads").mkdir()
        common = {
            "workflow": {"pilot_runs_per_species": 2, "mapq_thresholds": [0, 20], "primary_mapq": 20},
            "network": {"allow_read_download": False, "allow_reference_download": False},
            "local_inputs": {
                "read_source_policy": "local_required", "reference_source_policy": "local_required",
                "recursive_scan": True, "follow_symlinks": False, "validate_fastq_first_record": True,
                "deep_gzip_test": False, "ambiguity_policy": "deterministic_warn",
                "exclude_name_tokens": ["_trim", ".trim", "_paired", ".paired", "singletons"],
                "accepted_extensions": [".fastq.gz", ".fq.gz", ".fastq", ".fq"],
                "read_hash_mode": "none", "index_manifest_scope": "evaluable",
            },
            "resources": {"tmp_root": "tmp", "lock_root": "work/.locks"},
        }
        (td / "config/common.yaml").write_text(yaml.safe_dump(common, sort_keys=False))
        manifest = pd.DataFrame([
            {"Run": "SRR123456", "spots": 1, "LibraryStrategy": "RAD-Seq", "LibraryLayout": "PAIRED"},
            {"Run": "SRR123457", "spots": 1, "LibraryStrategy": "WGS", "LibraryLayout": "SINGLE"},
        ])
        manifest.to_csv(td / "manifests/test.tsv", sep="\t", index=False)
        override = pd.DataFrame([{"BioProject":"TESTPRJ","ArchiveLibraryStrategy":"WGS","CuratedLibraryStrategy":"RAD-Seq","curation_basis":"synthetic test","curation_note":"test"}])
        # Give the synthetic WGS row a project so curation can be tested.
        manifest.loc[manifest.Run.eq("SRR123457"), "BioProject"] = "TESTPRJ"
        manifest.to_csv(td / "manifests/test.tsv", sep="\t", index=False)
        override.to_csv(td / "manifests/override.tsv", sep="\t", index=False)
        species = {
            "species": {"id": "localtest", "scientific_name": "Test species", "mode": "plant_three_reference",
                        "exact_manifest": "manifests/test.tsv", "manifest_format": "tsv",
                        "expected_total_runs": 2, "expected_evaluable_runs": 2},
            "metadata_curation": {"strategy_overrides": "manifests/override.tsv"},
            "references": {"nuclear_assembly_accession": "TEST", "organelles": {"mitochondrial": "TESTMT"}},
            "local_inputs": {"read_roots": [str(td / "local_reads")], "references": {}},
        }
        (td / "config/test.yaml").write_text(yaml.safe_dump(species, sort_keys=False))
        write_fastq(td / "local_reads/SRR123456_1.fastq.gz", "SRR123456.1/1")
        write_fastq(td / "local_reads/SRR123456_2.fastq.gz", "SRR123456.1/2")
        write_fastq(td / "local_reads/SRR123457.fastq.gz", "SRR123457.1")
        run([sys.executable, str(WF / "manifest.py"), "--config", "config/test.yaml", "--mode", "full"], cwd=td)
        curated = pd.read_csv(td / "work/localtest/manifest/manifest.all.tsv", sep="\t", dtype=str).fillna("")
        cr = curated.loc[curated.Run.eq("SRR123457")].iloc[0]
        assert cr.ArchiveLibraryStrategy == "WGS" and cr.LibraryStrategy == "RAD-Seq"
        run([sys.executable, str(WF / "local_inputs.py"), "--config", "config/test.yaml"], cwd=td)
        resolved = pd.read_csv(td / "work/localtest/local_inputs/read_resolution.tsv", sep="\t", dtype=str).fillna("")
        assert len(resolved) == 2 and resolved["status"].str.strip().eq("resolved").all(), resolved.to_string()
        assert set(resolved["Run"]) == {"SRR123456", "SRR123457"}
        print("[OK] local FASTQ indexing: paired + single resolved")

        paf = td / "test.paf"
        paf.write_text("mt\t100\t0\t100\t+\tcontigA\t100\t0\t100\t98\t100\t60\n")
        run([sys.executable, str(WF / "paf_intervals.py"), "--paf", str(paf),
             "--whole-contigs", str(td / "whole.tsv"), "--mask-bed", str(td / "mask.bed"),
             "--whole-min-coverage", "0.8", "--whole-min-identity", "0.95",
             "--mask-min-length", "90", "--mask-min-identity", "0.9"], cwd=td)
        whole = pd.read_csv(td / "whole.tsv", sep="\t")
        assert str(whole.iloc[0]["remove"]).lower() == "true"
        print("[OK] PAF contig-removal and mask interval logic")

        sam = (
            "@HD\tVN:1.6\n"
            "r1\t0\tref\t1\t30\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\n"
            "r1\t256\tref\t1\t60\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\n"
            "r2\t0\tref\t1\t10\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\n"
        )
        qdir = td / "qsets"
        run([sys.executable, str(WF / "sam_qnames.py"), "--output-dir", str(qdir), "--tag", "nuclear",
             "--thresholds", "0", "20", "--summary", str(qdir / "summary.json")], cwd=td, input_text=sam)
        with gzip.open(qdir / "nuclear.q0.qnames.txt.gz", "rt") as fh: q0 = {x.strip() for x in fh if x.strip()}
        with gzip.open(qdir / "nuclear.q20.qnames.txt.gz", "rt") as fh: q20 = {x.strip() for x in fh if x.strip()}
        assert q0 == {"r1", "r2"} and q20 == {"r1"}
        print("[OK] one-pass SAM-to-MAPQ qname extraction")

        # Organelle-only mode: inclusive compatibility is available,
        # while strict/multi-reference metrics must remain NA without a nuclear reference.
        import importlib.util
        if str(WF) not in sys.path:
            sys.path.insert(0, str(WF))
        spec = importlib.util.spec_from_file_location("audit_worker", WF / "audit_worker.py")
        audit_worker = importlib.util.module_from_spec(spec); spec.loader.exec_module(audit_worker)
        org_run = td / "org_only"
        (org_run / "sets").mkdir(parents=True, exist_ok=True)
        for q, names in [(0, ["r1", "r2"]), (20, ["r1"])]:
            with gzip.open(org_run / "sets" / f"mitochondrial.q{q}.qnames.txt.gz", "wt") as fh:
                for name in names: fh.write(name + "\n")
        org_out = td / "organelle_only.tsv"
        audit_worker.classify("ORGONLY", 10, ["mitochondrial"], [0, 20], org_run, org_out, False)
        od = pd.read_csv(org_out, sep="\t")
        assert od.reference_mode.eq("organelle_only").all()
        assert abs(float(od.loc[od.MAPQ.eq(0), "frac_inclusive_total"].iloc[0]) - 0.2) < 1e-12
        assert abs(float(od.loc[od.MAPQ.eq(20), "frac_inclusive_total"].iloc[0]) - 0.1) < 1e-12
        assert od.frac_strict_total.isna().all() and od.frac_multi_total.isna().all()
        print("[OK] organelle-only mode: inclusive reported; strict/multi remain NA")

    print("[OK] all network-free smoke tests passed")


if __name__ == "__main__":
    main()
