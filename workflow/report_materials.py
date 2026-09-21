#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import pandas as pd
import yaml


def read_tsv(path: str):
    p = Path(path)
    return pd.read_csv(p, sep="\t") if p.exists() and p.stat().st_size else None


def fmt_pct(x):
    return f"{100 * float(x):.3g}%"


def main():
    out = Path("results/combined/report_materials")
    out.mkdir(parents=True, exist_ok=True)
    common = yaml.safe_load(Path("config/common.yaml").read_text())
    primary = common["workflow"]["primary_mapq"]
    pp = common["preprocess"]
    ref = common["reference_preparation"]

    methods = f"""# Workflow methods summary\n\nThe analysis uses pinned run manifests, publication-curated BioProject metadata, synchronized FASTQ preprocessing, independently prepared reference classes, and query-name-based compatibility sets. fastp uses a qualified-quality threshold of {pp['qualified_quality_phred']} and minimum retained length of {pp['length_required']} bp. Nuclear reference preparation first removes organelle-like records at the configured whole-contig thresholds and then masks residual organelle-like intervals. The primary mask requires at least {ref['main_mask_min_alignment_length']} bp at {ref['main_mask_min_identity']:.0%} identity; a higher-confidence sensitivity mask requires at least {ref['strict_mask_min_alignment_length']} bp at {ref['strict_mask_min_identity']:.0%} identity. Reads are evaluated across MAPQ {', '.join(map(str, common['workflow']['mapq_thresholds']))}, with MAPQ >= {primary} used as the primary reporting threshold.\n\nThe workflow reports inclusive organelle-compatible, strict organelle-compatible, and nuclear-organelle multi-reference fractions when a nuclear reference is available. In organelle-only mode, quantities requiring a nuclear comparator remain undefined rather than being reported as zero. Downstream alignment/depth and variant outputs are sensitivity analyses and are not treated as truth-set validation.\n"""
    (out / "methods_summary.md").write_text(methods, encoding="utf-8")

    lines = ["# Compact results summary", ""]
    for sid, label in [("silene", "Silene latifolia"), ("chicken", "Gallus gallus")]:
        d = read_tsv(f"results/{sid}/tables/audit_primary_mapq.tsv")
        if d is None:
            lines.append(f"- {label}: primary audit table not present in this run directory.")
            continue
        lines.append(f"- {label}: {d.Run.nunique()} evaluable runs at MAPQ >= {primary}; median strict fraction {fmt_pct(d.frac_strict_total.median())} and median multi-reference fraction {fmt_pct(d.frac_multi_total.median())}.")
    lines.append("")
    lines.append("Interpretation: these values quantify reference compatibility under the stated references and thresholds. They do not establish the biological origin of individual reads.")
    (out / "results_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    checks = []
    def add(item, path):
        p = Path(path)
        ok = p.exists() and p.stat().st_size > 0
        checks.append({"output": item, "status": "PASS" if ok else "MISSING", "path": str(path)})
    add("Pinned Silene manifest", "manifests/silene_SraRunInfo_retrieved_2025-12-25.csv")
    add("Pinned chicken manifest", "manifests/chicken_PRJNA573756_exact_141_runs.tsv")
    for sid in ["silene", "chicken"]:
        add(f"{sid} audit table", f"results/{sid}/tables/audit_primary_mapq.tsv")
        add(f"{sid} MAPQ summary", f"results/{sid}/tables/mapq_distribution_summary.tsv")
        add(f"{sid} protocol metadata", f"results/{sid}/tables/library_protocol_metadata.tsv")
        add(f"{sid} downstream sensitivity", f"results/{sid}/tables/downstream_variant_deltas.tsv")
        add(f"{sid} reference-mask sensitivity", f"results/{sid}/tables/reference_mask_sensitivity.tsv")
    add("Cross-species summary", "results/combined/species_validation_summary.tsv")
    pd.DataFrame(checks).to_csv(out / "computational_output_checklist.tsv", sep="\t", index=False)
    print(f"[OK] report materials written: {out}")


if __name__ == "__main__":
    main()
