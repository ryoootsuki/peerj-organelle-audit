#!/usr/bin/env python3
from __future__ import annotations

import sys
import argparse
import os
import shlex
import shutil
import subprocess
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from common import (
    atomic_json,
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


def find_sra(root: Path, run: str) -> Path:
    hits = list(root.rglob(f"{run}.sra"))
    if not hits:
        raise FileNotFoundError(f"SRA file missing after prefetch: {run}")
    return hits[0]


def write_names(path: Path, names: set[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for x in sorted(names):
            f.write(x + "\n")


def qnames_from_bam(bam: Path, outdir: Path, tag: str, thresholds: list[int],
                    threads: int, compression: int, log: Path) -> None:
    sam = shlex.join(["samtools", "view", "-@", str(threads), "-F", "2308", str(bam)])
    parser = shlex.join([
        sys.executable, "workflow/sam_qnames.py",
        "--output-dir", str(outdir), "--tag", tag,
        "--thresholds", *map(str, thresholds),
        "--summary", str(outdir / f"{tag}.mapping_summary.json"),
        "--compression-level", str(compression),
    ])
    shell(sam + " | " + parser, log)


def stream_organelle_mapping(ref: Path, fqs: list[Path], tag: str, run: str,
                              threads: int, thresholds: list[int], outdir: Path,
                              compression: int, log: Path) -> None:
    rg = f"@RG\tID:{run}.{tag}\tSM:{run}\tPL:ILLUMINA"
    bwa = shlex.join(["bwa-mem2", "mem", "-t", str(threads), "-R", rg, str(ref), *map(str, fqs)])
    parser = shlex.join([
        sys.executable, "workflow/sam_qnames.py",
        "--output-dir", str(outdir), "--tag", tag,
        "--thresholds", *map(str, thresholds),
        "--summary", str(outdir / f"{tag}.mapping_summary.json"),
        "--compression-level", str(compression),
    ])
    shell(bwa + " | " + parser, log)


def call_variants(bam: Path, ref: Path, out: Path, cfg: dict, threads: int, log: Path) -> None:
    v = cfg["variant_calling"]
    q = v["minimum_mapping_quality"]
    bq = v["minimum_base_quality"]
    qual = v["minimum_variant_quality"]
    dp = v["minimum_depth"]
    md = v["maximum_depth"]
    type_pipe = "" if v.get("include_indels", True) else f" | bcftools view --threads {threads} -v snps -Ou"
    ploidy = int(v.get("ploidy", 2))
    cmd = (
        f'bcftools mpileup --threads {threads} -f {shlex.quote(str(ref))} -q {q} -Q {bq} -d {md} -Ou {shlex.quote(str(bam))} '
        f'| bcftools call --threads {threads} -m -v --ploidy {ploidy} -Ou '
        f'| bcftools norm --threads {threads} -f {shlex.quote(str(ref))} -Ou'
        f'{type_pipe} '
        f'| bcftools filter --threads {threads} -i {shlex.quote(f"QUAL>={qual} && INFO/DP>={dp}")} '
        f'-Oz -o {shlex.quote(str(out))}'
    )
    shell(cmd, log)
    run_cmd(["bcftools", "index", "--threads", str(threads), "-f", str(out)], log=log)


def depth_metrics(bam: Path, cfg: dict, threads: int) -> dict[str, float | int]:
    v = cfg["variant_calling"]
    min_dp = int(v["minimum_depth"])
    bq = int(v["minimum_base_quality"])
    mq = int(v["minimum_mapping_quality"])
    p = subprocess.Popen(
        ["samtools", "depth", "-@", str(threads), "-q", str(bq), "-Q", str(mq), str(bam)],
        stdout=subprocess.PIPE, text=True,
    )
    covered = callable_pos = depth_sum = 0
    max_depth = 0
    assert p.stdout
    for line in p.stdout:
        d = int(line.rstrip().split("\t")[2])
        covered += 1
        depth_sum += d
        max_depth = max(max_depth, d)
        callable_pos += d >= min_dp
    rc = p.wait()
    if rc:
        raise RuntimeError(f"samtools depth failed: {bam}")
    return {
        "covered_positions": covered,
        "callable_positions_ge_min_depth": callable_pos,
        "mean_depth_covered": depth_sum / covered if covered else 0.0,
        "maximum_depth_observed": max_depth,
    }


def variant_metrics(vcf: Path) -> dict:
    p = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP[\t%GT]\n", str(vcf)],
        stdout=subprocess.PIPE, text=True, check=True,
    )
    sites: set[str] = set()
    n_snp = n_indel = n_het = n_hom_alt = n_missing_gt = 0
    quals: list[float] = []
    dps: list[float] = []
    for line in p.stdout.splitlines():
        parts = line.split("\t")
        chrom, pos, ref, alt, qual, dp = parts[:6]
        gt = parts[6] if len(parts) > 6 else "."
        sites.add(f"{chrom}:{pos}:{ref}:{alt}")
        alts = alt.split(",")
        if len(ref) == 1 and all(len(a) == 1 for a in alts):
            n_snp += 1
        else:
            n_indel += 1
        try:
            quals.append(float(qual))
        except ValueError:
            pass
        try:
            dps.append(float(dp))
        except ValueError:
            pass
        alleles = gt.replace("|", "/").split("/")
        if gt in {".", "./.", ".|."} or any(a == "." for a in alleles):
            n_missing_gt += 1
        elif len(set(alleles)) > 1:
            n_het += 1
        elif alleles and alleles[0] not in {"0", "."}:
            n_hom_alt += 1
    return {
        "variant_sites": sites,
        "n_variants": len(sites),
        "n_snps": n_snp,
        "n_indels": n_indel,
        "n_heterozygous": n_het,
        "n_homozygous_alternate": n_hom_alt,
        "n_missing_genotype": n_missing_gt,
        "mean_variant_qual": sum(quals) / len(quals) if quals else 0.0,
        "mean_variant_dp": sum(dps) / len(dps) if dps else 0.0,
    }


def hardlink_or_copy(src: Path, dst: Path) -> None:
    dst.unlink(missing_ok=True)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    sid = species_id(cfg)
    rd = workdir(cfg) / "downstream" / args.run
    tmp_dir = scratch_dir(cfg, "downstream", args.run)
    require_free_space(tmp_dir, float(resource_value(cfg, "minimum_free_tmp_gb", None, 20)))
    done = rd / ".done"
    log = rd / "downstream.log"
    status = rd / "status.json"
    if done.exists() and not args.force:
        print(f"[SKIP] downstream complete: {args.run}")
        return
    for d in ["sra", "raw", "paired", "singletons", "trimmed", "tmp", "bam", "sets", "conditions", "vcf", "fastp"]:
        (rd / d).mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(workdir(cfg) / "manifest/manifest.evaluable.tsv", sep="\t")
    m = manifest.loc[manifest.Run == args.run]
    if m.empty:
        raise SystemExit(f"[ERROR] selected run absent from manifest: {args.run}")
    row = m.iloc[0]
    layout = str(row.LibraryLayout).upper()
    threads = int(resource_value(cfg, "downstream_threads_per_run", "threads_per_run", 12))
    sort_threads = int(resource_value(cfg, "downstream_sort_threads", "sort_threads", 4))
    sort_mem = str(resource_value(cfg, "downstream_sort_mem", "sort_mem", "2G"))
    sort_level = int(resource_value(cfg, "downstream_sort_compression_level", None, 1))
    bcftools_threads = int(resource_value(cfg, "bcftools_threads", None, 2))
    fasterq_threads = int(resource_value(cfg, "fasterq_threads", "download_threads", 4))
    prefetch_slots = int(resource_value(cfg, "prefetch_slots", None, 2))
    fasterq_slots = int(resource_value(cfg, "fasterq_slots", None, 2))
    repair_mem = int(resource_value(cfg, "repair_memory_gb", None, 8))
    compression = int(cfg["preprocess"].get("compression_level", 1))
    primary = int(cfg["workflow"]["primary_mapq"])

    started = time.monotonic()
    st = {
        "Run": args.run, "species": sid, "status": "running", "stage": "input_resolution",
        "started_utc": utc_now(),
    }
    atomic_json(status, st)
    try:
        local_assignment = resolve_local_read_inputs(cfg, args.run, layout)
        read_policy = local_source_policy(cfg, "reads")
        source_r1: Path | None = None
        source_r2: Path | None = None
        source_single: Path | None = None
        if local_assignment is not None:
            st["stage"] = "local_fastq_reuse"
            st["read_source"] = "local_fastq"
            st["local_r1"] = local_assignment.get("r1", "")
            st["local_r2"] = local_assignment.get("r2", "")
            st["local_single"] = local_assignment.get("single", "")
            atomic_json(status, st)
            if layout == "PAIRED":
                source_r1 = Path(str(local_assignment["r1"]))
                source_r2 = Path(str(local_assignment["r2"]))
            else:
                source_single = Path(str(local_assignment["single"]))
            print(f"[INFO] reusing local FASTQ for downstream {args.run}")
        else:
            allow_download = bool(cfg.get("network", {}).get("allow_read_download", True))
            if read_policy == "local_required" or not allow_download:
                raise FileNotFoundError(
                    f"validated local FASTQ assignment missing for {args.run}; "
                    f"run workflow/local_inputs.py --config {args.config} and inspect "
                    f"{workdir(cfg) / 'local_inputs' / 'unresolved_runs.tsv'}"
                )
            st["stage"] = "prefetch"
            st["read_source"] = "SRA_download"
            atomic_json(status, st)
            with file_slot(cfg, "prefetch", prefetch_slots):
                run_cmd_retry([
                    "prefetch", args.run, "-O", str(rd / "sra"),
                    "--max-size", str(resource_value(cfg, "prefetch_max_size", None, "200G")),
                ], attempts=int(cfg["network"]["retries"]),
                    sleep_seconds=int(cfg["network"]["retry_sleep_seconds"]), log=log)
            sra = find_sra(rd / "sra", args.run)

            st["stage"] = "fasterq_dump"
            atomic_json(status, st)
            with file_slot(cfg, "fasterq", fasterq_slots):
                run_cmd_retry([
                    "fasterq-dump", str(sra), "--split-3", "-e", str(fasterq_threads),
                    "-t", str(tmp_dir), "-O", str(rd / "raw"), "--force",
                ], attempts=2, sleep_seconds=10, log=log)
            if layout == "PAIRED":
                source_r1 = rd / "raw" / f"{args.run}_1.fastq"
                source_r2 = rd / "raw" / f"{args.run}_2.fastq"
            else:
                candidates = [rd / "raw" / f"{args.run}.fastq", rd / "raw" / f"{args.run}_1.fastq"]
                source_single = next((p for p in candidates if p.exists() and p.stat().st_size > 0), None)

        st["stage"] = "repair_preprocess"
        atomic_json(status, st)
        pp = cfg["preprocess"]
        fj = rd / "fastp" / f"{args.run}.json"
        fh = rd / "fastp" / f"{args.run}.html"
        if layout == "PAIRED":
            r1 = source_r1
            r2 = source_r2
            if r1 is None or r2 is None or not (r1.exists() and r2.exists() and r1.stat().st_size > 0 and r2.stat().st_size > 0):
                raise FileNotFoundError(f"paired FASTQ inputs missing for {args.run}: {r1}, {r2}")
            p1 = rd / "paired" / f"{args.run}_1.fastq.gz"
            p2 = rd / "paired" / f"{args.run}_2.fastq.gz"
            sing = rd / "singletons" / f"{args.run}.fastq.gz"
            run_cmd([
                "repair.sh", f"-Xmx{repair_mem}g", f"in1={r1}", f"in2={r2}",
                f"out1={p1}", f"out2={p2}", f"outs={sing}",
                "repair=t", "overwrite=t", "zl=1",
            ], log=log)
            t1 = rd / "trimmed" / f"{args.run}_1.trim.fastq.gz"
            t2 = rd / "trimmed" / f"{args.run}_2.trim.fastq.gz"
            cmd = [
                "fastp", "--thread", str(threads), "--in1", str(p1), "--in2", str(p2),
                "--out1", str(t1), "--out2", str(t2), "--compression", str(compression),
                "--qualified_quality_phred", str(pp["qualified_quality_phred"]),
                "--length_required", str(pp["length_required"]),
                "--unqualified_percent_limit", str(pp["unqualified_percent_limit"]),
                "--n_base_limit", str(pp["n_base_limit"]),
                "--cut_right_window_size", str(pp["cut_right_window_size"]),
                "--cut_right_mean_quality", str(pp["cut_right_mean_quality"]),
                "--json", str(fj), "--html", str(fh),
            ]
            if pp.get("detect_adapter_for_pe", True):
                cmd.append("--detect_adapter_for_pe")
            if pp.get("trim_poly_g", True):
                cmd.append("--trim_poly_g")
            if pp.get("cut_right", True):
                cmd.append("--cut_right")
            run_cmd(cmd, log=log)
            run_cmd([sys.executable, "workflow/pair_check.py", "--r1", str(t1), "--r2", str(t2)], log=log)
            fqs = [t1, t2]
        else:
            raw = source_single
            if raw is None or not raw.exists() or raw.stat().st_size <= 0:
                raise FileNotFoundError(f"single-end FASTQ input missing for {args.run}: {raw}")
            t1 = rd / "trimmed" / f"{args.run}.trim.fastq.gz"
            cmd = [
                "fastp", "--thread", str(threads), "--in1", str(raw), "--out1", str(t1),
                "--compression", str(compression),
                "--qualified_quality_phred", str(pp["qualified_quality_phred"]),
                "--length_required", str(pp["length_required"]),
                "--unqualified_percent_limit", str(pp["unqualified_percent_limit"]),
                "--n_base_limit", str(pp["n_base_limit"]),
                "--cut_right_window_size", str(pp["cut_right_window_size"]),
                "--cut_right_mean_quality", str(pp["cut_right_mean_quality"]),
                "--json", str(fj), "--html", str(fh),
            ]
            if pp.get("trim_poly_g", True):
                cmd.append("--trim_poly_g")
            if pp.get("cut_right", True):
                cmd.append("--cut_right")
            run_cmd(cmd, log=log)
            fqs = [t1]

        st["stage"] = "mapping"
        atomic_json(status, st)
        prep = workdir(cfg) / "references" / "prepared"
        organelle_refs = {c: prep / f"{c}.fa" for c in cfg["references"]["organelles"]}
        nuclear_refs = {
            "nuclear_contig_filtered": prep / "nuclear.contig_filtered.fa",
            "nuclear_masked_main": prep / "nuclear.masked.main.fa",
            "nuclear_masked_strict": prep / "nuclear.masked.strict.fa",
        }
        thresholds = [0, primary]

        # Organelle references are needed only for qname compatibility sets, so
        # avoid writing/sorting BAMs and stream SAM directly into one-pass sets.
        for tag, ref in organelle_refs.items():
            stream_organelle_mapping(ref, fqs, tag, args.run, threads, thresholds,
                                     rd / "sets", compression, log)

        # Nuclear references require coordinate-sorted BAMs for depth and variants.
        for tag, ref in nuclear_refs.items():
            bam = rd / "bam" / f"{tag}.bam"
            rg = f"@RG\tID:{args.run}.{tag}\tSM:{args.run}\tPL:ILLUMINA"
            bwa = shlex.join(["bwa-mem2", "mem", "-t", str(threads), "-R", rg, str(ref), *map(str, fqs)])
            sort = shlex.join([
                "samtools", "sort", "-@", str(sort_threads), "-m", sort_mem,
                "-l", str(sort_level), "-T", str(tmp_dir / tag), "-o", str(bam), "-",
            ])
            shell(bwa + " | " + sort, log)
            run_cmd(["samtools", "index", "-@", str(sort_threads), str(bam)], log=log)
            run_cmd(["samtools", "quickcheck", "-v", str(bam)], log=log)
            qnames_from_bam(bam, rd / "sets", tag, thresholds, sort_threads, compression, log)

        organelle_primary: set[str] = set()
        organelle_any: set[str] = set()
        for c in cfg["references"]["organelles"]:
            organelle_primary |= read_qname_set(rd / "sets" / f"{c}.q{primary}.qnames.txt.gz")
            organelle_any |= read_qname_set(rd / "sets" / f"{c}.q0.qnames.txt.gz")
        nuclear_main = read_qname_set(rd / "sets" / f"nuclear_masked_main.q{primary}.qnames.txt.gz")
        multi = organelle_primary & nuclear_main
        write_names(rd / "sets" / "organelle.primary.txt", organelle_primary)
        write_names(rd / "sets" / "organelle.any_mapq.txt", organelle_any)
        write_names(rd / "sets" / "multi.primary.txt", multi)

        st["stage"] = "reference_sensitivity"
        atomic_json(status, st)
        sensitivity: list[dict] = []
        reference_site_sets: dict[str, set[str]] = {}
        reference_depth: dict[str, dict] = {}
        reference_variant_no_sites: dict[str, dict] = {}
        for tag, ref_fa in nuclear_refs.items():
            n0 = read_qname_set(rd / "sets" / f"{tag}.q0.qnames.txt.gz")
            nq = read_qname_set(rd / "sets" / f"{tag}.q{primary}.qnames.txt.gz")
            ref_bam = rd / "bam" / f"{tag}.bam"
            dm = depth_metrics(ref_bam, cfg, sort_threads)
            ref_vcf = rd / "vcf" / f"reference_{tag}.vcf.gz"
            call_variants(ref_bam, ref_fa, ref_vcf, cfg, bcftools_threads, log)
            with (rd / "vcf" / f"reference_{tag}.bcftools.stats.txt").open("w", encoding="utf-8") as out_stats:
                subprocess.run(["bcftools", "stats", str(ref_vcf)], check=True, stdout=out_stats, text=True)
            vm = variant_metrics(ref_vcf)
            reference_site_sets[tag] = vm.pop("variant_sites")
            reference_depth[tag] = dict(dm)
            reference_variant_no_sites[tag] = dict(vm)
            rec = {
                "Run": args.run, "nuclear_reference": tag,
                "nuclear_qnames_q0": len(n0), "nuclear_qnames_primary": len(nq),
                "organelle_qnames_primary": len(organelle_primary),
                "organelle_qnames_any_mapq": len(organelle_any),
                "multi_qnames_primary": len(organelle_primary & nq),
                "strict_organelle_qnames_primary": len(organelle_primary - nq),
            }
            rec.update(dm)
            rec.update(vm)
            sensitivity.append(rec)

        reference_baseline = reference_site_sets["nuclear_masked_main"]
        for rec in sensitivity:
            ss = reference_site_sets[rec["nuclear_reference"]]
            rec["variant_jaccard_vs_masked_main"] = len(ss & reference_baseline) / len(ss | reference_baseline) if ss | reference_baseline else 1.0
            rec["variant_sites_lost_vs_masked_main"] = len(reference_baseline - ss)
            rec["variant_sites_gained_vs_masked_main"] = len(ss - reference_baseline)
        pd.DataFrame(sensitivity).to_csv(rd / "reference_sensitivity.tsv", sep="\t", index=False)

        st["stage"] = "condition_analysis"
        atomic_json(status, st)
        conditions = {
            "baseline": set(),
            "remove_multi_primary": multi,
            "remove_all_organelle_any_mapq": organelle_any & nuclear_main,
        }
        metric_rows: list[dict] = []
        variant_rows: list[dict] = []
        site_sets: dict[str, set[str]] = {}
        basebam = rd / "bam" / "nuclear_masked_main.bam"
        ref = prep / "nuclear.masked.main.fa"

        for condition, remove in conditions.items():
            cbam = rd / "conditions" / f"{condition}.bam"
            if condition == "baseline":
                hardlink_or_copy(basebam, cbam)
                hardlink_or_copy(Path(str(basebam) + ".bai"), Path(str(cbam) + ".bai"))
                dm = dict(reference_depth["nuclear_masked_main"])
                vm = dict(reference_variant_no_sites["nuclear_masked_main"])
                site_sets[condition] = set(reference_site_sets["nuclear_masked_main"])
                source_vcf = rd / "vcf" / "reference_nuclear_masked_main.vcf.gz"
                vcf = rd / "vcf" / f"{condition}.vcf.gz"
                hardlink_or_copy(source_vcf, vcf)
                for ext in [".csi", ".tbi"]:
                    src_i = Path(str(source_vcf) + ext)
                    if src_i.exists():
                        hardlink_or_copy(src_i, Path(str(vcf) + ext))
            else:
                keep = nuclear_main - remove
                keepfile = rd / "sets" / f"{condition}.keep.txt"
                write_names(keepfile, keep)
                if keep:
                    run_cmd([
                        "samtools", "view", "-@", str(sort_threads), "-N", str(keepfile),
                        "-b", "-1", "-o", str(cbam), str(basebam),
                    ], log=log)
                else:
                    shell(f'samtools view -H {shlex.quote(str(basebam))} | samtools view -b -1 -o {shlex.quote(str(cbam))} -', log)
                run_cmd(["samtools", "index", "-@", str(sort_threads), str(cbam)], log=log)
                run_cmd(["samtools", "quickcheck", "-v", str(cbam)], log=log)
                dm = depth_metrics(cbam, cfg, sort_threads)
                vcf = rd / "vcf" / f"{condition}.vcf.gz"
                call_variants(cbam, ref, vcf, cfg, bcftools_threads, log)
                vm_full = variant_metrics(vcf)
                site_sets[condition] = vm_full.pop("variant_sites")
                vm = vm_full

            dm.update({
                "Run": args.run, "condition": condition,
                "removed_primary_qnames": len(remove),
                "retained_primary_nuclear_qnames": len(nuclear_main - remove),
                "baseline_primary_nuclear_qnames": len(nuclear_main),
            })
            metric_rows.append(dm)
            vm.update({"Run": args.run, "condition": condition})
            variant_rows.append(vm)
            with (rd / "vcf" / f"{condition}.bcftools.stats.txt").open("w", encoding="utf-8") as out_stats:
                subprocess.run(["bcftools", "stats", str(rd / "vcf" / f"{condition}.vcf.gz")], check=True, stdout=out_stats, text=True)

        for vr in variant_rows:
            s = site_sets[vr["condition"]]
            b = site_sets["baseline"]
            vr["jaccard_vs_baseline"] = len(s & b) / len(s | b) if s | b else 1.0
            vr["sites_lost_vs_baseline"] = len(b - s)
            vr["sites_gained_vs_baseline"] = len(s - b)
        pd.DataFrame(metric_rows).to_csv(rd / "downstream_metrics.tsv", sep="\t", index=False)
        pd.DataFrame(variant_rows).to_csv(rd / "variant_metrics.tsv", sep="\t", index=False)

        st.update({
            "status": "ready", "stage": "complete", "completed_utc": utc_now(),
            "wall_seconds": round(time.monotonic() - started, 3),
            "organelle_mapping_mode": "BAM-less_streaming",
            "baseline_variant_reused_from_reference_sensitivity": True,
        })
        atomic_json(status, st)
        done.write_text("ok\n")
        for d in ["sra", "raw", "paired", "singletons", "tmp"]:
            shutil.rmtree(rd / d, ignore_errors=True)
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print(f"[OK] downstream run complete: {args.run} wall={st['wall_seconds']}s")
    except Exception as exc:
        st.update({
            "status": "failed", "error": str(exc), "traceback": traceback.format_exc(),
            "failed_utc": utc_now(), "wall_seconds": round(time.monotonic() - started, 3),
        })
        atomic_json(status, st)
        raise


if __name__ == "__main__":
    main()
