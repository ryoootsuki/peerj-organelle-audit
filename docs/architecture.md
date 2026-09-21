# Architecture and robustness design

## User interface

`./run` invokes `workflow/cli.py`. The CLI owns stage order, species selection, logging boundaries, and safe cleanup. All other modules have one narrowly defined responsibility.

## Efficiency

- Local inputs are indexed once; workers do not recursively rescan multi-terabyte roots.
- Archive auditing does not create coordinate-sorted BAM files.
- A single SAM pass stores the maximum primary MAPQ per normalized query name and writes all configured threshold sets.
- Dynamic largest-first scheduling lets workers request a new run immediately after completion.
- Selected-run downstream analysis alone writes sorted BAMs and VCFs.
- Reference indexes and preparation results are reused when the configuration and local-source signatures are unchanged.

## Robustness

- Exact run manifests and expected counts are validated before analysis.
- Run IDs require a complete SRR/ERR/DRR accession pattern.
- FASTQ resolution is deterministic and records ambiguous candidates.
- Paired reads are repaired before fastp and checked after trimming.
- Shell pipelines use `pipefail` and Python workers write atomic JSON status files.
- Temporary disk free space is checked before run-level work.
- Java heap for BBTools is capped per process.
- CPU oversubscription from BLAS/OpenMP is disabled.
- External local inputs are never deleted.

## Reproducibility

- Input manifests are included in the release.
- Local reference provenance includes paths, sizes, mtimes, and SHA-256 hashes.
- Software versions and runtime resource values are written per species.
- Reference thresholds and variant filters come from YAML.
- Missing public metadata are explicit rather than imputed.
