# Organelle compatibility audit workflow for public sequencing archives

This repository provides a reproducible workflow for auditing organelle-compatible reads and cross-reference ambiguity in public short-read datasets. The benchmark implementation uses *Silene latifolia* in three-reference mode (nuclear, chloroplast, mitochondrial) and *Gallus gallus* in nuclear–mitochondrial two-reference mode.

## What the workflow reports

The workflow independently maps reads to the available references and summarizes three complementary quantities across MAPQ 0, 10, 20, 30, 40, and 50: inclusive organelle compatibility, strict organelle compatibility after excluding nuclear-compatible query names, and multi-reference ambiguity. MAPQ 20 is the primary reporting threshold, but the full sweep is retained. Operational triage rules are provided only as auditable screening heuristics; the included 27-combination sensitivity grid demonstrates that categorical assignments depend on the selected cutoffs. Continuous fractions remain the primary evidence.

## Benchmark data

- *Silene latifolia*: pinned SRA snapshot of 2,686 records, with 2,679 evaluable runs. Publication-based curation treats all records as RAD-seq or related reduced-representation data. PRJNA285775 is SbfI RAD-seq; the raw archive label is preserved in the pinned manifest and corrected only in derived metadata.
- *Gallus gallus*: 141 RAD-seq runs from PRJNA573756, analyzed in two-reference mode.

Curated protocol metadata are in `metadata/bioproject_protocol_metadata_curated.tsv`. Source-unavailable fields are left unresolved rather than inferred.

## Verified compact results

`examples/results_summary/` contains lightweight summaries from the completed analysis, including the MAPQ distributions, downstream variant sensitivity, nuclear-mask sensitivity, the full 27-row triage-threshold grid, and plant/animal validation summary. At the prespecified 5%/20%/1% triage setting, 208/2,679 Silene runs change category; across the full threshold grid, the count ranges from 20 to 1,146.

## Installation and network-free tests

```bash
mamba env create -f environment.yml
conda activate peerj-organelle-audit
make test
```

The tests validate the pinned manifests, the PRJNA285775 metadata correction, local FASTQ resolution, PAF masking logic, query-name extraction, one-reference behavior, key numerical results, and public-package safety checks.

## Configure local inputs

Edit only the `local_inputs` sections of `config/species/silene.yaml` and `config/species/chicken.yaml` to point to local FASTQ/reference locations. Public configuration files contain no machine-specific analysis paths.

## Typical execution

```bash
make plan
make autotune
make doctor
make local-check
make references
make pilot
make full
make readiness
```

The pipeline is resume-safe at the run level. Large raw FASTQ, BAM, VCF, reference FASTA, and index files are intentionally excluded from this source repository.

## Primary parameters

- fastp qualified Phred threshold: 20
- minimum read length: 30 bp
- primary MAPQ: 20
- MAPQ sweep: 0, 10, 20, 30, 40, 50
- primary nuclear mask: alignment length ≥200 bp and identity ≥90%
- high-confidence mask sensitivity setting: alignment length ≥500 bp and identity ≥95%
- primary triage heuristics: low 5%, high 20%, ambiguity 1%
- triage sensitivity grid: low 2.5/5/10%, high 15/20/25%, ambiguity 0.5/1/2%

## Reference modes

- three-reference mode: nuclear + chloroplast + mitochondrial
- two-reference mode: nuclear + one organelle reference
- organelle-only mode: inclusive compatibility is available; strict and nuclear–organelle ambiguity are reported as unavailable rather than zero

See `docs/reference_modes.md` and `docs/reproducibility_scope.md`.

## Reproducibility scope

This source release contains fixed manifests, configuration, workflow code, curated metadata, compact verified result summaries, and recorded software/hardware summaries recoverable from retained outputs. It does not redistribute large sequence/alignment intermediates or claim byte-identical recovery of run-specific intermediate mask-report files that were not retained in the compact result bundle. New executions generate reference-preparation reports and hashes under `work/<species>/references/reports/`; users should archive those outputs for exact run provenance.

## License

MIT License applies to original workflow code in this repository. Third-party sequence data, reference assemblies, and publications retain their own terms.
