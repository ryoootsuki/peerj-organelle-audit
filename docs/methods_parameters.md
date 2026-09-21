# Methods parameters recorded by the workflow

## Pre-processing
- `fastp --qualified_quality_phred 20`
- `--length_required 30`
- `--unqualified_percent_limit 40`
- `--n_base_limit 5`
- `--cut_right_window_size 4`
- `--cut_right_mean_quality 20`
- `--trim_poly_g --cut_right`
- paired-end libraries: `--detect_adapter_for_pe`
- paired inputs are synchronized with BBTools `repair.sh` before trimming.

## Independent mapping
Reads are mapped independently to each available reference. Query-name sets are generated from primary mapped alignments after the workflow's fixed SAM flag handling and evaluated at MAPQ 0, 10, 20, 30, 40 and 50. MAPQ 20 is the prespecified primary reporting value.

## Nuclear-reference audit
- whole-contig removal: merged organelle alignment coverage ≥80% and weighted identity ≥95% (plus explicit organelle records)
- primary local mask: alignment length ≥200 bp and identity ≥90%
- high-confidence sensitivity mask: alignment length ≥500 bp and identity ≥95%
- masking uses `N`; lowercase soft masking is not assumed to prevent BWA-MEM2 mapping.

## Triage rules
Primary screening: low=5%, high=20%, ambiguity=1%. These are operational categories, not natural biological thresholds. Sensitivity output spans low 2.5/5/10%, high 15/20/25%, ambiguity 0.5/1/2%.

## Variant sensitivity
- mapping quality ≥20
- base quality ≥20
- variant QUAL ≥30
- depth ≥3
- configured ploidy: 2

Variant outputs are sensitivity demonstrations, not truth-set validation.
