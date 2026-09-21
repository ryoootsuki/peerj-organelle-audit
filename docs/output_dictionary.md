# Principal outputs

## Per species

- `work/<species>/local_inputs/read_resolution.tsv`: deterministic accession-to-local-FASTQ assignments used by all workers.
- `work/<species>/local_inputs/read_inventory_summary.tsv`: resolved/missing/invalid counts and scanned roots.
- `work/<species>/local_inputs/unresolved_runs.tsv`: runs that must be fixed before local-required execution.
- `work/<species>/local_inputs/read_candidates.tsv`: all raw-looking local candidates considered by the resolver.
- `work/<species>/references/reports/reference_source_provenance.tsv`: local reference source paths, sizes, mtimes and SHA-256 values.
- `results/<species>/provenance/software_versions.tsv`: executable versions and input hashes.
- `results/<species>/provenance/run_completion_status.tsv`: success/failure state for every pinned run.
- `results/<species>/provenance/analysis_completion_summary.tsv`: expected versus observed run counts.
- `work/<species>/references/reports/`: removed contigs, BED masks, SeqKit statistics, reference hashes and thresholds.
- `results/<species>/tables/fastp_run_metrics.tsv`: before/after fastp metrics.
- `work/<species>/runs/<Run>/sets/*.mapping_summary.json`: BAM-less audit mapping summaries produced by the one-pass SAM parser.
- `results/combined/provenance/performance_summary.tsv`: observed run-level wall-time and throughput summary.
- `results/combined/provenance/run_wall_times.tsv`: per-run audit and downstream wall times.
- `results/<species>/tables/audit_mapq_sweep.tsv`: all completed runs at MAPQ 0, 10, 20, 30, 40 and 50.
- `results/<species>/tables/audit_primary_mapq.tsv`: main MAPQ 20 table.
- `results/<species>/tables/mapq_distribution_summary.tsv`: median, quartiles, IQR, P95 and maximum; provides numerical distribution summaries corresponding to the MAPQ figure.
- `results/<species>/tables/grouped_by_*.tsv`: project, model, layout, strategy and platform summaries.
- `results/<species>/tables/library_protocol_metadata.tsv`: run-level archive metadata merged with publication-curated BioProject protocol fields; unresolved values remain explicit.
- `results/<species>/tables/technical_group_comparisons.tsv`: applicable descriptive technical-group comparisons with effect sizes. A WGS/RRS row is produced only if a genuinely curated WGS group exists; it is absent for the current Silene benchmark.
- `results/<species>/tables/downstream_subset.tsv`: prespecified run selection and rationale.
- `results/<species>/tables/downstream_alignment_depth_deltas.tsv`: alignment/depth changes after filtering.
- `results/<species>/tables/downstream_variant_deltas.tsv`: variant count, gain/loss and Jaccard changes.
- `results/<species>/tables/reference_mask_sensitivity.tsv`: mapping, depth and variant outcomes for contig-filtered and two masked nuclear references.
- `results/<species>/tables/reference_mask_variant_deltas.tsv`: depth/variant changes and Jaccard overlap relative to the primary mask.
- `results/<species>/figures/`: manuscript and supplementary figures.

## Combined

- `results/combined/species_validation_summary.tsv`: plant and animal validation summary.
- `results/combined/report_materials/`: generated methods/results summaries and a computational-output checklist.



## Compact retained summaries

- `examples/results_summary/`: verified compact summaries from the completed analysis, including MAPQ distributions, threshold sensitivity, downstream sensitivity, reference-mask sensitivity, and species validation.
