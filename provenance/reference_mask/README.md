# Reference-mask provenance required for final public release

The completed analysis generated these files under `work/<species>/references/reports/`, but they were not included in the supplied `results.zip`. They must be copied here **without reconstruction or guessing** before the final `v1.0.0` release.

For each of `silene/` and `chicken/`, retain at least:

- `reference_provenance.tsv`
- `reference_source_provenance.tsv`
- `removed_organelle_contigs.ids`
- `whole_contig_alignment_report.tsv`
- `mask.main.merged.bed`
- `mask.strict.merged.bed`
- `mask.main.summary.tsv`
- `mask.strict.summary.tsv`
- `reference_seqkit_stats.tsv`
- `reference_sha256.txt`
- `reference_config.sha256`

These artifacts are required because the manuscript explicitly describes the reference-removal/masking procedure and states that it is auditable. They should be copied from the exact completed analysis, not regenerated with altered references.
