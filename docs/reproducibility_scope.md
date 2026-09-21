# Reproducibility scope

This source package distinguishes three levels of reproducibility evidence.

1. **Pinned inputs and methods.** Exact archive manifests, reference accessions, configuration files, thresholds, and workflow code are included.
2. **Verified retained summaries.** Compact result tables under `examples/results_summary/` support the numerical claims used in the revised manuscript.
3. **Run-specific intermediates.** Large FASTQ/BAM/VCF/reference files and completed reference-mask report directories are not redistributed in this compact release. The workflow regenerates them and writes provenance records during a new run.

Recorded public-safe software summaries are under `provenance/software/`. The environment YAML is a reproducible specification; it should not be interpreted as a byte-identical solver lock for the historical completed run. The example 512-GB EPYC profile is a planning configuration, whereas `provenance/hardware/hardware_observed_from_results.tsv` records the approximately 1-TB memory detected in the retained completed-run summary.

These distinctions prevent the repository from claiming exact provenance artifacts that are not actually present while preserving a complete, testable workflow and the retained evidence used for the reported sensitivity results.
