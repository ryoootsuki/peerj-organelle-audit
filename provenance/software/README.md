# Software/environment provenance

`software_versions_silene_recorded_public.tsv` and `software_versions_chicken_recorded_public.tsv` are public-safe copies of version information recovered from retained completed-analysis outputs. Machine-specific path strings were redacted; see `SANITIZATION_NOTE.md`.

The retained records support, among other items, Python 3.12.13, fastp 1.3.6, SeqKit 2.13.0, minimap2 2.31-r1302, bedtools 2.31.1, samtools 1.22.1, and bcftools 1.22.1. The historical `repair.sh` line does not provide a clean BBTools package version, so no exact BBTools version is inferred.

`environment.yml` and `environment.flexible.yml` are reproducible environment specifications for new executions. They are not presented as byte-identical exports of the historical solver state.
