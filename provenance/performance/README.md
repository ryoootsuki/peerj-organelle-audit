# Performance provenance

The original performance summary in the supplied result bundle counted stale/pilot downstream runs (33 Silene and 17 chicken), whereas the final scientific downstream subsets contain 26 and 9 unique runs.

`performance_summary_final_subset.tsv` recomputes downstream performance only for accession IDs present in the final `downstream_unique_runs.tsv` files. It does **not** change any scientific results; it only prevents stale/pilot timing records from being presented as final-subset provenance.
