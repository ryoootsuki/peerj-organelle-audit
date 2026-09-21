# Changelog

## 1.0.0

- Preserves the pinned raw SRA manifest while applying publication-based BioProject curation as a separate overlay.
- Corrects PRJNA285775 from the archive-derived WGS label to SbfI RAD-seq for derived analyses and removes the invalid WGS-versus-RRS contrast.
- Adds publication-curated BioProject protocol metadata with explicit source-limited missing fields.
- Includes MAPQ 0/10/20/30/40/50 sensitivity and the complete 27-combination triage-threshold grid.
- Includes downstream alignment/depth/variant sensitivity, nuclear-mask sensitivity, and a 141-run chicken two-reference validation.
- Implements three-reference, two-reference, and organelle-only behavior, with unavailable strict/multi-reference values represented as NA when no nuclear reference exists.
- Includes compact verified result summaries, color-accessible manuscript figures, offline tests, and public-safe provenance summaries.
- Separates reproducible method specifications from historical run-specific intermediate artifacts that are not redistributed in the compact source release.
