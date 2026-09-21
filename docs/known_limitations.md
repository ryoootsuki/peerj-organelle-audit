# Known limitations

1. **Compatibility is not origin assignment.** Independent mapping identifies reference compatibility. It does not prove the cellular or evolutionary origin of each read and is not a locus-level NUMT/NUPT discovery method.
2. **Downstream analysis is a sensitivity test.** Without an external truth set, changed variant calls cannot be labeled more accurate merely because organelle-compatible reads were removed.
3. **Public datasets are confounded.** BioProject, tissue, extraction, enzyme, fragment selection, layout, instrument, depth, and archive representation are not independently randomized. Group differences are descriptive.
4. **No verified WGS comparator remains in the curated Silene panel.** PRJNA285775 is SbfI RAD-seq despite the archive-derived `WGS` field. The invalid WGS–RRS comparison was removed.
5. **Protocol metadata are source-limited.** PRJDB16402 and PRJEB34584 retain unresolved fields rather than borrowing values from related studies.
6. **PRJNA907022 paper/archive layout differs.** Paper-reported paired-end design and analyzed SINGLE archive records are retained as separate fields; the cause is not inferred.
7. **Nuclear masking is itself an analytical choice.** Main and high-confidence masks are sensitivity settings, not validated complete catalogs of organelle-derived nuclear inserts.
8. **Organelle-only mode is necessarily limited.** Inclusive compatibility can be reported, but strict and nuclear–organelle ambiguity metrics are undefined without a nuclear reference and remain NA.
9. **Historical byte-level execution provenance is incomplete in the compact retained result bundle.** The exact generated mask-report files and a byte-identical solved Conda lock were not retained. This source release therefore does not claim that those historical generated files are included; it provides the code, parameters, manifests, compact verified outputs, and documented provenance that are available.
10. **Analysis and source-release identifiers are distinct.** The completed result bundle recorded an internal executed workflow identifier (`2.0.0-lean-epyc7713-local`), whereas this public source tree is versioned `1.0.0`. The scientific outputs in `examples/results_summary/` are retained results from the completed analysis, and the public source tree is the cleaned, portable implementation distributed with them.
