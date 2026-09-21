# BioProject protocol curation policy

Archive metadata and publication-derived metadata are retained separately. The raw SRA/ENA manifest is never rewritten. The workflow may apply an explicit curated overlay to derived grouping variables when a primary publication contradicts an archive label.

The key example is PRJNA285775: the archive-derived `LibraryStrategy` field is `WGS`, but Qiu et al. (2016; DOI 10.1111/mec.13297) describes an SbfI RAD mapping experiment. The workflow therefore preserves `ArchiveLibraryStrategy=WGS` for provenance and uses `LibraryStrategy=RAD-Seq` for curated scientific comparisons.

Missing tissue, extraction, enzyme, fragment-selection, or PCR details are never copied from another BioProject. For PRJDB16402, the associated article and author-supplied supplementary files did not provide PCR-cycle, restriction-enzyme, extraction, or size-selection details; those fields are reported as unavailable and are not inferred.

PRJNA907022 is handled differently: the paper reports paired-end sequencing whereas the retrieved archive records analyzed here are SINGLE. Because the reason could not be established, no layout override is applied. Paper-reported design and archive-analyzed layout are stored as separate fields, and layout comparisons are interpreted descriptively.
