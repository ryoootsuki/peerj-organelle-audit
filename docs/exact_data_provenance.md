# Exact data provenance

## *Silene latifolia*

- Raw archive metadata: `manifests/silene_SraRunInfo_retrieved_2025-12-25.csv`
- Retrieval date recorded for the submitted benchmark: 25 December 2025
- Total records: 2,686
- Evaluable records (non-zero SRA spot count): 2,679
- Nuclear assembly accession: GCF_048544455.1
- Chloroplast accession: NC_016730.1
- Selected mitochondrial reference: HM562727.1

The raw manifest is preserved unchanged. It contains 26 PRJNA285775 records with the archive-derived `LibraryStrategy=WGS` field. Publication-level curation (Qiu et al. 2016, DOI 10.1111/mec.13297) identifies these as SbfI RAD-seq. The workflow therefore writes both `ArchiveLibraryStrategy` and the curated `LibraryStrategy`, and all derived strategy comparisons use the curated value. No valid WGS-versus-RRS comparison is claimed for this curated Silene collection.

Public example paths in `config/species/silene.yaml` are placeholders; users must configure their own local FASTQ and reference files.

## *Gallus gallus*

- Exact run manifest: `manifests/chicken_PRJNA573756_exact_141_runs.tsv`
- BioProject: PRJNA573756
- Runs: 141
- Nuclear assembly accession: GCF_016699485.2
- Mitochondrial accession: NC_053523.1

The mitochondrial record can be extracted from the same local assembly FASTA before preparation of the nuclear reference.

## Runtime provenance

At execution, the workflow records resolved input paths, file sizes, modification times, layouts/candidate counts, reference hashes, software versions, completion status and reference-mask reports. Machine-specific paths are runtime provenance and are intentionally not committed to the public repository.
