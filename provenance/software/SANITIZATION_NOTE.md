# Path redaction in recorded software provenance

The completed result bundle recorded tool-version command output containing machine-specific absolute paths. For the public source release, **only machine-specific path strings** were replaced with neutral placeholders such as `<HOME>` or `<LOCAL_STORAGE>`.

No software version, analytical parameter, hash, run count, or scientific result was altered. The files are therefore named `*_recorded_public.tsv` rather than `*_raw.tsv` to avoid implying byte-for-byte identity with the private execution record.

The `repair.sh` row is not a reliable BBMap version report: it records the Java launcher command. The package therefore treats the exact BBMap/BBTools version as unresolved until verified from the completed Conda environment (for example, `conda list bbmap`).
