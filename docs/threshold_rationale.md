# Threshold rationale

MAPQ 20 and the 5%/20%/1% triage cutoffs are prespecified reporting/screening rules rather than biological discontinuities. The workflow therefore retains continuous run-level fractions, performs a full MAPQ 0/10/20/30/40/50 sweep, and writes `triage_threshold_sensitivity.tsv` across alternative cutoffs.

This design prevents a run immediately above a threshold from being interpreted as biologically categorically different from a run immediately below it. Triage labels are an operational aid for prioritizing downstream sensitivity analyses; they are not evidence of cellular origin or data quality by themselves.
