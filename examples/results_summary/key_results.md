# Verified summary of the completed analysis run

These files are small, compact verified summaries extracted from the completed analysis. They are included for transparency; the workflow remains the authoritative way to regenerate outputs.

## Silene archive audit
- evaluable runs at primary MAPQ: **2679**
- MAPQ 20 strict fraction: median **2.381%**, 95th percentile **21.248%**, maximum **32.155%**
- MAPQ 20 multi-reference fraction: median **0.458%**
- triage category changed after separating strict signal from ambiguity: **208/2679 (7.76%)**
- publication-based curation means there is **no valid WGS comparator** in this Silene set; PRJNA285775 is SbfI RAD-seq.

## Downstream sensitivity
- Silene selected-run variant comparison: **26 unique runs**, **52 filtered-condition comparisons**
- largest absolute Silene variant-count change: **4644**
- minimum Silene variant-site Jaccard versus baseline: **0.98999**
- chicken selected-run variant comparison: **9 unique runs**
- largest absolute chicken variant-count change: **1**
- minimum chicken variant-site Jaccard versus baseline: **0.99980**

## Interpretation
These are sensitivity results, not truth-set validation. Large or small changes indicate how strongly the selected analysis depends on organelle-compatible/multi-reference read handling under the stated references and thresholds.
