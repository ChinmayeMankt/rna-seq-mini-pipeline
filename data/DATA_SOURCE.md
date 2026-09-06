# Data source

## Default: airway dexamethasone subset

| Field | Value |
|-------|--------|
| GEO | **GSE52778** |
| SRA study | **SRP033351** |
| Paper / design | Himes et al. — human airway smooth muscle, dexamethasone vs untreated |
| Count table origin | **recount2** gene counts (`counts_gene.tsv.gz` for SRP033351) |
| What is in git | `airway_dex_counts.csv` — **600** Ensembl genes × **6** runs (3 untreated, 3 dex), top expressed genes from the recount matrix |
| Metadata | `airway_sample_metadata.csv` (SRR ↔ condition / cell line / GSM) |

Samples used: `SRR1039508`, `SRR1039516`, `SRR1039520` (untreated / control) and `SRR1039509`, `SRR1039513`, `SRR1039517` (dexamethasone / treated). Labels follow the Bioconductor `airway` package convention (dex vs untreated; albuterol arms omitted). Classic airway vignette sample `SRR1039512` is not present in this recount2 matrix, so the subset is balanced 3 vs 3 without it.

**Honesty:** This is a **trimmed learning table**, not a full reanalysis of GSE52778. Counts are from recount2 (Rail-RNA / recount pipeline), not my own STAR/Salmon quantification. Do not treat the Welch t-test results as publication DE.

## Offline toy

`toy_counts.csv` + `sample_metadata.csv` remain available via `python src/analyze.py --dataset toy`.

## Note on count magnitude

recount2 gene-level values for this study are **Rail-RNA coverage-derived** and often have very large per-sample totals compared with featureCounts/Salmon integer read counts. Log-CPM still puts samples on a comparable scale for this learning PCA/t-test demo. For publication-style DE, use the Bioconductor `airway` SummarizedExperiment or re-quantify from FASTQ — not this trimmed CSV alone.
