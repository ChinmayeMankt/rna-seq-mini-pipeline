# LEARNING.md

I (Chinmaye Mankatalia) am using this repository to **practice**, not to claim RNA-seq or production NGS expertise.

## What this repo is meant to show

- Comfort with pandas on a gene×sample count table from a **documented public accession**
- Library-size / CPM thinking and a simple log transform
- That “differential expression” in production is a statistical model (negative binomial, empirical Bayes), not just a t-test
- How to cite GEO/SRA/recount when swapping toy data for a real table

## What this repo is not

- Not alignment, quantification, or FASTQ QC (FastQC, MultiQC, STAR, salmon)
- Not DESeq2, edgeR, or limma-voom
- Not a clinical or publication-grade analysis of GSE52778
- Not an internship, course credit, or collaboration with a lab

## Public data used

- **GSE52778** / **SRP033351** (airway smooth muscle, dex vs untreated), counts from **recount2**, trimmed to 150 genes × 6 samples for a git-friendly demo. See `data/DATA_SOURCE.md`.

## Next steps I want to take

1. Repeat the same table in R with DESeq2 on the full airway object
2. Read a short nf-core/rnaseq overview so I know where this notebook sits in a real workflow
3. Learn why batch, outliers, and gene-length matter before claiming I “do RNA-seq”
