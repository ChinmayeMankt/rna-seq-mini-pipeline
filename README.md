# rna-seq-mini-pipeline

**Learning project** (in progress). This is a small, fully local walkthrough of a *count-table* RNA-seq analysis: QC, library-size normalization, and a simple differential-expression *style* comparison with pandas and scipy.

It is **not** a production NGS pipeline. It does **not** start from FASTQ, does **not** run STAR/Salmon/HISAT2, and does **not** replace DESeq2/edgeR. I am Chinmaye Mankatalia, a B.Tech Biotechnology student, working through these concepts so I can talk honestly about RNA-seq internships.

## What it does

1. Loads a **tiny toy count matrix** checked into `data/toy_counts.csv` (40 genes × 8 samples: 4 control, 4 treated).
2. Reports per-sample library sizes, detects genes with all-zero counts, and plots a PCA of log-CPM values.
3. Computes log2 fold-change (treated vs control) and a Welch t-test on log-CPM, with Benjamini–Hochberg FDR.

The toy matrix is synthetic. A handful of genes (`GENE001`–`GENE003` up, `GENE011`–`GENE013` down) were planted as “DE-like” so the script has something to recover. Treat results as a **methods demo**, not biology.

## How to run

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
python src/analyze.py
```

Outputs land in `results/` (`de_results.csv`, `qc_metrics.csv`, optional PNG plots).

Optional notebook (same workflow, more commentary):

```bash
jupyter notebook notebooks/rna_seq_mini.ipynb
```

## Data

| File | Role |
|------|------|
| `data/toy_counts.csv` | Gene × sample integer counts (synthetic) |
| `data/sample_metadata.csv` | Sample → condition |

If you later swap in a real GEO count table, keep it small and cite the accession in this README. Do not commit FASTQs.

## Stack

Python, pandas, numpy, scipy, matplotlib. No Bioconductor / scanpy here on purpose.

See [LEARNING.md](LEARNING.md) for what I still do not know.
