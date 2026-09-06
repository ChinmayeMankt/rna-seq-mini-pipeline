# rna-seq-mini-pipeline

**Learning project** (in progress). A small, fully local walkthrough of a *count-table* RNA-seq analysis: QC, library-size normalization, and a simple differential-expression *style* comparison with pandas and scipy.

**Default input is a real public count table** (airway / **GSE52778** / **SRP033351**, recount2 gene counts, small subset). A synthetic toy matrix remains as `--dataset toy`.

It is **not** a production NGS pipeline. It does **not** start from FASTQ, does **not** run STAR/Salmon/HISAT2, and does **not** replace DESeq2/edgeR. I am Chinmaye Mankatalia, a B.Tech Biotechnology student, working through these concepts so I can talk honestly about RNA-seq internships.

## What it does

1. Loads `data/airway_dex_counts.csv` (150 genes × 6 samples: 3 untreated, 3 dexamethasone) plus metadata.
2. Reports per-sample library sizes, detects genes with all-zero counts, and plots a PCA of log-CPM values.
3. Computes log2 fold-change (treated vs control) and a Welch t-test on log-CPM, with Benjamini–Hochberg FDR.

Treat statistical hits as a **methods demo**. Real DE for this design is usually modeled with DESeq2/edgeR on the full transcriptome.

## How to run

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
python src/analyze.py                 # real airway subset (default)
python src/analyze.py --dataset toy   # synthetic fallback
```

Outputs land in `results/` (`de_results.csv`, `qc_metrics.csv`, `pca_samples.png`).

Optional notebook:

```bash
jupyter notebook notebooks/rna_seq_mini.ipynb
```

## Data

| File | Role |
|------|------|
| `data/airway_dex_counts.csv` | Real public subset — GSE52778 / SRP033351 via recount2 |
| `data/airway_sample_metadata.csv` | SRR → condition (control/treated), cell line, GSM |
| `data/DATA_SOURCE.md` | Accession + trimming notes |
| `data/toy_counts.csv` | Synthetic (optional) |
| `data/sample_metadata.csv` | Toy sample → condition |

Do not commit FASTQs.

## Stack

Python, pandas, numpy, scipy, matplotlib. No Bioconductor / scanpy here on purpose.

See [LEARNING.md](LEARNING.md) for what I still do not know.
