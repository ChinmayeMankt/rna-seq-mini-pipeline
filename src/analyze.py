#!/usr/bin/env python3
"""RNA-seq count-table QC + DE-style comparison on a real public subset.

Default data: airway smooth-muscle RNA-seq (GSE52778 / SRP033351), gene counts
from recount2, dexamethasone vs untreated (6 samples × 600 genes subset).

Learning script only. Not a substitute for DESeq2/edgeR. Not FASTQ→counts.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

# Documented public source (see data/DATA_SOURCE.md)
DEFAULT_COUNTS = DATA / "airway_dex_counts.csv"
DEFAULT_META = DATA / "airway_sample_metadata.csv"
TOY_COUNTS = DATA / "toy_counts.csv"
TOY_META = DATA / "sample_metadata.csv"


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR (simple implementation)."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def log_cpm(counts: pd.DataFrame, prior: float = 0.5) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    cpm = (counts + prior).div(lib, axis=1) * 1e6
    return np.log2(cpm)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--dataset",
        choices=("airway", "toy"),
        default="airway",
        help="airway = real GSE52778/SRP033351 recount2 subset; toy = synthetic CSV",
    )
    args = p.parse_args()

    if args.dataset == "airway":
        counts_path, meta_path = DEFAULT_COUNTS, DEFAULT_META
        title = "PCA of log-CPM (airway GSE52778 subset)"
    else:
        counts_path, meta_path = TOY_COUNTS, TOY_META
        title = "PCA of log-CPM (toy data)"

    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path).set_index("sample")
    counts = counts.loc[:, meta.index]

    qc = pd.DataFrame(
        {
            "library_size": counts.sum(axis=0),
            "detected_genes": (counts > 0).sum(axis=0),
            "condition": meta["condition"],
        }
    )
    qc.to_csv(OUT / "qc_metrics.csv")
    print(f"Dataset: {args.dataset}  counts={counts_path.name}")
    print("QC (per sample)")
    print(qc.to_string())
    print(f"\nAll-zero genes: {(counts.sum(axis=1) == 0).sum()}")
    print(f"Matrix shape: {counts.shape[0]} genes × {counts.shape[1]} samples")

    lcpm = log_cpm(counts)

    X = lcpm.T.to_numpy()
    X = X - X.mean(axis=0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    pcs = X @ vt[:2].T
    fig, ax = plt.subplots(figsize=(5, 4))
    for cond, color in [("control", "#4C78A8"), ("treated", "#F58518")]:
        mask = (meta["condition"] == cond).to_numpy()
        ax.scatter(pcs[mask, 0], pcs[mask, 1], label=cond, c=color, s=60)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "pca_samples.png", dpi=120)
    plt.close(fig)

    ctrl = lcpm.loc[:, meta["condition"] == "control"]
    trt = lcpm.loc[:, meta["condition"] == "treated"]
    log2fc = trt.mean(axis=1) - ctrl.mean(axis=1)
    pvals = []
    for gene in lcpm.index:
        t = stats.ttest_ind(trt.loc[gene], ctrl.loc[gene], equal_var=False)
        pvals.append(t.pvalue)
    res = pd.DataFrame({"log2fc": log2fc, "pvalue": pvals})
    res["fdr"] = bh_fdr(res["pvalue"].to_numpy())
    res = res.sort_values("fdr")
    res.to_csv(OUT / "de_results.csv")
    print("\nTop genes by FDR (Welch t-test on log-CPM — learning demo, not DESeq2):")
    print(res.head(8).to_string())
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
