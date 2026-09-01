#!/usr/bin/env python3
"""Toy RNA-seq count-table QC + DE-style comparison.

Learning script only. Not a substitute for DESeq2/edgeR.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


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
    counts = pd.read_csv(DATA / "toy_counts.csv", index_col=0)
    meta = pd.read_csv(DATA / "sample_metadata.csv").set_index("sample")
    counts = counts.loc[:, meta.index]

    qc = pd.DataFrame(
        {
            "library_size": counts.sum(axis=0),
            "detected_genes": (counts > 0).sum(axis=0),
            "condition": meta["condition"],
        }
    )
    qc.to_csv(OUT / "qc_metrics.csv")
    print("QC (per sample)")
    print(qc.to_string())
    print(f"\nAll-zero genes: {(counts.sum(axis=1) == 0).sum()}")

    lcpm = log_cpm(counts)

    X = lcpm.T.to_numpy()
    X = X - X.mean(axis=0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    pcs = X @ vt[:2].T
    fig, ax = plt.subplots(figsize=(5, 4))
    for cond, color in [("control", "#4C78A8"), ("treated", "#F58518")]:
        mask = meta["condition"] == cond
        ax.scatter(pcs[mask, 0], pcs[mask, 1], label=cond, c=color, s=60)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of log-CPM (toy data)")
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
    print("\nTop genes by FDR (toy DE-style t-test on log-CPM):")
    print(res.head(8).to_string())
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
