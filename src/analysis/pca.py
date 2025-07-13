#!/usr/bin/env python3
"""
pca.py
──────
Create a side-by-side figure with:

A.  2-D PCA scatter of radiomic features, coloured by hippocampal side
B.  Bar-plot of the top PC1/PC2 loadings (|loading| ≥ 0.1)

The plotting code (Seaborn + Matplotlib) and all atlas/legend tweaks are
_verbatim_ from the original 2D PCA Loadings script so the figure will look
identical.  Default CLI arguments reproduce the manuscript figure out-of-the-box.

Example
-------
python -m rad_pipeline.analysis.pca \
    --csv data/Manually_preprocessed_and_verified/radiomics_output/combined_features_all_mice_hippocampus_V.1.2_binsize_32.csv \
    --out pca_scatter_and_loadings.png
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


# ──────────────────────────────────────────────────────────────────────────
# Core function (logic unchanged)
# ──────────────────────────────────────────────────────────────────────────
def make_pca_figure(csv_path: Path, out_path: Path) -> None:
    # 1. load + pivot
    data = pd.read_csv(csv_path)
    data_pivot = (
        data.pivot(index=["FileName", "Image"], columns="Feature", values="Value")
        .dropna()
        .reset_index()
    )

    # 2. PCA
    features = data_pivot.drop(columns=["FileName", "Image"]).values
    pca = PCA(n_components=2)
    comps = pca.fit_transform(features)
    data_pivot["PCA1"], data_pivot["PCA2"] = comps[:, 0], comps[:, 1]

    # 3. loadings
    feat_cols = data_pivot.drop(columns=["FileName", "Image", "PCA1", "PCA2"]).columns
    feat_cols = [c.replace("original_", "") for c in feat_cols]
    loadings = pd.DataFrame(
        {
            "Feature": feat_cols,
            "PC1 Loading": pca.components_[0],
            "PC2 Loading": pca.components_[1],
        }
    ).set_index("Feature")
    filt = loadings[
        (loadings["PC1 Loading"].abs() >= 0.1) | (loadings["PC2 Loading"].abs() >= 0.1)
    ]

    # 4. plot
    fig, axes = plt.subplots(1, 2, figsize=(22, 12))

    # scatter
    sns.scatterplot(
        x=data_pivot["PCA1"],
        y=data_pivot["PCA2"],
        hue=data_pivot["Image"],
        palette={
            "Hippocampal Region Left": "blue",
            "Hippocampal Region Right": "red",
        },
        s=400,
        edgecolor="k",
        ax=axes[0],
    )
    axes[0].set_title(
        "A. 2D PCA of Radiomic Features by Hippocampal Region",
        fontsize=22,
        fontweight="bold",
    )
    axes[0].set_xlabel("PCA Component 1", fontsize=22, fontweight="bold")
    axes[0].set_ylabel("PCA Component 2", fontsize=22, fontweight="bold")
    axes[0].tick_params(axis="both", which="major", labelsize=18)

    # point labels “M_n”
    for i, txt in enumerate(
        data_pivot["FileName"].str.replace("_Verified", "", regex=True)
    ):
        match = re.search(r"\d+", txt)
        number = match.group() if match else txt
        label = f"M_{number}"
        axes[0].annotate(
            label,
            (data_pivot["PCA1"][i], data_pivot["PCA2"][i]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=16,
            color="black",
            fontweight="bold",
        )

    leg_scatter = axes[0].legend(
        handles=axes[0].get_legend_handles_labels()[0],
        labels=[
            "Hippocampal Region Left",
            "Hippocampal Region Right (irradiated)",
        ],
        title="Hippocampal Region",
        title_fontsize=16,
        fontsize=14,
    )
    for t in leg_scatter.get_texts():
        t.set_fontweight("bold")
    leg_scatter.get_title().set_fontweight("bold")

    # bar-plot
    filt[["PC1 Loading", "PC2 Loading"]].plot(
        kind="bar",
        width=0.75,
        edgecolor="black",
        color=["crimson", "royalblue"],
        ax=axes[1],
    )
    axes[1].set_title(
        "B. Top Feature Contributions to PC1 and PC2",
        fontsize=22,
        fontweight="bold",
    )
    axes[1].set_xlabel("Features", fontsize=18, fontweight="bold")
    axes[1].set_ylabel("Loading", fontsize=18, fontweight="bold")
    axes[1].tick_params(axis="x", rotation=45, labelsize=16)
    axes[1].tick_params(axis="y", labelsize=16)
    for lbl in axes[1].get_xticklabels():
        lbl.set_fontweight("bold")
    leg_bar = axes[1].legend(
        title="Principal Components", title_fontsize=16, fontsize=14, loc="upper right"
    )
    for t in leg_bar.get_texts():
        t.set_fontweight("bold")
    leg_bar.get_title().set_fontweight("bold")

    plt.tight_layout()
    plt.savefig(out_path, transparent=True)
    print(f"✓ figure written to {out_path}")
    plt.show()


# ──────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot 2-D PCA scatter + loadings bar.")
    p.add_argument(
        "--csv",
        default="data/Manually_preprocessed_and_verified/radiomics_output/combined_features_all_mice_hippocampus_V.1.2_binsize_32.csv",
        help="Combined-features CSV (default: manuscript dataset)",
    )
    p.add_argument(
        "--out",
        default="pca_scatter_and_loadings.png",
        help="Output PNG file (default: pca_scatter_and_loadings.png)",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    make_pca_figure(Path(args.csv), Path(args.out))


if __name__ == "__main__":
    main()
