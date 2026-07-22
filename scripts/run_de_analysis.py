"""Reproducible APOEhi vs Control differential expression pipeline.

Generates figures/ and figures/apoehi_vs_control_de_genes.csv without Jupyter.
"""
from __future__ import annotations

import os
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

CONTROL_COLS = [
    "Control-1 (FPKM)",
    "Control-2 (FPKM)",
    "Control-3 (FPKM)",
]
APOEHI_COLS = [
    "APOEhi-1 (FPKM)",
    "APOEhi-2 (FPKM)",
    "APOEhi-3 (FPKM)",
]
APOELO_COLS = [
    "APOElo-1 (FPKM)",
    "APOElo-2 (FPKM)",
    "APOElo-3 (FPKM)",
]
NUMERIC_COLS = CONTROL_COLS + APOEHI_COLS + APOELO_COLS
PSEUDOCOUNT = 0.001


def load_and_clean(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path)
    df.columns = df.iloc[0]
    df = df[1:].reset_index(drop=True)
    df = df.loc[:, ~df.columns.isna()]
    df[NUMERIC_COLS] = df[NUMERIC_COLS].astype(float)
    df["Exonic.gene.sizes"] = pd.to_numeric(df["Exonic.gene.sizes"])
    return df


def add_means_and_log2fc(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Control_mean"] = df[CONTROL_COLS].mean(axis=1)
    df["APOEhi_mean"] = df[APOEHI_COLS].mean(axis=1)
    df["APOElo_mean"] = df[APOELO_COLS].mean(axis=1)
    df["log2FC_APOEhi_vs_Control"] = np.log2(
        (df["APOEhi_mean"] + PSEUDOCOUNT) / (df["Control_mean"] + PSEUDOCOUNT)
    )
    return df


def differential_expression(df: pd.DataFrame) -> pd.DataFrame:
    expression_cols = CONTROL_COLS + APOEHI_COLS
    df_volcano = df[df[expression_cols].mean(axis=1) > 1].copy()
    df_volcano = df_volcano[
        df_volcano["GeneSymbol"].notna() & (df_volcano["GeneSymbol"] != "-")
    ].copy()

    df_volcano["Control_mean"] = df_volcano[CONTROL_COLS].mean(axis=1)
    df_volcano["APOEhi_mean"] = df_volcano[APOEHI_COLS].mean(axis=1)
    df_volcano["log2FC_APOEhi_vs_Control"] = np.log2(
        (df_volcano["APOEhi_mean"] + PSEUDOCOUNT)
        / (df_volcano["Control_mean"] + PSEUDOCOUNT)
    )

    control = np.log2(df_volcano[CONTROL_COLS] + 1)
    apoehi = np.log2(df_volcano[APOEHI_COLS] + 1)

    pvals = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for i in range(len(df_volcano)):
            _, p = stats.ttest_ind(
                apoehi.iloc[i], control.iloc[i], equal_var=False
            )
            pvals.append(p)

    df_volcano["p_value"] = np.nan_to_num(pvals, nan=1.0)

    p = df_volcano["p_value"].values
    m = len(p)
    order = np.argsort(p)
    rank = np.empty(m)
    rank[order] = np.arange(1, m + 1)
    fdr = p * m / rank
    fdr_sorted = fdr[order]
    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    df_volcano["FDR"] = np.empty(m)
    df_volcano.loc[df_volcano.index[order], "FDR"] = np.minimum(fdr_sorted, 1)
    df_volcano["neg_log10_FDR"] = -np.log10(df_volcano["FDR"] + 1e-300)

    df_volcano["significance"] = "Not significant"
    df_volcano.loc[
        (df_volcano["log2FC_APOEhi_vs_Control"] > 1) & (df_volcano["FDR"] < 0.05),
        "significance",
    ] = "Upregulated"
    df_volcano.loc[
        (df_volcano["log2FC_APOEhi_vs_Control"] < -1) & (df_volcano["FDR"] < 0.05),
        "significance",
    ] = "Downregulated"
    return df_volcano


def plot_volcano(df_volcano: pd.DataFrame) -> None:
    colors = {
        "Not significant": "lightgray",
        "Upregulated": "red",
        "Downregulated": "blue",
    }
    plt.figure(figsize=(8, 6))
    for group, color in colors.items():
        subset = df_volcano[df_volcano["significance"] == group]
        plt.scatter(
            subset["log2FC_APOEhi_vs_Control"],
            subset["neg_log10_FDR"],
            c=color,
            label=group,
            alpha=0.6,
            s=12,
        )
    plt.axvline(1, color="black", linestyle="--", linewidth=1)
    plt.axvline(-1, color="black", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)

    top_up = df_volcano.sort_values(
        "log2FC_APOEhi_vs_Control", ascending=False
    ).head(1)
    top_down = df_volcano.sort_values(
        "log2FC_APOEhi_vs_Control", ascending=True
    ).head(1)
    for _, row in pd.concat([top_up, top_down]).iterrows():
        x = row["log2FC_APOEhi_vs_Control"]
        y = row["neg_log10_FDR"]
        offset = (40, 15) if x > 0 else (-55, 18)
        ha = "left" if x > 0 else "right"
        plt.annotate(
            row["GeneSymbol"],
            xy=(x, y),
            xytext=offset,
            textcoords="offset points",
            fontsize=10,
            fontweight="bold",
            color="black",
            ha=ha,
            va="center",
            arrowprops=dict(arrowstyle="-", color="black", linewidth=0.8),
        )

    plt.xlabel("log2 Fold Change (APOEhi vs Control)")
    plt.ylabel("-log10(FDR)")
    plt.title("Volcano Plot - APOEhi vs Control")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "volcano_corrected_labeled_top1.png", dpi=300)
    plt.close()


def plot_pca(df: pd.DataFrame) -> None:
    pca_cols = CONTROL_COLS + APOELO_COLS + APOEHI_COLS
    df_pca = df[df[pca_cols].mean(axis=1) > 1].copy()
    df_pca["variance"] = df_pca[pca_cols].var(axis=1)
    df_pca = df_pca.sort_values("variance", ascending=False).head(2000)

    X = np.log2(df_pca[pca_cols].T + 1)
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2)
    components = pca.fit_transform(X_scaled)
    pca_df = pd.DataFrame(components, columns=["PC1", "PC2"])
    pca_df["condition"] = [
        "Control",
        "Control",
        "Control",
        "APOElo",
        "APOElo",
        "APOElo",
        "APOEhi",
        "APOEhi",
        "APOEhi",
    ]

    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="condition", s=120)
    plt.title("PCA - Transcriptomic Profiles")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    plt.tight_layout()
    plt.savefig(FIGURES / "pca.png", dpi=300)
    plt.close()


def plot_top_genes(df_volcano: pd.DataFrame) -> pd.DataFrame:
    sig_genes = df_volcano[
        (df_volcano["log2FC_APOEhi_vs_Control"].abs() > 1)
        & (df_volcano["FDR"] < 0.05)
    ].copy()
    top_up = sig_genes.sort_values(
        "log2FC_APOEhi_vs_Control", ascending=False
    ).head(10)
    top_down = sig_genes.sort_values(
        "log2FC_APOEhi_vs_Control", ascending=True
    ).head(10)
    top_genes = pd.concat([top_up, top_down])

    plt.figure(figsize=(8, 6))
    plt.barh(top_genes["GeneSymbol"], top_genes["log2FC_APOEhi_vs_Control"])
    plt.axvline(0, color="black", linewidth=1)
    plt.xlabel("log2 Fold Change (APOEhi vs Control)")
    plt.title("Top Differentially Expressed Genes")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(FIGURES / "top_genes_clean.png", dpi=300)
    plt.close()

    heatmap_data = top_genes[
        ["GeneSymbol", "log2FC_APOEhi_vs_Control"]
    ].set_index("GeneSymbol")
    plt.figure(figsize=(6, 8))
    sns.heatmap(
        heatmap_data, cmap="coolwarm", annot=True, fmt=".2f", linewidths=0.5
    )
    plt.title("Heatmap of Top Differentially Expressed Genes")
    plt.tight_layout()
    plt.savefig(FIGURES / "heatmap_clean.png", dpi=300)
    plt.close()
    return top_genes


def export_results(df_volcano: pd.DataFrame) -> Path:
    out = FIGURES / "apoehi_vs_control_de_genes.csv"
    cols = [
        "GeneSymbol",
        "Control_mean",
        "APOEhi_mean",
        "log2FC_APOEhi_vs_Control",
        "p_value",
        "FDR",
        "significance",
    ]
    export = df_volcano[cols].sort_values("FDR")
    export.to_csv(out, index=False)
    return out


def main() -> None:
    data_path = ROOT / "rnaseq-01.xlsx"
    if not data_path.exists():
        raise FileNotFoundError(
            f"Missing {data_path}. See docs/DOWNLOAD_DATA.md"
        )

    df = load_and_clean(data_path)
    df = add_means_and_log2fc(df)
    df_volcano = differential_expression(df)

    counts = df_volcano["significance"].value_counts()
    print(counts.to_string())

    plot_volcano(df_volcano)
    plot_pca(df)
    plot_top_genes(df_volcano)
    out = export_results(df_volcano)
    print(f"Wrote results table: {out}")
    print(f"Wrote figures to: {FIGURES}")


if __name__ == "__main__":
    main()
