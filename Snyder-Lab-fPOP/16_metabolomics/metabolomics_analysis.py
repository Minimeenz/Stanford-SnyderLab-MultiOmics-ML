"""
Snyder Lab fPOP — Metabolomics Analysis

Processes untargeted/targeted metabolomics data (LC-MS/MS) across the
proband and family members, completing the classical four-omics panel
(genomics, transcriptomics, proteomics, metabolomics) and adding direct
clinical relevance given known metabolic contributions to seizure
etiology (e.g., ketone body metabolism, neurotransmitter precursor
pathways).

NOTE: No sample identifiers or patient data are included.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


def load_metabolomics_data(peak_table_path: str) -> pd.DataFrame:
    """Load LC-MS/MS peak table: samples x metabolite features (normalized intensities)."""
    df = pd.read_csv(peak_table_path, index_col=0)
    print(f"[fPOP:metabolomics] Loaded {df.shape[0]} samples x {df.shape[1]} metabolite features")
    return df


def normalize_and_scale(df: pd.DataFrame) -> pd.DataFrame:
    """Log-transform and Pareto-scale, standard practice for untargeted metabolomics."""
    log_transformed = np.log2(df.replace(0, np.nan))
    log_transformed = log_transformed.fillna(log_transformed.min().min())

    scaler = StandardScaler()
    scaled = pd.DataFrame(
        scaler.fit_transform(log_transformed),
        index=log_transformed.index,
        columns=log_transformed.columns,
    )
    return scaled


def pca_overview(scaled_df: pd.DataFrame, n_components: int = 5) -> pd.DataFrame:
    """PCA to visualize overall metabolomic separation between proband and family."""
    pca = PCA(n_components=n_components, random_state=42)
    components = pca.fit_transform(scaled_df)
    print(f"[fPOP:metabolomics] PCA explained variance (top {n_components}): "
          f"{pca.explained_variance_ratio_.sum():.3f}")
    return pd.DataFrame(components, index=scaled_df.index,
                         columns=[f"PC{i+1}" for i in range(n_components)])


def differential_metabolites(scaled_df: pd.DataFrame, group_labels: pd.Series) -> pd.DataFrame:
    """Compare metabolite levels between proband and unaffected family members."""
    proband_samples = group_labels[group_labels == "proband"].index
    family_samples = group_labels[group_labels == "family"].index

    results = []
    for metabolite in scaled_df.columns:
        proband_vals = scaled_df.loc[proband_samples, metabolite].dropna()
        family_vals = scaled_df.loc[family_samples, metabolite].dropna()

        if len(proband_vals) == 0 or len(family_vals) == 0:
            continue

        stat, pval = stats.mannwhitneyu(proband_vals, family_vals, alternative="two-sided")
        results.append({
            "metabolite": metabolite,
            "proband_mean": proband_vals.mean(),
            "family_mean": family_vals.mean(),
            "diff": proband_vals.mean() - family_vals.mean(),
            "pvalue": pval,
        })

    results_df = pd.DataFrame(results)
    results_df["padj"] = multipletests(results_df["pvalue"], method="fdr_bh")[1]
    return results_df.sort_values("padj")


def main():
    parser = argparse.ArgumentParser(description="fPOP metabolomics differential analysis")
    parser.add_argument("--peak-table", default="./metabolomics_peak_table.csv")
    parser.add_argument("--groups", default="../metadata/proband_vs_family_groups.csv")
    parser.add_argument("--output", default="./results")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw = load_metabolomics_data(args.peak_table)
    scaled = normalize_and_scale(raw)

    pca_scores = pca_overview(scaled)
    pca_scores.to_csv(output_dir / "pca_scores.csv")

    group_labels = pd.read_csv(args.groups, index_col=0)["group"]
    results = differential_metabolites(scaled, group_labels)
    results.to_csv(output_dir / "differential_metabolites.csv", index=False)

    n_significant = (results["padj"] < 0.05).sum()
    print(f"[fPOP:metabolomics] {n_significant} metabolites differentially abundant (padj < 0.05)")
    print(results.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
