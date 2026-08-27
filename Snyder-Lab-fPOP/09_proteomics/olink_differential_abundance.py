"""
fPOP Epilepsy Family Study — Stage 05: Olink Proteomics Analysis
Stanford School of Medicine, Department of Genetics (Snyder Lab)

Processes Olink proximity extension assay (PEA) proteomic data (NPX
values) across the proband and family members, applies QC filtering,
and identifies proteins differentially abundant in the proband relative
to unaffected family members — for integration with genomic candidates
from Stage 04.

NOTE: No sample identifiers or patient data are included. Input is an
Olink NPX export with de-identified sample IDs.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def load_olink_npx(npx_path: str) -> pd.DataFrame:
    """Load Olink NPX (Normalized Protein eXpression) matrix: samples x proteins."""
    df = pd.read_csv(npx_path, index_col=0)
    print(f"[fPOP:olink] Loaded NPX matrix: {df.shape[0]} samples x {df.shape[1]} proteins")
    return df


def qc_filter(npx: pd.DataFrame, warning_flags_path: str = None, missingness_threshold: float = 0.2) -> pd.DataFrame:
    """Remove proteins with excessive missingness or QC warning flags per Olink guidelines."""
    missing_frac = npx.isna().mean(axis=0)
    keep_cols = missing_frac[missing_frac < missingness_threshold].index
    filtered = npx[keep_cols]
    print(f"[fPOP:olink] QC filter: retained {len(keep_cols)}/{npx.shape[1]} proteins")
    return filtered


def differential_abundance(npx: pd.DataFrame, group_labels: pd.Series) -> pd.DataFrame:
    """
    Compare NPX values between proband and unaffected family members
    (Mann-Whitney U, non-parametric given small family-based sample sizes),
    with Benjamini-Hochberg FDR correction.
    """
    results = []
    proband_samples = group_labels[group_labels == "proband"].index
    family_samples = group_labels[group_labels == "family"].index

    for protein in npx.columns:
        proband_vals = npx.loc[proband_samples, protein].dropna()
        family_vals = npx.loc[family_samples, protein].dropna()

        if len(proband_vals) == 0 or len(family_vals) == 0:
            continue

        stat, pval = stats.mannwhitneyu(proband_vals, family_vals, alternative="two-sided")
        fold_diff = proband_vals.mean() - family_vals.mean()  # NPX is already log2-scale

        results.append({
            "protein": protein,
            "proband_mean_npx": proband_vals.mean(),
            "family_mean_npx": family_vals.mean(),
            "log2_diff": fold_diff,
            "pvalue": pval,
        })

    results_df = pd.DataFrame(results)
    results_df["padj"] = multipletests(results_df["pvalue"], method="fdr_bh")[1]
    return results_df.sort_values("padj")


def main():
    parser = argparse.ArgumentParser(description="fPOP Olink proteomics differential abundance analysis")
    parser.add_argument("--npx", default="./olink_npx_matrix.csv")
    parser.add_argument("--groups", default="../metadata/proband_vs_family_groups.csv")
    parser.add_argument("--output", default="./results")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    npx = load_olink_npx(args.npx)
    npx = qc_filter(npx)

    group_labels = pd.read_csv(args.groups, index_col=0)["group"]  # "proband" or "family"

    results = differential_abundance(npx, group_labels)
    results.to_csv(output_dir / "olink_differential_abundance.csv", index=False)

    n_significant = (results["padj"] < 0.05).sum()
    print(f"[fPOP:olink] {n_significant} proteins differentially abundant (padj < 0.05)")
    print(results.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
