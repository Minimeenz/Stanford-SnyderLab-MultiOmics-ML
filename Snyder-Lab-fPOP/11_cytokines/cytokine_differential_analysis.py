"""
fPOP Epilepsy Family Study - Stage 07: Cytokine / Immune Profiling Analysis
Stanford School of Medicine, Department of Genetics (Snyder Lab)

Analyzes multiplex cytokine panel data (e.g., Luminex-style immunoassay
output) across the proband and family members to assess whether immune/
inflammatory signaling differs in the proband, consistent with growing
evidence of neuroinflammatory involvement in epilepsy.

NOTE: No sample identifiers or patient data are included.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def load_cytokine_panel(panel_path: str) -> pd.DataFrame:
    """Load multiplex cytokine concentration matrix: samples x analytes (pg/mL)."""
    df = pd.read_csv(panel_path, index_col=0)
    print(f"[fPOP:cytokine] Loaded panel: {df.shape[0]} samples x {df.shape[1]} analytes")
    return df


def log_transform_and_impute(df: pd.DataFrame) -> pd.DataFrame:
    """Log2-transform concentrations and impute below-detection-limit values at half the panel minimum."""
    imputed = df.copy()
    for col in imputed.columns:
        min_val = imputed[col][imputed[col] > 0].min()
        imputed[col] = imputed[col].fillna(min_val / 2).clip(lower=min_val / 2)
    return np.log2(imputed)


def differential_cytokines(log_panel: pd.DataFrame, group_labels: pd.Series) -> pd.DataFrame:
    """Compare log2 cytokine concentrations between proband and unaffected family members."""
    proband_samples = group_labels[group_labels == "proband"].index
    family_samples = group_labels[group_labels == "family"].index

    results = []
    for analyte in log_panel.columns:
        proband_vals = log_panel.loc[proband_samples, analyte].dropna()
        family_vals = log_panel.loc[family_samples, analyte].dropna()

        if len(proband_vals) == 0 or len(family_vals) == 0:
            continue

        stat, pval = stats.mannwhitneyu(proband_vals, family_vals, alternative="two-sided")
        results.append({
            "analyte": analyte,
            "proband_mean_log2": proband_vals.mean(),
            "family_mean_log2": family_vals.mean(),
            "log2_diff": proband_vals.mean() - family_vals.mean(),
            "pvalue": pval,
        })

    results_df = pd.DataFrame(results)
    results_df["padj"] = multipletests(results_df["pvalue"], method="fdr_bh")[1]
    return results_df.sort_values("padj")


def main():
    parser = argparse.ArgumentParser(description="fPOP cytokine panel differential analysis")
    parser.add_argument("--panel", default="./cytokine_panel_matrix.csv")
    parser.add_argument("--groups", default="../metadata/proband_vs_family_groups.csv")
    parser.add_argument("--output", default="./results")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    panel = load_cytokine_panel(args.panel)
    log_panel = log_transform_and_impute(panel)

    group_labels = pd.read_csv(args.groups, index_col=0)["group"]

    results = differential_cytokines(log_panel, group_labels)
    results.to_csv(output_dir / "cytokine_differential_analysis.csv", index=False)

    n_significant = (results["padj"] < 0.05).sum()
    print(f"[fPOP:cytokine] {n_significant} analytes differentially abundant (padj < 0.05)")
    print(results.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
