"""
fPOP Epilepsy Family Study - Multiome Analysis (Single-Cell RNA + ATAC)
Stanford School of Medicine, Department of Genetics (Snyder Lab)

Processes single-cell multiome (joint RNA + ATAC) data to characterize
cell-type-specific transcriptional and chromatin accessibility changes in
the proband relative to family members, integrating with bulk-level
findings from earlier pipeline stages.

NOTE: No sample identifiers or patient data are included.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def load_multiome_data(rna_path: str, atac_path: str) -> ad.AnnData:
    """Load paired single-cell RNA and ATAC count matrices into a joint AnnData object."""
    rna = sc.read_10x_h5(rna_path)
    atac = sc.read_10x_h5(atac_path)

    rna.var_names_make_unique()
    print(f"[fPOP:multiome] RNA: {rna.n_obs} cells x {rna.n_vars} genes")
    print(f"[fPOP:multiome] ATAC: {atac.n_obs} cells x {atac.n_vars} peaks")

    rna.obsm["ATAC"] = atac.X
    return rna


def preprocess_rna(adata: ad.AnnData, min_genes: int = 200, max_pct_mt: float = 15.0) -> ad.AnnData:
    """Standard scRNA-seq QC: filter low-quality cells, normalize, log-transform."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[adata.obs.n_genes_by_counts > min_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    print(f"[fPOP:multiome] Post-QC: {adata.n_obs} cells retained")
    return adata


def cluster_and_annotate(adata: ad.AnnData, resolution: float = 0.8) -> ad.AnnData:
    """PCA -> neighbors -> Leiden clustering for cell-type identification."""
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.leiden(adata, resolution=resolution)
    sc.tl.umap(adata)

    print(f"[fPOP:multiome] Identified {adata.obs['leiden'].nunique()} cell clusters")
    return adata


def differential_expression_by_group(adata: ad.AnnData, group_col: str = "proband_status") -> pd.DataFrame:
    """
    Per-cluster differential expression between proband and family cells,
    to identify cell-type-specific transcriptional signatures.
    """
    results = []
    for cluster in adata.obs["leiden"].unique():
        cluster_data = adata[adata.obs["leiden"] == cluster]
        sc.tl.rank_genes_groups(cluster_data, group_col, groups=["proband"], reference="family", method="wilcoxon")

        de_df = sc.get.rank_genes_groups_df(cluster_data, group="proband")
        de_df["cluster"] = cluster
        results.append(de_df)

    return pd.concat(results, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description="fPOP single-cell multiome analysis")
    parser.add_argument("--rna", default="./raw_data/multiome_rna.h5")
    parser.add_argument("--atac", default="./raw_data/multiome_atac.h5")
    parser.add_argument("--output", default="./results")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_multiome_data(args.rna, args.atac)
    adata = preprocess_rna(adata)
    adata = cluster_and_annotate(adata)

    de_results = differential_expression_by_group(adata)
    de_results.to_csv(output_dir / "cluster_differential_expression.csv", index=False)

    adata.write(output_dir / "processed_multiome.h5ad")

    print(f"[fPOP:multiome] Complete. Results written to {output_dir}")


if __name__ == "__main__":
    main()
