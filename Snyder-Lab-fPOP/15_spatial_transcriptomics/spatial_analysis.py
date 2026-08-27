"""
Snyder Lab fPOP — Spatial Transcriptomics Analysis

Processes spatial transcriptomics data (e.g., 10x Visium / Xenium) to map
candidate gene expression (from Stage 08 annotation and Stage 09/10
proteomic/methylation findings) onto tissue spatial context, identifying
whether candidate genes show spatially localized expression patterns
relevant to disease-affected regions.

NOTE: No sample identifiers or patient data are included.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


def load_spatial_data(visium_path: str) -> "sc.AnnData":
    """Load 10x Visium spatial transcriptomics data with spatial coordinates."""
    adata = sc.read_visium(visium_path)
    adata.var_names_make_unique()
    print(f"[fPOP:spatial] Loaded {adata.n_obs} spots x {adata.n_vars} genes")
    return adata


def preprocess_spatial(adata: "sc.AnnData", min_counts: int = 500) -> "sc.AnnData":
    """Standard spatial transcriptomics QC and normalization."""
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    return adata


def spatial_cluster(adata: "sc.AnnData", resolution: float = 0.5) -> "sc.AnnData":
    """Cluster spots by expression and compute a spatial neighbors graph."""
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=resolution)

    sq.gr.spatial_neighbors(adata)
    return adata


def candidate_gene_spatial_profile(adata: "sc.AnnData", candidate_genes: list[str]) -> pd.DataFrame:
    """
    For each candidate gene, compute spatial autocorrelation (Moran's I) to
    assess whether its expression is spatially structured rather than
    randomly distributed across the tissue — a signal of regional disease
    relevance.
    """
    genes_present = [g for g in candidate_genes if g in adata.var_names]
    print(f"[fPOP:spatial] {len(genes_present)}/{len(candidate_genes)} candidate genes found in spatial data")

    sq.gr.spatial_autocorr(adata, mode="moran", genes=genes_present)
    moran_results = adata.uns["moranI"].loc[genes_present]
    return moran_results.sort_values("I", ascending=False)


def main():
    parser = argparse.ArgumentParser(description="fPOP spatial transcriptomics analysis")
    parser.add_argument("--visium-path", default="./raw_data/visium_sample")
    parser.add_argument("--candidates", default="../08_annotation/annotated/candidate_disease_genes.csv")
    parser.add_argument("--output", default="./results")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_spatial_data(args.visium_path)
    adata = preprocess_spatial(adata)
    adata = spatial_cluster(adata)

    candidates_df = pd.read_csv(args.candidates)
    candidate_genes = candidates_df["Gene.refGene"].dropna().unique().tolist()

    moran_results = candidate_gene_spatial_profile(adata, candidate_genes)
    moran_results.to_csv(output_dir / "candidate_gene_spatial_autocorrelation.csv")

    adata.write(output_dir / "processed_spatial.h5ad")

    print(f"[fPOP:spatial] Complete. Top spatially structured candidate genes:")
    print(moran_results.head(10).to_string())


if __name__ == "__main__":
    main()
