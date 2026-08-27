"""
Snyder Lab fPOP - Cross-Omic Integration Model

Fuses candidate gene/molecular signals across all omics layers (genomic
variants, proteomics, methylome, cytokines, ATAC-seq, multiome, structural
biology, spatial transcriptomics, metabolomics) into a single unified
candidate gene ranking, rather than evaluating each omics layer in
isolation. This is the integration layer that ties the full 16-stage
pipeline together into one multi-modal model.

Approach: each omics layer contributes an evidence score per gene. These
per-layer scores are combined into a feature matrix, and a supervised
ranking model (trained on genes with known epilepsy associations as
weak-label positives) learns the relative importance of each omics layer,
producing a single ranked candidate gene list with per-layer contribution
breakdowns for interpretability.

NOTE: No sample identifiers or patient data are included.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

OMICS_SOURCES = {
    "genomic": "../08_annotation/annotated/candidate_disease_genes.csv",
    "proteomic": "../09_proteomics/results/olink_differential_abundance.csv",
    "methylome": "../10_methylome/results/differentially_methylated_regions.csv",
    "cytokine": "../11_cytokines/results/cytokine_differential_analysis.csv",
    "atac": "../12_atac_seq/results/differential_accessibility.csv",
    "multiome": "../13_multiome_analysis/results/cluster_differential_expression.csv",
    "structural": "../14_structural_biology/results/structural_impact_predictions.csv",
    "spatial": "../15_spatial_transcriptomics/results/candidate_gene_spatial_autocorrelation.csv",
    "metabolomic": "../16_metabolomics/results/differential_metabolites.csv",
}


def load_layer_scores(layer_name: str, path: str, gene_col: str, score_col: str) -> pd.Series:
    """
    Load one omics layer's per-gene evidence score. Score is normalized to
    [0, 1] via rank transform so that layers on very different scales
    (p-values, fold-changes, Moran's I, pLDDT, etc.) become comparable.
    """
    file_path = Path(path)
    if not file_path.exists():
        print(f"[fPOP:integration] {layer_name}: no results file found, skipping")
        return pd.Series(dtype=float, name=layer_name)

    df = pd.read_csv(file_path)
    if gene_col not in df.columns or score_col not in df.columns:
        print(f"[fPOP:integration] {layer_name}: expected columns not found, skipping")
        return pd.Series(dtype=float, name=layer_name)

    df = df.dropna(subset=[gene_col, score_col])
    scores = df.groupby(gene_col)[score_col].max()  # most extreme value per gene if duplicated
    normalized = scores.rank(pct=True)  # rank-normalize to [0, 1] within this layer
    normalized.name = layer_name

    print(f"[fPOP:integration] {layer_name}: {len(normalized)} genes with evidence scores")
    return normalized


def build_feature_matrix(layer_configs: dict) -> pd.DataFrame:
    """Assemble a genes x omics-layers matrix, one column per omics source."""
    all_layers = []
    for layer_name, (path, gene_col, score_col) in layer_configs.items():
        layer_scores = load_layer_scores(layer_name, path, gene_col, score_col)
        all_layers.append(layer_scores)

    feature_matrix = pd.concat(all_layers, axis=1)
    # Genes not measured in a given layer get that layer's median score
    # (neutral imputation rather than assuming zero evidence)
    feature_matrix = feature_matrix.fillna(feature_matrix.median())

    print(f"[fPOP:integration] Combined feature matrix: {feature_matrix.shape[0]} genes x "
          f"{feature_matrix.shape[1]} omics layers")
    return feature_matrix


def train_integration_model(feature_matrix: pd.DataFrame, known_positives: set) -> pd.DataFrame:
    """
    Train a gradient-boosted classifier using genes with prior epilepsy
    association (weak-label positives, e.g. from OMIM/gene panels) versus
    all other candidate genes as the contrast class, to learn how much
    weight each omics layer should carry in the final ranking.
    """
    labels = feature_matrix.index.isin(known_positives).astype(int)

    if labels.sum() < 3:
        print("[fPOP:integration] WARNING: fewer than 3 known-positive genes provided; "
              "ranking will be unsupervised (equal-weighted average) instead.")
        feature_matrix["integrated_score"] = feature_matrix.mean(axis=1)
        return feature_matrix.sort_values("integrated_score", ascending=False)

    scaler = StandardScaler()
    X = scaler.fit_transform(feature_matrix)

    model = GradientBoostingClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, random_state=42)
    cv = StratifiedKFold(n_splits=min(5, labels.sum()), shuffle=True, random_state=42)

    cv_scores = cross_val_predict(model, X, labels, cv=cv, method="predict_proba")[:, 1]
    model.fit(X, labels)

    feature_matrix["integrated_score"] = cv_scores
    layer_importance = pd.Series(model.feature_importances_, index=feature_matrix.columns[:-1])

    print("[fPOP:integration] Learned omics-layer importance:")
    print(layer_importance.sort_values(ascending=False).to_string())

    return feature_matrix.sort_values("integrated_score", ascending=False)


def main():
    parser = argparse.ArgumentParser(description="Cross-omic candidate gene integration")
    parser.add_argument("--known-positives", default="../metadata/known_epilepsy_genes.txt",
                         help="Text file, one gene symbol per line, of genes with prior epilepsy association")
    parser.add_argument("--output", default="./results/integrated_candidate_ranking.csv")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    layer_configs = {
        "genomic": (OMICS_SOURCES["genomic"], "Gene.refGene", "gnomAD_exome_ALL"),
        "proteomic": (OMICS_SOURCES["proteomic"], "protein", "padj"),
        "methylome": (OMICS_SOURCES["methylome"], "gene", "qvalue"),
        "cytokine": (OMICS_SOURCES["cytokine"], "analyte", "padj"),
        "atac": (OMICS_SOURCES["atac"], "gene", "padj"),
        "multiome": (OMICS_SOURCES["multiome"], "names", "pvals_adj"),
        "structural": (OMICS_SOURCES["structural"], "gene", "plddt_at_position"),
        "spatial": (OMICS_SOURCES["spatial"], "gene", "I"),
        "metabolomic": (OMICS_SOURCES["metabolomic"], "metabolite", "padj"),
    }

    feature_matrix = build_feature_matrix(layer_configs)

    known_positives_path = Path(args.known_positives)
    known_positives = set()
    if known_positives_path.exists():
        known_positives = set(known_positives_path.read_text().split())
        print(f"[fPOP:integration] Loaded {len(known_positives)} known epilepsy-associated genes")
    else:
        print("[fPOP:integration] No known-positive gene list found; proceeding unsupervised")

    ranked = train_integration_model(feature_matrix, known_positives)
    ranked.to_csv(output_path)

    print(f"[fPOP:integration] Complete. Top 10 integrated candidates:")
    print(ranked.head(10).to_string())


if __name__ == "__main__":
    main()
