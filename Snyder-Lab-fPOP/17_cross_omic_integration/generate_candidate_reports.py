"""
Snyder Lab fPOP - Cross-Omic Integration: Candidate Report Generation

Generates a human-readable per-candidate-gene report explaining *why* each
gene was ranked highly by the cross-omic integration model, showing its
evidence score in each omics layer, so results are interpretable rather
than a single opaque ranking number. This directly supports the kind of
evaluation-beyond-accuracy transparency emphasized in biological ML model
assessment.
"""

import argparse
from pathlib import Path

import pandas as pd


def generate_gene_report(ranked_df: pd.DataFrame, gene: str, top_n_layers: int = 3) -> str:
    """Build a plain-language summary of which omics layers support a given candidate gene."""
    if gene not in ranked_df.index:
        return f"{gene}: not found in integrated candidate list."

    row = ranked_df.loc[gene]
    layer_scores = row.drop("integrated_score").sort_values(ascending=False)
    top_layers = layer_scores.head(top_n_layers)

    lines = [f"Gene: {gene}"]
    lines.append(f"  Integrated score: {row['integrated_score']:.3f}")
    lines.append(f"  Strongest supporting omics layers:")
    for layer_name, score in top_layers.items():
        lines.append(f"    - {layer_name}: percentile rank {score:.2f}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate interpretable per-gene cross-omic reports")
    parser.add_argument("--ranked-candidates", default="./results/integrated_candidate_ranking.csv")
    parser.add_argument("--top-n-genes", type=int, default=15)
    parser.add_argument("--output", default="./results/candidate_gene_reports.txt")
    args = parser.parse_args()

    ranked_df = pd.read_csv(args.ranked_candidates, index_col=0)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    top_genes = ranked_df.head(args.top_n_genes).index

    reports = []
    for gene in top_genes:
        reports.append(generate_gene_report(ranked_df, gene))

    full_report = "\n\n".join(reports)
    output_path.write_text(full_report)

    print(f"[fPOP:integration] Generated reports for top {len(top_genes)} candidate genes")
    print(full_report)


if __name__ == "__main__":
    main()
