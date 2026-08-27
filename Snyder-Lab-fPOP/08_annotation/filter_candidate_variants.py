"""
fPOP Epilepsy Family Study — Candidate Variant Filtering
Stanford School of Medicine, Department of Genetics (Snyder Lab)

Filters ANNOVAR-annotated variant tables to rare, predicted-damaging
candidates for downstream cross-omic prioritization. No sample
identifiers or patient data are included.
"""

import argparse

import pandas as pd


def filter_candidates(input_path: str, output_path: str, max_af: float = 0.001):
    df = pd.read_csv(input_path, low_memory=False)

    # Rare in population databases
    df["gnomAD_exome_ALL"] = pd.to_numeric(df.get("gnomAD_exome_ALL", 0), errors="coerce").fillna(0)
    rare = df["gnomAD_exome_ALL"] < max_af

    # Predicted damaging by at least one in-silico tool (CADD, REVEL as examples)
    damaging_cols = [c for c in df.columns if c in ("CADD_phred", "REVEL_score")]
    damaging = pd.Series(False, index=df.index)
    for col in damaging_cols:
        scores = pd.to_numeric(df[col], errors="coerce")
        threshold = 20 if col == "CADD_phred" else 0.5
        damaging |= scores > threshold

    # Exonic / splicing consequence only
    coding = df.get("Func.refGene", "").isin(["exonic", "splicing", "exonic;splicing"])

    candidates = df[rare & damaging & coding].copy()
    candidates = candidates.sort_values("gnomAD_exome_ALL")

    candidates.to_csv(output_path, index=False)
    print(f"[fPOP:annotation] {len(candidates)} candidate variants written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter fPOP candidate disease variants")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-af", type=float, default=0.001)
    args = parser.parse_args()

    filter_candidates(args.input, args.output, max_af=args.max_af)
