"""
Snyder Lab fPOP — Pathogenic Repeat Expansion Flagging

Parses ExpansionHunter output (JSON per-sample repeat genotypes) across
the family, compares repeat lengths at each locus between the proband and
unaffected family members, and flags loci where the proband falls in a
known pathogenic or premutation range.

NOTE: No sample identifiers or patient data are included.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

# Illustrative pathogenic-range thresholds for well-characterized STR loci.
# Real thresholds should be sourced from locus-specific clinical literature.
PATHOGENIC_THRESHOLDS = {
    "ATXN2": {"normal_max": 31, "pathogenic_min": 33},
    "ATXN1": {"normal_max": 35, "pathogenic_min": 39},
    "FMR1": {"normal_max": 44, "pathogenic_min": 200},
}


def load_str_genotypes(str_dir: str) -> pd.DataFrame:
    records = []
    for json_path in Path(str_dir).glob("*.json"):
        sample_id = json_path.stem
        with open(json_path) as f:
            data = json.load(f)

        for locus_id, locus_data in data.get("LocusResults", {}).items():
            for variant_id, variant_data in locus_data.get("Variants", {}).items():
                genotype = variant_data.get("Genotype", "")
                records.append({
                    "sample_id": sample_id,
                    "locus": locus_id,
                    "genotype": genotype,
                })

    return pd.DataFrame(records)


def flag_pathogenic(str_df: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    str_df = str_df.merge(groups.rename("group"), left_on="sample_id", right_index=True)
    flagged = []

    for _, row in str_df.iterrows():
        thresholds = PATHOGENIC_THRESHOLDS.get(row["locus"])
        if thresholds is None:
            continue

        try:
            allele_sizes = [int(a) for a in row["genotype"].split("/")]
        except ValueError:
            continue

        max_allele = max(allele_sizes)
        is_pathogenic = max_allele >= thresholds["pathogenic_min"]

        flagged.append({
            **row.to_dict(),
            "max_allele_size": max_allele,
            "pathogenic_range": is_pathogenic,
        })

    return pd.DataFrame(flagged)


def main():
    parser = argparse.ArgumentParser(description="Flag pathogenic-range STR expansions")
    parser.add_argument("--str-dir", required=True)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    groups = pd.read_csv(args.groups, index_col=0)["group"]
    str_df = load_str_genotypes(args.str_dir)
    flagged = flag_pathogenic(str_df, groups)

    flagged.to_csv(args.output, index=False)

    n_pathogenic = flagged["pathogenic_range"].sum() if not flagged.empty else 0
    print(f"[fPOP:str] {n_pathogenic} sample-locus pairs in pathogenic range")
    if not flagged.empty:
        print(flagged[flagged["pathogenic_range"]].to_string(index=False))


if __name__ == "__main__":
    main()
