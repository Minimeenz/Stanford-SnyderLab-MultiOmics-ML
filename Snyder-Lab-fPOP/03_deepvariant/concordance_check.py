"""
Snyder Lab fPOP — DeepVariant / GATK Caller Concordance Check

Compares variant calls between the deep-learning-based DeepVariant caller
and the GATK HaplotypeCaller pipeline to identify high-confidence
concordant calls, which are prioritized as candidates for downstream
annotation and family segregation analysis.

NOTE: No sample identifiers or patient data are included.
"""

import argparse
from pathlib import Path

import pandas as pd
from cyvcf2 import VCF


def load_variant_positions(vcf_path: str) -> set:
    """Extract (chrom, pos, ref, alt) tuples from a VCF for set-based comparison."""
    positions = set()
    for variant in VCF(vcf_path):
        for alt in variant.ALT:
            positions.add((variant.CHROM, variant.POS, variant.REF, alt))
    return positions


def compute_concordance(deepvariant_dir: str, gatk_dir: str) -> pd.DataFrame:
    dv_path = Path(deepvariant_dir)
    gatk_path = Path(gatk_dir)

    results = []
    for dv_vcf in dv_path.glob("*.deepvariant.vcf.gz"):
        sample_id = dv_vcf.name.replace(".deepvariant.vcf.gz", "")
        gatk_vcf = gatk_path / f"{sample_id}.g.vcf.gz"

        if not gatk_vcf.exists():
            print(f"[fPOP:concordance] No matching GATK call set for {sample_id}, skipping")
            continue

        dv_variants = load_variant_positions(str(dv_vcf))
        gatk_variants = load_variant_positions(str(gatk_vcf))

        concordant = dv_variants & gatk_variants
        dv_only = dv_variants - gatk_variants
        gatk_only = gatk_variants - dv_variants

        results.append({
            "sample_id": sample_id,
            "deepvariant_total": len(dv_variants),
            "gatk_total": len(gatk_variants),
            "concordant": len(concordant),
            "deepvariant_only": len(dv_only),
            "gatk_only": len(gatk_only),
            "concordance_rate": len(concordant) / max(len(dv_variants | gatk_variants), 1),
        })

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="DeepVariant / GATK concordance check")
    parser.add_argument("--deepvariant-dir", required=True)
    parser.add_argument("--gatk-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    summary = compute_concordance(args.deepvariant_dir, args.gatk_dir)
    summary.to_csv(args.output, index=False)

    print(f"[fPOP:concordance] Mean concordance rate: {summary['concordance_rate'].mean():.3f}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
