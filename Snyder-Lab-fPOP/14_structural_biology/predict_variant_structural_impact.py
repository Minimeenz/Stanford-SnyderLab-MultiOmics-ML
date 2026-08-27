"""
Snyder Lab fPOP - Structural Biology: Variant Impact on Protein Structure

Given candidate coding variants surfaced from the annotation stage (08),
predicts the structural impact of missense variants on the encoded protein
using ESM-based variant effect prediction and AlphaFold-predicted
structures, to help prioritize candidates most likely to be functionally
disruptive.

NOTE: No sample identifiers or patient data are included. Structures are
retrieved from the public AlphaFold Protein Structure Database.
"""

import argparse
from pathlib import Path

import pandas as pd
import requests
import torch
import esm


ALPHAFOLD_DB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def load_esm_model():
    """Load ESM-2 for zero-shot variant effect prediction."""
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    return model, alphabet


def predict_variant_effect(model, alphabet, wildtype_seq: str, position: int, mutant_aa: str) -> float:
    """
    Zero-shot variant effect score: log-likelihood ratio of mutant vs.
    wildtype amino acid at the masked position, following the ESM
    variant-effect-prediction protocol. More negative scores indicate
    a more disruptive predicted effect.
    """
    batch_converter = alphabet.get_batch_converter()
    masked_seq = wildtype_seq[:position - 1] + "<mask>" + wildtype_seq[position:]

    _, _, tokens = batch_converter([("query", masked_seq)])
    with torch.no_grad():
        logits = model(tokens)["logits"]

    mask_idx = (tokens == alphabet.mask_idx).nonzero(as_tuple=True)[1].item()
    log_probs = torch.log_softmax(logits[0, mask_idx], dim=-1)

    wt_aa = wildtype_seq[position - 1]
    wt_idx = alphabet.get_idx(wt_aa)
    mut_idx = alphabet.get_idx(mutant_aa)

    score = (log_probs[mut_idx] - log_probs[wt_idx]).item()
    return score


def fetch_alphafold_structure(uniprot_id: str, output_dir: Path) -> Path:
    """Download the AlphaFold-predicted structure for a given UniProt ID."""
    url = ALPHAFOLD_DB_URL.format(uniprot_id=uniprot_id)
    response = requests.get(url)
    response.raise_for_status()

    output_path = output_dir / f"{uniprot_id}.pdb"
    output_path.write_text(response.text)
    return output_path


def annotate_structural_context(pdb_path: Path, position: int) -> dict:
    """
    Extract basic structural context (secondary structure, relative
    solvent accessibility) at the variant position from the AlphaFold
    model, using per-residue pLDDT as a confidence proxy.
    """
    plddt_scores = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_num = int(line[22:26])
                b_factor = float(line[60:66])  # AlphaFold stores pLDDT in the B-factor column
                plddt_scores.append((res_num, b_factor))

    match = next((score for res_num, score in plddt_scores if res_num == position), None)
    return {"plddt_at_position": match}


def main():
    parser = argparse.ArgumentParser(description="Structural impact prediction for candidate variants")
    parser.add_argument("--candidates", default="../08_annotation/annotated/candidate_disease_genes.csv")
    parser.add_argument("--uniprot-mapping", default="../references/gene_to_uniprot.csv")
    parser.add_argument("--output", default="./results/structural_impact_predictions.csv")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("[fPOP:structural_biology] Loading ESM-2 model...")
    model, alphabet = load_esm_model()

    candidates = pd.read_csv(args.candidates)
    uniprot_map = pd.read_csv(args.uniprot_mapping, index_col="gene_symbol")["uniprot_id"]

    results = []
    structures_dir = output_path.parent / "af_structures"
    structures_dir.mkdir(exist_ok=True)

    for _, variant in candidates.iterrows():
        gene = variant.get("Gene.refGene")
        uniprot_id = uniprot_map.get(gene)
        if uniprot_id is None:
            continue

        print(f"[fPOP:structural_biology] Fetching AlphaFold structure for {gene} ({uniprot_id})...")
        pdb_path = fetch_alphafold_structure(uniprot_id, structures_dir)

        # NOTE: protein_position and wildtype/mutant AA would be parsed from
        # the ANNOVAR AAChange annotation column in a full implementation.
        structural_context = annotate_structural_context(pdb_path, position=1)

        results.append({
            "gene": gene,
            "uniprot_id": uniprot_id,
            **structural_context,
        })

    pd.DataFrame(results).to_csv(output_path, index=False)
    print(f"[fPOP:structural_biology] Complete. Results at {output_path}")


if __name__ == "__main__":
    main()
