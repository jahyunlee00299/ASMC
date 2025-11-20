#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced pocket prediction using multiple geometric and chemical features
Inspired by P2RANK methodology

This uses:
1. Burial depth (neighbor count)
2. Protrusion index (how much residue sticks out locally)
3. Chemical properties (hydrophobicity, charge)
4. Secondary structure context

Author: ASMC
Date: 2025-11-20
"""

import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    from Bio.PDB import PDBParser, NeighborSearch
except ImportError as e:
    print(f"Error: BioPython required. Install: pip install biopython")
    print(f"Import error: {e}")
    sys.exit(1)


# Hydrophobicity scale (Kyte-Doolittle)
HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

# Charge at pH 7
CHARGE = {
    'ARG': 1, 'LYS': 1, 'HIS': 0.1,
    'ASP': -1, 'GLU': -1,
}


def calculate_protrusion_index(structure, chain_id='A', local_radius=8.0, global_radius=12.0):
    """
    Calculate protrusion index for each residue

    Protrusion = how much a residue sticks out relative to its neighbors
    Higher protrusion = more likely to be in a pocket

    Args:
        structure: BioPython structure
        chain_id: Chain to analyze
        local_radius: Radius for local environment (Angstroms)
        global_radius: Radius for global environment (Angstroms)

    Returns:
        dict: {residue_number: protrusion_index}
    """
    chain = structure[0][chain_id]
    all_atoms = list(structure.get_atoms())
    ns = NeighborSearch(all_atoms)

    protrusion_indices = {}

    for residue in chain:
        if residue.id[0] != ' ':
            continue
        if 'CA' not in residue:
            continue

        res_num = residue.id[1]
        ca_coord = residue['CA'].get_coord()

        # Count neighbors in local vs global radius
        local_neighbors = len(ns.search(ca_coord, local_radius))
        global_neighbors = len(ns.search(ca_coord, global_radius))

        # Protrusion = inverse of local density relative to global
        if global_neighbors > 0:
            protrusion = 1.0 - (local_neighbors / global_neighbors)
        else:
            protrusion = 0.0

        protrusion_indices[res_num] = protrusion

    return protrusion_indices


def calculate_pocket_score(structure, chain_id='A'):
    """
    Calculate comprehensive pocket score for each residue

    Combines multiple features:
    - Burial depth (moderate burial preferred)
    - Protrusion index (pockets are cavities)
    - Hydrophobicity (mixed hydrophobic/hydrophilic)
    - Charge distribution

    Returns:
        dict: {residue_number: (score, residue_name, feature_dict)}
    """
    chain = structure[0][chain_id]
    all_atoms = list(structure.get_atoms())
    ns = NeighborSearch(all_atoms)

    # Calculate features
    protrusion = calculate_protrusion_index(structure, chain_id)

    pocket_scores = {}

    for residue in chain:
        if residue.id[0] != ' ':
            continue
        if 'CA' not in residue:
            continue

        res_num = residue.id[1]
        res_name = residue.get_resname()
        ca_coord = residue['CA'].get_coord()

        # Feature 1: Burial depth (10Å radius)
        neighbors_10 = len(ns.search(ca_coord, 10.0))

        # Normalize to 0-1 scale (based on empirical range 50-250)
        burial_norm = (neighbors_10 - 50) / 200.0
        burial_norm = max(0, min(1, burial_norm))

        # Prefer moderate burial (peak at 0.6-0.7)
        burial_score = 1.0 - abs(burial_norm - 0.65) / 0.65

        # Feature 2: Protrusion index
        prot_index = protrusion.get(res_num, 0.0)

        # Feature 3: Hydrophobicity
        hydro = HYDROPHOBICITY.get(res_name, 0.0)
        # Normalize to 0-1 (range -4.5 to 4.5)
        hydro_norm = (hydro + 4.5) / 9.0
        # Prefer moderate hydrophobicity (mixed environments)
        hydro_score = 1.0 - abs(hydro_norm - 0.5) / 0.5

        # Feature 4: Charge
        charge = abs(CHARGE.get(res_name, 0))
        charge_score = charge  # Charged residues often in active sites

        # Combine features with weights
        # These weights are approximations of P2RANK-like scoring
        score = (
            0.35 * burial_score +
            0.25 * prot_index +
            0.20 * hydro_score +
            0.20 * charge_score
        )

        features = {
            'neighbors': neighbors_10,
            'burial_score': burial_score,
            'protrusion': prot_index,
            'hydro_score': hydro_score,
            'charge_score': charge_score
        }

        pocket_scores[res_num] = (score, res_name, features)

    return pocket_scores


def predict_pockets(pocket_scores, top_n=30):
    """
    Predict pocket residues based on scores

    Args:
        pocket_scores: dict from calculate_pocket_score
        top_n: Number of top-scoring residues to return

    Returns:
        list: [(res_num, res_name, score, features), ...]
    """
    # Sort by score
    sorted_residues = sorted(
        [(num, name, score, feat) for num, (score, name, feat) in pocket_scores.items()],
        key=lambda x: x[2],
        reverse=True
    )

    return sorted_residues[:top_n]


def compare_with_known(predicted, known_sites):
    """Compare predictions with known active sites"""
    pred_nums = [r[0] for r in predicted]
    predicted_set = set(pred_nums)
    known_set = set(known_sites)

    true_positive = predicted_set & known_set
    false_positive = predicted_set - known_set
    false_negative = known_set - predicted_set

    precision = len(true_positive) / len(predicted_set) if predicted_set else 0
    recall = len(true_positive) / len(known_set) if known_set else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {
        'true_positive': sorted(true_positive),
        'false_positive': sorted(false_positive),
        'false_negative': sorted(false_negative),
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'n_predicted': len(predicted_set),
        'n_known': len(known_set),
        'n_correct': len(true_positive)
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python predict_pockets_enhanced.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]

    # Known AtUdh active sites
    known_active_sites = [10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137,
                         139, 161, 163, 165, 187, 189, 211, 233, 235, 257, 259]

    print("="*80)
    print("Enhanced Pocket Prediction for AtUdh (3RFV Chain A)")
    print("P2RANK-inspired multi-feature approach")
    print("="*80)
    print(f"PDB file: {pdb_file}")
    print(f"Known active sites: {len(known_active_sites)} residues")
    print()

    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Calculate pocket scores
    print("Calculating multi-feature pocket scores...")
    print("  - Burial depth analysis")
    print("  - Protrusion index calculation")
    print("  - Chemical property scoring")
    print()

    pocket_scores = calculate_pocket_score(structure, chain_id='A')

    # Predict pockets
    predictions = predict_pockets(pocket_scores, top_n=30)

    # Compare with known sites
    results = compare_with_known(predictions, known_active_sites)

    # Display results
    print("Prediction Results:")
    print("-" * 80)
    print(f"Total predicted:  {results['n_predicted']} residues")
    print(f"Total known:      {results['n_known']} residues")
    print(f"Correctly found:  {results['n_correct']} residues")
    print()
    print(f"Precision: {results['precision']:.2%}")
    print(f"Recall:    {results['recall']:.2%}")
    print(f"F1-score:  {results['f1_score']:.2%}")
    print()

    print("Correctly predicted active sites:")
    print(f"  {results['true_positive']}")
    print()

    print("Missed active sites (false negatives):")
    print(f"  {results['false_negative']}")
    print()

    print(f"False positives: {len(results['false_positive'])} residues")
    print()

    # Detailed predictions
    print("Top 30 Pocket Predictions:")
    print("-" * 80)
    print(f"{'Rank':<6} {'Res':<8} {'Name':<6} {'Score':<8} {'Burial':<8} {'Protru':<8} {'Status'}")
    print("-" * 80)

    known_set = set(known_active_sites)
    for i, (res_num, res_name, score, features) in enumerate(predictions, 1):
        status = "[*] ACTIVE" if res_num in known_set else "    pred"
        print(f"{i:<6} {res_num:<8} {res_name:<6} {score:>6.3f}  "
              f"{features['neighbors']:<8} {features['protrusion']:>6.3f}  {status}")

    print("="*80)
    print("\nConclusion:")
    if results['f1_score'] > 0.5:
        print("[GOOD] Strong agreement with known active sites")
    elif results['f1_score'] > 0.3:
        print("[MODERATE] Partial agreement - captures some active sites")
    else:
        print("[POOR] Limited agreement - true P2RANK would perform better")

    print("\nNote: This implements P2RANK-inspired heuristics.")
    print("      Real P2RANK uses machine learning trained on thousands of structures.")
    print("="*80)

    # Save predictions
    output_file = pdb_file.replace('.pdb', '_predicted_pockets.txt')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("# Enhanced pocket prediction results\n")
        f.write(f"# PDB: {pdb_file}\n")
        f.write(f"# F1-score: {results['f1_score']:.3f}\n")
        f.write(f"# Precision: {results['precision']:.3f}\n")
        f.write(f"# Recall: {results['recall']:.3f}\n")
        f.write("#\n")
        f.write("# Residue_Number\tResidue_Name\tPocket_Score\tIs_Active_Site\n")

        for res_num, res_name, score, features in predictions:
            is_active = "YES" if res_num in known_set else "NO"
            f.write(f"{res_num}\t{res_name}\t{score:.4f}\t{is_active}\n")

    print(f"\nPredictions saved to: {output_file}")


if __name__ == "__main__":
    main()
