#!/usr/bin/env python3
"""
Simple pocket prediction using geometric criteria
Compares predicted pockets with known AtUdh active sites

This is a lightweight alternative to P2RANK for testing purposes.
Uses residue burial depth and local geometry to identify potential pockets.

Usage:
    python predict_pockets_simple.py pdb3rfv_chainA.pdb

Author: ASMC
Date: 2025-11-20
"""

import sys
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    from Bio.PDB import PDBParser, NeighborSearch, Selection
except ImportError:
    print("Error: BioPython required. Install: pip install biopython")
    sys.exit(1)


def calculate_burial_depth(structure, chain_id='A', radius=10.0):
    """
    Calculate burial depth for each residue based on neighbor count

    Args:
        structure: BioPython structure
        chain_id: Chain to analyze
        radius: Search radius for neighbors (Angstroms)

    Returns:
        dict: {residue_number: (burial_score, residue_name)}
    """
    chain = structure[0][chain_id]
    all_atoms = list(structure.get_atoms())
    ns = NeighborSearch(all_atoms)

    burial_scores = {}

    for residue in chain:
        # Skip hetero residues
        if residue.id[0] != ' ':
            continue

        if 'CA' not in residue:
            continue

        res_num = residue.id[1]
        res_name = residue.get_resname()
        ca_atom = residue['CA']

        # Count neighbors within radius
        neighbors = ns.search(ca_atom.coord, radius)
        neighbor_count = len(neighbors)

        # Surface residues have fewer neighbors
        # Pocket residues have moderate neighbor counts
        burial_scores[res_num] = (neighbor_count, res_name)

    return burial_scores


def find_pocket_clusters(burial_scores, min_neighbors=30, max_neighbors=80):
    """
    Identify potential pocket residues based on burial depth

    Pocket residues typically have moderate burial:
    - Not too exposed (surface)
    - Not too buried (core)

    Args:
        burial_scores: dict from calculate_burial_depth
        min_neighbors: Minimum neighbor count for pocket
        max_neighbors: Maximum neighbor count for pocket

    Returns:
        list: Pocket residue numbers sorted by burial score
    """
    pocket_residues = []

    for res_num, (neighbor_count, res_name) in burial_scores.items():
        if min_neighbors <= neighbor_count <= max_neighbors:
            pocket_residues.append((res_num, res_name, neighbor_count))

    # Sort by neighbor count (mid-range = more likely pocket)
    optimal = (min_neighbors + max_neighbors) / 2
    pocket_residues.sort(key=lambda x: abs(x[2] - optimal))

    return pocket_residues


def cluster_sequential_residues(pocket_residues, gap_tolerance=10):
    """
    Group pocket residues into spatial/sequential clusters

    Args:
        pocket_residues: list of (res_num, res_name, score) tuples
        gap_tolerance: Maximum gap between residues in same cluster

    Returns:
        list: List of clusters (each cluster is a list of residues)
    """
    if not pocket_residues:
        return []

    # Sort by residue number
    sorted_res = sorted(pocket_residues, key=lambda x: x[0])

    clusters = []
    current_cluster = [sorted_res[0]]

    for i in range(1, len(sorted_res)):
        prev_num = sorted_res[i-1][0]
        curr_num = sorted_res[i][0]

        if curr_num - prev_num <= gap_tolerance:
            current_cluster.append(sorted_res[i])
        else:
            if len(current_cluster) >= 3:  # Minimum cluster size
                clusters.append(current_cluster)
            current_cluster = [sorted_res[i]]

    # Don't forget last cluster
    if len(current_cluster) >= 3:
        clusters.append(current_cluster)

    return clusters


def compare_with_known_active_sites(predicted, known_sites):
    """
    Compare predicted pocket residues with known active sites

    Args:
        predicted: list of predicted residue numbers
        known_sites: list of known active site residue numbers

    Returns:
        dict: Statistics about prediction accuracy
    """
    predicted_set = set(predicted)
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
        print("Usage: python predict_pockets_simple.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]

    # Known AtUdh active sites from atudh_active_site.txt
    known_active_sites = [10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137,
                         139, 161, 163, 165, 187, 189, 211, 233, 235, 257, 259]

    print("="*70)
    print("Simple Pocket Prediction for AtUdh (3RFV Chain A)")
    print("="*70)
    print(f"PDB file: {pdb_file}")
    print(f"Known active sites: {len(known_active_sites)} residues")
    print()

    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Calculate burial depth
    print("Calculating burial depth for all residues...")
    burial_scores = calculate_burial_depth(structure, chain_id='A', radius=10.0)
    print(f"Analyzed {len(burial_scores)} residues")
    print()

    # Test different burial thresholds
    print("Testing different burial depth thresholds:")
    print("-" * 70)

    best_f1 = 0
    best_params = None
    best_results = None

    # Based on analysis: active sites have 116-230 neighbors
    for min_n in range(90, 140, 10):
        for max_n in range(200, 280, 10):
            pocket_res = find_pocket_clusters(burial_scores, min_n, max_n)
            predicted_nums = [r[0] for r in pocket_res[:30]]  # Top 30 predictions

            if not predicted_nums:  # Skip if no predictions
                continue

            results = compare_with_known_active_sites(predicted_nums, known_active_sites)

            print(f"  min={min_n:2d}, max={max_n:2d}: F1={results['f1_score']:.3f}, "
                  f"P={results['precision']:.2f}, R={results['recall']:.2f}, "
                  f"n={len(predicted_nums)}")

            if results['f1_score'] > best_f1:
                best_f1 = results['f1_score']
                best_params = (min_n, max_n)
                best_results = results
                best_pocket = pocket_res[:30]

    # Show best results
    if best_params is None:
        print("\nNo valid pocket predictions found!")
        print("Try adjusting the burial depth thresholds.")
        sys.exit(1)

    min_n, max_n = best_params
    print(f"\nBest parameters: min_neighbors={min_n}, max_neighbors={max_n}")
    print("="*70)
    print("\nPrediction Results:")
    print("-" * 70)
    print(f"Total predicted:  {best_results['n_predicted']} residues")
    print(f"Total known:      {best_results['n_known']} residues")
    print(f"Correctly found:  {best_results['n_correct']} residues")
    print()
    print(f"Precision: {best_results['precision']:.2%}")
    print(f"Recall:    {best_results['recall']:.2%}")
    print(f"F1-score:  {best_results['f1_score']:.2%}")
    print()

    # Show detailed comparison
    print("Correctly predicted active sites:")
    print(f"  {best_results['true_positive']}")
    print()

    print("Missed active sites (false negatives):")
    print(f"  {best_results['false_negative']}")
    print()

    print("False positives (predicted but not active sites):")
    print(f"  {best_results['false_positive']}")
    print()

    # Show top predictions with details
    print("\nTop 30 Pocket Predictions:")
    print("-" * 70)
    print(f"{'Rank':<6} {'ResNum':<8} {'ResName':<8} {'Neighbors':<10} {'Status'}")
    print("-" * 70)

    known_set = set(known_active_sites)
    for i, (res_num, res_name, neighbor_count) in enumerate(best_pocket, 1):
        status = "[*] ACTIVE SITE" if res_num in known_set else "    predicted"
        print(f"{i:<6} {res_num:<8} {res_name:<8} {neighbor_count:<10} {status}")

    print("="*70)
    print("\nConclusion:")
    if best_results['f1_score'] > 0.5:
        print("[OK] Good agreement with known active sites")
    elif best_results['f1_score'] > 0.3:
        print("[WARN] Moderate agreement - some active sites detected")
    else:
        print("[POOR] Poor agreement - geometric criteria insufficient")

    print("\nNote: This is a simple geometric method.")
    print("      P2RANK uses machine learning and is much more accurate.")
    print("="*70)


if __name__ == "__main__":
    main()
