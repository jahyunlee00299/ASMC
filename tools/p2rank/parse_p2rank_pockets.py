#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse P2RANK predictions.csv and compare with known active sites

P2RANK format:
name,rank,score,probability,sas_points,surf_atoms,center_x,center_y,center_z,residue_ids,surf_atom_ids
pocket1,1,14.75,0.733,102,53,75.96,49.16,16.64,A_109 A_110 A_111 ...

Author: ASMC
Date: 2025-11-20
"""

import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import csv
from pathlib import Path


# Known AtUdh active sites
KNOWN_ACTIVE_SITES = [10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137,
                      139, 161, 163, 165, 187, 189, 211, 233, 235, 257, 259]


def parse_p2rank_csv(csv_file):
    """
    Parse P2RANK predictions.csv file

    Returns:
        list: [(pocket_name, rank, score, residues), ...]
        where residues is list of residue numbers
    """
    pockets = []

    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            # Strip whitespace from keys and values
            row = {k.strip(): v.strip() if isinstance(v, str) else v for k, v in row.items()}

            pocket_name = row.get('name', 'unknown')
            rank = int(row.get('rank', 0))
            score = float(row.get('score', 0.0))
            probability = float(row.get('probability', 0.0))

            # Parse residue_ids: "A_109 A_110 A_111 ..."
            residue_ids_str = row.get('residue_ids', '')
            residues = []

            for res_id in residue_ids_str.split():
                # Format: A_109 → residue number 109
                if '_' in res_id:
                    parts = res_id.split('_')
                    try:
                        res_num = int(parts[1])
                        residues.append(res_num)
                    except (IndexError, ValueError):
                        continue

            pockets.append({
                'name': pocket_name,
                'rank': rank,
                'score': score,
                'probability': probability,
                'residues': sorted(residues),
                'n_residues': len(residues)
            })

    return sorted(pockets, key=lambda x: x['rank'])


def compare_with_known(pocket_residues, known_sites):
    """Compare pocket residues with known active sites"""
    predicted_set = set(pocket_residues)
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
        print("Usage: python parse_p2rank_pockets.py <predictions.csv>")
        sys.exit(1)

    csv_file = sys.argv[1]

    if not Path(csv_file).exists():
        print(f"Error: File not found: {csv_file}")
        sys.exit(1)

    print("="*80)
    print("P2RANK Pocket Analysis for AtUdh (3RFV)")
    print("="*80)
    print(f"File: {csv_file}")
    print(f"Known active sites: {len(KNOWN_ACTIVE_SITES)} residues")
    print()

    # Parse pockets
    pockets = parse_p2rank_csv(csv_file)

    print(f"P2RANK found {len(pockets)} pockets:")
    print()

    # Show all pockets
    print("="*80)
    print("All Predicted Pockets")
    print("="*80)
    print(f"{'Rank':<6} {'Name':<10} {'Score':<8} {'Prob':<8} {'Residues':<10} {'Preview'}")
    print("-"*80)

    for pocket in pockets:
        residues_preview = ' '.join(str(r) for r in pocket['residues'][:5])
        if len(pocket['residues']) > 5:
            residues_preview += f" ... +{len(pocket['residues'])-5}"

        print(f"{pocket['rank']:<6} {pocket['name']:<10} {pocket['score']:>6.2f}  "
              f"{pocket['probability']:>6.3f}  {pocket['n_residues']:<10} {residues_preview}")

    print()

    # Analyze each pocket
    print("="*80)
    print("Active Site Overlap Analysis")
    print("="*80)
    print()

    best_f1 = 0
    best_pocket = None
    best_comparison = None

    for pocket in pockets:
        comparison = compare_with_known(pocket['residues'], KNOWN_ACTIVE_SITES)

        print(f"{pocket['name']} (Rank {pocket['rank']}, Score {pocket['score']:.2f}):")
        print(f"  Predicted: {comparison['n_predicted']} residues")
        print(f"  Overlap:   {comparison['n_correct']}/{comparison['n_known']} active sites")
        print(f"  Precision: {comparison['precision']:.1%}")
        print(f"  Recall:    {comparison['recall']:.1%}")
        print(f"  F1-score:  {comparison['f1_score']:.1%}")

        if comparison['n_correct'] > 0:
            print(f"  Found:     {comparison['true_positive'][:10]}")
            if len(comparison['true_positive']) > 10:
                print(f"             ... and {len(comparison['true_positive'])-10} more")

        print()

        if comparison['f1_score'] > best_f1:
            best_f1 = comparison['f1_score']
            best_pocket = pocket
            best_comparison = comparison

    # Summary
    print("="*80)
    print("Summary")
    print("="*80)
    print()

    if best_pocket:
        print(f"Best pocket: {best_pocket['name']} (Rank {best_pocket['rank']})")
        print(f"  Score: {best_pocket['score']:.2f}")
        print(f"  F1-score: {best_f1:.1%}")
        print(f"  Coverage: {best_comparison['n_correct']}/{len(KNOWN_ACTIVE_SITES)} active sites")
        print()

        print("Known active sites found in best pocket:")
        print(f"  {best_comparison['true_positive']}")
        print()

        print("Known active sites NOT found:")
        print(f"  {best_comparison['false_negative']}")
        print()

        # Combined analysis (top 2 pockets)
        if len(pockets) >= 2:
            combined_residues = set(pockets[0]['residues']) | set(pockets[1]['residues'])
            combined_comparison = compare_with_known(list(combined_residues), KNOWN_ACTIVE_SITES)

            print(f"Combined Top 2 Pockets:")
            print(f"  Total residues: {len(combined_residues)}")
            print(f"  Active sites found: {combined_comparison['n_correct']}/{len(KNOWN_ACTIVE_SITES)}")
            print(f"  F1-score: {combined_comparison['f1_score']:.1%}")
            print()

    print("="*80)
    print("Conclusion")
    print("="*80)

    if best_f1 > 0.6:
        print("[EXCELLENT] P2RANK accurately predicted the active site!")
    elif best_f1 > 0.4:
        print("[GOOD] P2RANK found the active site region")
    elif best_f1 > 0.2:
        print("[MODERATE] P2RANK partially identified the active site")
    else:
        print("[POOR] P2RANK did not identify the active site well")

    print()
    print(f"P2RANK is {'RELIABLE' if best_f1 > 0.5 else 'UNRELIABLE'} for AtUdh active site prediction")
    print("="*80)


if __name__ == "__main__":
    main()
