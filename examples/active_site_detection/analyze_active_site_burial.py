#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze burial depth of known active site residues
to understand what parameters we should use for prediction
"""

import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

# Known active sites
known_active_sites = [10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137,
                     139, 161, 163, 165, 187, 189, 211, 233, 235, 257, 259]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "pdb3rfv_chainA.pdb")

chain = structure[0]['A']
all_atoms = list(structure.get_atoms())
ns = NeighborSearch(all_atoms)

print("Active Site Residue Burial Analysis")
print("="*70)
print(f"{'ResNum':<8} {'ResName':<8} {'Neighbors(10A)':<15} {'Status'}")
print("-"*70)

burial_counts = []

for residue in chain:
    if residue.id[0] != ' ':
        continue
    if 'CA' not in residue:
        continue

    res_num = residue.id[1]
    res_name = residue.get_resname()
    ca_atom = residue['CA']

    neighbors = ns.search(ca_atom.coord, 10.0)
    neighbor_count = len(neighbors)

    is_active = res_num in known_active_sites

    if is_active:
        burial_counts.append(neighbor_count)
        status = "← ACTIVE SITE"
        print(f"{res_num:<8} {res_name:<8} {neighbor_count:<15} {status}")

print("="*70)
print(f"\nStatistics for {len(burial_counts)} active site residues:")
print(f"  Min neighbors:    {min(burial_counts)}")
print(f"  Max neighbors:    {max(burial_counts)}")
print(f"  Mean neighbors:   {np.mean(burial_counts):.1f}")
print(f"  Median neighbors: {np.median(burial_counts):.1f}")
print(f"  Std dev:          {np.std(burial_counts):.1f}")

print(f"\nRecommended range: {int(min(burial_counts) * 0.8)} - {int(max(burial_counts) * 1.2)}")
