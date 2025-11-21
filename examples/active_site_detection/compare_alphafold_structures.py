#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare experimental structure with AlphaFold2 and AlphaFold3 predictions

Compares:
1. 3RFV Chain A (experimental)
2. AlphaFold2 prediction (AF-Q7CRQ0)
3. AlphaFold3 prediction (REPORTED_AtUDH)

Analyzes:
- RMSD alignment
- Structural similarity
- Active site residue positions
- Predicted pockets (if available)

Author: ASMC
Date: 2025-11-20
"""

import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from pathlib import Path

try:
    from Bio.PDB import PDBParser, Superimposer, Selection
except ImportError:
    print("Error: BioPython required. Install: pip install biopython")
    sys.exit(1)


# Known AtUdh active sites from experimental structure
KNOWN_ACTIVE_SITES = [10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137,
                      139, 161, 163, 165, 187, 189, 211, 233, 235, 257, 259]


class StructureComparator:
    """Compare multiple protein structures"""

    def __init__(self, reference_pdb, chain_id='A'):
        self.parser = PDBParser(QUIET=True)
        self.reference_pdb = reference_pdb
        self.chain_id = chain_id

        # Load reference
        self.ref_structure = self.parser.get_structure("reference", reference_pdb)
        self.ref_chain = self.ref_structure[0][chain_id]

    def get_ca_atoms(self, structure, chain_id):
        """Extract CA atoms from structure"""
        atoms = []
        try:
            chain = structure[0][chain_id]
            for residue in chain:
                if residue.id[0] == ' ' and 'CA' in residue:
                    atoms.append(residue['CA'])
        except KeyError:
            print(f"Warning: Chain {chain_id} not found")
        return atoms

    def align_structures(self, target_pdb, target_chain='A'):
        """
        Align target structure to reference

        Returns:
            dict: {
                'rmsd': float,
                'aligned_structure': structure object,
                'n_atoms': int
            }
        """
        # Load target
        target_structure = self.parser.get_structure("target", target_pdb)

        # Get CA atoms
        ref_ca = self.get_ca_atoms(self.ref_structure, self.chain_id)
        target_ca = self.get_ca_atoms(target_structure, target_chain)

        if not ref_ca or not target_ca:
            return {'rmsd': None, 'error': 'No CA atoms found'}

        # Use minimum length
        n_atoms = min(len(ref_ca), len(target_ca))
        ref_ca = ref_ca[:n_atoms]
        target_ca = target_ca[:n_atoms]

        # Superimpose
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_ca, target_ca)
        super_imposer.apply(target_structure.get_atoms())

        return {
            'rmsd': super_imposer.rms,
            'aligned_structure': target_structure,
            'n_atoms': n_atoms
        }

    def get_residue_info(self, structure, chain_id='A'):
        """Get residue information from structure"""
        residues = []
        try:
            chain = structure[0][chain_id]
            for residue in chain:
                if residue.id[0] == ' ':
                    residues.append({
                        'number': residue.id[1],
                        'name': residue.get_resname(),
                        'ca_coord': residue['CA'].get_coord() if 'CA' in residue else None
                    })
        except KeyError:
            pass
        return residues

    def compare_active_sites(self, target_structure, target_chain='A',
                           known_sites=None, distance_cutoff=4.0):
        """
        Compare active site positions between reference and target

        Args:
            target_structure: Aligned target structure
            target_chain: Target chain ID
            known_sites: List of known active site residue numbers
            distance_cutoff: Distance cutoff for matching (Angstroms)

        Returns:
            dict: Comparison statistics
        """
        if known_sites is None:
            known_sites = KNOWN_ACTIVE_SITES

        # Get reference active site CA atoms
        ref_chain = self.ref_structure[0][self.chain_id]
        ref_active_coords = []

        for residue in ref_chain:
            if residue.id[0] == ' ' and residue.id[1] in known_sites:
                if 'CA' in residue:
                    ref_active_coords.append({
                        'number': residue.id[1],
                        'coord': residue['CA'].get_coord()
                    })

        # Find corresponding residues in target
        target_chain_obj = target_structure[0][target_chain]
        found_sites = []

        for ref_site in ref_active_coords:
            ref_coord = ref_site['coord']
            min_dist = float('inf')
            closest_res = None

            for target_res in target_chain_obj:
                if target_res.id[0] != ' ' or 'CA' not in target_res:
                    continue

                target_coord = target_res['CA'].get_coord()
                dist = np.linalg.norm(ref_coord - target_coord)

                if dist < min_dist:
                    min_dist = dist
                    closest_res = target_res

            if closest_res and min_dist <= distance_cutoff:
                found_sites.append({
                    'ref_number': ref_site['number'],
                    'target_number': closest_res.id[1],
                    'distance': min_dist
                })

        return {
            'n_known': len(known_sites),
            'n_found': len(found_sites),
            'found_sites': found_sites,
            'coverage': len(found_sites) / len(known_sites) if known_sites else 0
        }


def main():
    print("="*80)
    print("AtUdh Structure Comparison: Experimental vs AlphaFold2 vs AlphaFold3")
    print("="*80)
    print()

    # File paths
    experimental = "../../data/pdb_structures/pdb3rfv_chainA.pdb"
    alphafold2 = "../../data/pdb_structures/AF-Q7CRQ0-F1-model_v6.pdb"
    alphafold3 = "../../data/pdb_structures/REPORTED_AtUDH.pdb"

    structures = {
        'Experimental (3RFV)': experimental,
        'AlphaFold2': alphafold2,
        'AlphaFold3': alphafold3
    }

    # Initialize comparator with experimental structure as reference
    print("Loading experimental structure (3RFV Chain A) as reference...")
    comparator = StructureComparator(experimental, chain_id='A')

    # Get reference info
    ref_residues = comparator.get_residue_info(comparator.ref_structure, 'A')
    print(f"Reference structure: {len(ref_residues)} residues")
    print(f"Known active sites: {len(KNOWN_ACTIVE_SITES)} residues")
    print()

    # Compare structures
    results = {}

    print("="*80)
    print("RMSD Alignment Analysis")
    print("="*80)
    print()

    for name, pdb_file in structures.items():
        if name == 'Experimental (3RFV)':
            results[name] = {
                'rmsd': 0.0,
                'n_atoms': len(ref_residues),
                'active_site_coverage': 1.0,
                'n_found_sites': len(KNOWN_ACTIVE_SITES)
            }
            print(f"{name}:")
            print(f"  RMSD: 0.00 A (reference)")
            print(f"  Atoms aligned: {len(ref_residues)}")
            print(f"  Active site coverage: 100%")
            print()
            continue

        print(f"{name}:")

        # Align
        alignment = comparator.align_structures(pdb_file, target_chain='A')

        if 'error' in alignment:
            print(f"  Error: {alignment['error']}")
            results[name] = {'error': alignment['error']}
            print()
            continue

        rmsd = alignment['rmsd']
        n_atoms = alignment['n_atoms']

        print(f"  RMSD: {rmsd:.2f} A")
        print(f"  Atoms aligned: {n_atoms}")

        # Compare active sites
        active_comparison = comparator.compare_active_sites(
            alignment['aligned_structure'],
            target_chain='A',
            known_sites=KNOWN_ACTIVE_SITES,
            distance_cutoff=4.0
        )

        coverage = active_comparison['coverage']
        n_found = active_comparison['n_found']

        print(f"  Active site coverage: {coverage:.1%} ({n_found}/{len(KNOWN_ACTIVE_SITES)})")

        results[name] = {
            'rmsd': rmsd,
            'n_atoms': n_atoms,
            'active_site_coverage': coverage,
            'n_found_sites': n_found,
            'found_sites': active_comparison['found_sites']
        }

        print()

    # Summary comparison
    print("="*80)
    print("Summary Comparison")
    print("="*80)
    print()

    print(f"{'Structure':<25} {'RMSD (A)':<12} {'Aligned':<10} {'Active Sites':<15} {'Coverage'}")
    print("-"*80)

    for name, result in results.items():
        if 'error' in result:
            print(f"{name:<25} {'ERROR':<12} {'-':<10} {'-':<15} {'-'}")
        else:
            rmsd = result['rmsd']
            n_atoms = result['n_atoms']
            n_found = result['n_found_sites']
            coverage = result['active_site_coverage']

            print(f"{name:<25} {rmsd:>6.2f}      {n_atoms:<10} {n_found}/{len(KNOWN_ACTIVE_SITES):<13} {coverage:>6.1%}")

    print()

    # Detailed active site comparison
    print("="*80)
    print("Active Site Residue Mapping")
    print("="*80)
    print()

    for name, result in results.items():
        if name == 'Experimental (3RFV)' or 'error' in result:
            continue

        print(f"{name}:")
        print(f"  {'Ref Res':<10} {'Target Res':<12} {'Distance (A)'}")
        print(f"  {'-'*35}")

        found_sites = result.get('found_sites', [])
        for site in found_sites[:10]:  # Show first 10
            ref_num = site['ref_number']
            target_num = site['target_number']
            dist = site['distance']
            print(f"  {ref_num:<10} {target_num:<12} {dist:>6.2f}")

        if len(found_sites) > 10:
            print(f"  ... and {len(found_sites) - 10} more")

        print()

    # Conclusions
    print("="*80)
    print("Conclusions")
    print("="*80)
    print()

    if 'AlphaFold2' in results and 'rmsd' in results['AlphaFold2']:
        af2_rmsd = results['AlphaFold2']['rmsd']
        af2_coverage = results['AlphaFold2']['active_site_coverage']

        print(f"AlphaFold2 (v6):")
        print(f"  - RMSD: {af2_rmsd:.2f} A")
        print(f"  - Active site prediction: {af2_coverage:.1%}")

        if af2_rmsd < 2.0:
            print(f"  - Excellent structural accuracy")
        elif af2_rmsd < 5.0:
            print(f"  - Good structural accuracy")
        else:
            print(f"  - Moderate structural deviation")
        print()

    if 'AlphaFold3' in results and 'rmsd' in results['AlphaFold3']:
        af3_rmsd = results['AlphaFold3']['rmsd']
        af3_coverage = results['AlphaFold3']['active_site_coverage']

        print(f"AlphaFold3:")
        print(f"  - RMSD: {af3_rmsd:.2f} A")
        print(f"  - Active site prediction: {af3_coverage:.1%}")

        if af3_rmsd < 2.0:
            print(f"  - Excellent structural accuracy")
        elif af3_rmsd < 5.0:
            print(f"  - Good structural accuracy")
        else:
            print(f"  - Moderate structural deviation")
        print()

    # Compare AF2 vs AF3
    if ('AlphaFold2' in results and 'rmsd' in results['AlphaFold2'] and
        'AlphaFold3' in results and 'rmsd' in results['AlphaFold3']):

        af2_rmsd = results['AlphaFold2']['rmsd']
        af3_rmsd = results['AlphaFold3']['rmsd']

        print("AlphaFold2 vs AlphaFold3:")

        if abs(af2_rmsd - af3_rmsd) < 0.5:
            print("  - Similar accuracy")
        elif af2_rmsd < af3_rmsd:
            print(f"  - AlphaFold2 is more accurate ({af2_rmsd:.2f} vs {af3_rmsd:.2f} A)")
        else:
            print(f"  - AlphaFold3 is more accurate ({af3_rmsd:.2f} vs {af2_rmsd:.2f} A)")

    print()
    print("="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == "__main__":
    main()
