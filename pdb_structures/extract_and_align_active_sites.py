#!/usr/bin/env python3
"""
Extract and superpose active sites from UDH structures

This script:
1. Extracts active site residues from reference AtUdh and target UDH structures
2. Performs superposition using ONLY active site residues
3. Calculates active site RMSD
4. Generates alignment results and statistics

Target structures (sorted by Rossmann RMSD):
- A0A1I2LZE5 (15.78 Å) - Most similar
- A0A1I1AQL9 (17.27 Å)
- A0A2Z4AB20 (19.85 Å)
- A0A1N6JJV2 (20.89 Å)
- A0A1U7CQA8 (25.19 Å) - Most different

Author: ASMC
Date: 2025-11-20
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np

try:
    from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

# Target UDH structures to analyze
TARGET_STRUCTURES = [
    ('A0A1I2LZE5', 15.78, 'Most similar'),
    ('A0A1I1AQL9', 17.27, ''),
    ('A0A2Z4AB20', 19.85, ''),
    ('A0A1N6JJV2', 20.89, ''),
    ('A0A1U7CQA8', 25.19, 'Most different'),
]


class ActiveSiteSelector(Select):
    """BioPython Select class to save only active site residues"""

    def __init__(self, chain_id: str, residue_numbers: List[int]):
        self.chain_id = chain_id
        self.residue_numbers = set(residue_numbers)

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # Only accept standard residues in active site
        return (residue.id[0] == ' ' and
                residue.id[1] in self.residue_numbers)


class ActiveSiteExtractorAligner:
    """Extract and align active sites from UDH structures"""

    def __init__(self, reference_pdb: str, reference_sites_file: str):
        """
        Initialize the active site extractor and aligner

        Args:
            reference_pdb: Path to reference AtUdh PDB file
            reference_sites_file: Path to file containing reference active site residues
        """
        self.reference_pdb = reference_pdb
        self.reference_sites_file = reference_sites_file
        self.parser = PDBParser(QUIET=True)

        # Load reference structure
        self.ref_structure = self.parser.get_structure("reference", reference_pdb)

        # Load reference active site information
        self.ref_chain, self.ref_residues = self.load_active_sites(reference_sites_file)

        print(f"Reference: {os.path.basename(reference_pdb)}")
        print(f"Reference chain: {self.ref_chain}")
        print(f"Reference active sites: {len(self.ref_residues)} residues")
        print(f"Residue numbers: {','.join(map(str, sorted(self.ref_residues)))}\n")

    def load_active_sites(self, sites_file: str) -> Tuple[str, List[int]]:
        """Load active site residue information from file"""
        with open(sites_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 3:
                    chain = parts[1]
                    residues = [int(x.strip()) for x in parts[2].split(',')]
                    return chain, residues

        raise ValueError(f"No valid active site definition found in {sites_file}")

    def get_active_site_atoms(self, structure, chain_id: str,
                              residue_numbers: List[int],
                              atom_type: str = 'CA'):
        """
        Get specific atoms from active site residues

        Args:
            structure: BioPython structure
            chain_id: Chain identifier
            residue_numbers: List of residue numbers
            atom_type: Atom type to extract (default: CA for C-alpha)

        Returns:
            List of atoms, List of residue info (num, name)
        """
        atoms = []
        residue_info = []

        try:
            chain = structure[0][chain_id]

            for residue in chain:
                if residue.id[0] != ' ':  # Skip hetero atoms
                    continue

                res_num = residue.id[1]
                if res_num not in residue_numbers:
                    continue

                if atom_type in residue:
                    atoms.append(residue[atom_type])
                    residue_info.append((res_num, residue.get_resname()))

        except KeyError as e:
            print(f"Warning: Chain {chain_id} not found or incomplete: {e}")

        return atoms, residue_info

    def align_by_active_sites(self, target_pdb: str, target_chain: str = None,
                             distance_cutoff: float = 4.0) -> Dict:
        """
        Align target structure to reference using ONLY active site residues

        Args:
            target_pdb: Path to target PDB file
            target_chain: Target chain ID (auto-detect if None)
            distance_cutoff: Distance cutoff for finding corresponding residues

        Returns:
            Dictionary with alignment results
        """
        target_name = Path(target_pdb).stem
        print(f"\n{'='*70}")
        print(f"Processing: {target_name}")
        print(f"{'='*70}")

        # Load target structure
        target_structure = self.parser.get_structure("target", target_pdb)

        # Auto-detect chain
        if target_chain is None:
            chains = list(target_structure[0].get_chains())
            if len(chains) == 0:
                raise ValueError(f"No chains found in {target_pdb}")
            target_chain = chains[0].id
            print(f"Auto-detected chain: {target_chain}")

        # Step 1: Find corresponding active site residues in target
        # First, do a global alignment to get approximate positions
        ref_all_ca = []
        target_all_ca = []

        for res in self.ref_structure[0][self.ref_chain]:
            if res.id[0] == ' ' and 'CA' in res:
                ref_all_ca.append(res['CA'])

        for res in target_structure[0][target_chain]:
            if res.id[0] == ' ' and 'CA' in res:
                target_all_ca.append(res['CA'])

        # Align using minimum number of residues
        n_atoms = min(len(ref_all_ca), len(target_all_ca))
        print(f"Initial alignment using {n_atoms} CA atoms...")

        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_all_ca[:n_atoms], target_all_ca[:n_atoms])
        super_imposer.apply(target_structure.get_atoms())

        initial_rmsd = super_imposer.rms
        print(f"Initial global RMSD: {initial_rmsd:.2f} Å")

        # Step 2: Find corresponding active site residues in target
        ref_active_atoms, ref_active_info = self.get_active_site_atoms(
            self.ref_structure, self.ref_chain, self.ref_residues
        )

        print(f"\nFinding corresponding active sites (cutoff: {distance_cutoff} Å)...")

        target_active_atoms = []
        target_active_info = []
        correspondences = []

        target_chain_obj = target_structure[0][target_chain]

        for i, ref_atom in enumerate(ref_active_atoms):
            ref_coord = ref_atom.get_coord()
            ref_num, ref_name = ref_active_info[i]

            min_dist = float('inf')
            closest_atom = None
            closest_res_info = None

            for target_res in target_chain_obj:
                if target_res.id[0] != ' ' or 'CA' not in target_res:
                    continue

                target_coord = target_res['CA'].get_coord()
                dist = np.linalg.norm(ref_coord - target_coord)

                if dist < min_dist:
                    min_dist = dist
                    closest_atom = target_res['CA']
                    closest_res_info = (target_res.id[1], target_res.get_resname())

            if closest_atom and min_dist <= distance_cutoff:
                target_active_atoms.append(closest_atom)
                target_active_info.append(closest_res_info)
                correspondences.append({
                    'ref_num': ref_num,
                    'ref_name': ref_name,
                    'target_num': closest_res_info[0],
                    'target_name': closest_res_info[1],
                    'distance': min_dist
                })
                print(f"  {ref_name}{ref_num:3d} -> {closest_res_info[1]}{closest_res_info[0]:3d} "
                      f"(distance: {min_dist:.2f} Å)")

        print(f"\nFound {len(target_active_atoms)}/{len(ref_active_atoms)} "
              f"corresponding active site residues")

        # Step 3: Align using ONLY active site residues
        if len(target_active_atoms) < 3:
            print("Warning: Too few active site residues for alignment")
            return {
                'target': target_name,
                'chain': target_chain,
                'n_active_sites': len(target_active_atoms),
                'global_rmsd': initial_rmsd,
                'active_site_rmsd': None,
                'correspondences': correspondences,
                'success': False
            }

        print(f"\nAligning structures using {len(target_active_atoms)} active site CA atoms...")

        active_super = Superimposer()
        active_super.set_atoms(ref_active_atoms[:len(target_active_atoms)],
                              target_active_atoms)

        # Apply transformation to entire target structure
        active_super.apply(target_structure.get_atoms())

        active_site_rmsd = active_super.rms

        print(f"Active site RMSD: {active_site_rmsd:.2f} Å")
        print(f"Improvement: {initial_rmsd - active_site_rmsd:.2f} Å")

        return {
            'target': target_name,
            'chain': target_chain,
            'n_active_sites': len(target_active_atoms),
            'ref_active_sites': len(ref_active_atoms),
            'global_rmsd': initial_rmsd,
            'active_site_rmsd': active_site_rmsd,
            'improvement': initial_rmsd - active_site_rmsd,
            'correspondences': correspondences,
            'aligned_structure': target_structure,
            'success': True
        }

    def process_multiple_targets(self, target_dir: str,
                                 target_ids: List[str],
                                 output_dir: str = "active_site_results",
                                 pattern: str = "*.pdb") -> List[Dict]:
        """
        Process multiple target structures

        Args:
            target_dir: Directory containing target PDB files
            target_ids: List of target structure IDs to process
            output_dir: Output directory
            pattern: File pattern

        Returns:
            List of result dictionaries
        """
        os.makedirs(output_dir, exist_ok=True)

        results = []

        for struct_id, rossmann_rmsd, note in TARGET_STRUCTURES:
            # Find matching PDB file
            pdb_files = list(Path(target_dir).glob(f"*{struct_id}*.pdb"))

            if not pdb_files:
                print(f"\nWarning: No PDB file found for {struct_id}")
                continue

            pdb_file = pdb_files[0]

            try:
                result = self.align_by_active_sites(str(pdb_file))
                result['struct_id'] = struct_id
                result['rossmann_rmsd'] = rossmann_rmsd
                result['note'] = note
                results.append(result)

                # Save aligned structure
                if result['success'] and result['aligned_structure']:
                    self.save_aligned_structure(
                        result['aligned_structure'],
                        struct_id,
                        output_dir
                    )

                # Save individual result
                self.save_individual_result(result, output_dir)

            except Exception as e:
                print(f"Error processing {struct_id}: {e}")
                continue

        # Create summary
        if results:
            self.create_summary(results, output_dir)

        return results

    def save_aligned_structure(self, structure, struct_id: str, output_dir: str):
        """Save aligned structure to PDB file"""
        output_file = os.path.join(output_dir, f"{struct_id}_aligned.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        print(f"Saved aligned structure: {output_file}")

    def save_individual_result(self, result: Dict, output_dir: str):
        """Save individual alignment result"""
        output_file = os.path.join(output_dir,
                                   f"{result['struct_id']}_active_site_alignment.txt")

        with open(output_file, 'w') as f:
            f.write(f"# Active Site Alignment Results\n")
            f.write(f"# Target: {result['struct_id']}\n")
            f.write(f"# Rossmann RMSD: {result.get('rossmann_rmsd', 'N/A')} Å\n")
            f.write(f"# Note: {result.get('note', '')}\n")
            f.write(f"#\n")
            f.write(f"# Chain: {result['chain']}\n")
            f.write(f"# Active sites found: {result['n_active_sites']}/{result.get('ref_active_sites', 'N/A')}\n")
            f.write(f"# Global RMSD: {result['global_rmsd']:.2f} Å\n")
            f.write(f"# Active site RMSD: {result['active_site_rmsd']:.2f} Å\n")
            f.write(f"# Improvement: {result.get('improvement', 0):.2f} Å\n")
            f.write(f"#\n")
            f.write(f"# Correspondences:\n")
            f.write(f"# Ref_Res\tTarget_Res\tDistance(Å)\n")

            for corr in result['correspondences']:
                f.write(f"{corr['ref_name']}{corr['ref_num']}\t"
                       f"{corr['target_name']}{corr['target_num']}\t"
                       f"{corr['distance']:.2f}\n")

        print(f"Saved result: {output_file}")

    def create_summary(self, results: List[Dict], output_dir: str):
        """Create summary file with all results"""
        summary_file = os.path.join(output_dir, "active_site_alignment_summary.tsv")

        with open(summary_file, 'w') as f:
            # Header
            f.write("Structure_ID\tRossmann_RMSD(Å)\tNote\tChain\t"
                   "N_Active_Sites\tGlobal_RMSD(Å)\tActive_Site_RMSD(Å)\t"
                   "Improvement(Å)\tSuccess\n")

            # Data
            for result in results:
                f.write(f"{result['struct_id']}\t")
                f.write(f"{result.get('rossmann_rmsd', 'N/A')}\t")
                f.write(f"{result.get('note', '')}\t")
                f.write(f"{result['chain']}\t")
                f.write(f"{result['n_active_sites']}/{result.get('ref_active_sites', 'N/A')}\t")
                f.write(f"{result['global_rmsd']:.2f}\t")

                if result['active_site_rmsd'] is not None:
                    f.write(f"{result['active_site_rmsd']:.2f}\t")
                    f.write(f"{result.get('improvement', 0):.2f}\t")
                else:
                    f.write("N/A\tN/A\t")

                f.write(f"{result['success']}\n")

        print(f"\n{'='*70}")
        print(f"Summary saved: {summary_file}")
        print(f"{'='*70}")

        # Print summary table
        print("\nAlignment Summary:")
        print(f"{'Structure':<12} {'Rossmann':<10} {'Active Site':<15} {'Global':<10} {'Note'}")
        print(f"{'ID':<12} {'RMSD (Å)':<10} {'RMSD (Å)':<15} {'RMSD (Å)':<10}")
        print("-" * 70)

        for result in results:
            struct_id = result['struct_id']
            ross_rmsd = result.get('rossmann_rmsd', 'N/A')
            active_rmsd = f"{result['active_site_rmsd']:.2f}" if result['active_site_rmsd'] else "N/A"
            global_rmsd = f"{result['global_rmsd']:.2f}"
            note = result.get('note', '')

            print(f"{struct_id:<12} {ross_rmsd:<10} {active_rmsd:<15} {global_rmsd:<10} {note}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract and align active sites from UDH structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Target structures (sorted by Rossmann RMSD):
  A0A1I2LZE5 (15.78 Å) - Most similar
  A0A1I1AQL9 (17.27 Å)
  A0A2Z4AB20 (19.85 Å)
  A0A1N6JJV2 (20.89 Å)
  A0A1U7CQA8 (25.19 Å) - Most different

Example:
  python extract_and_align_active_sites.py \\
      --reference AtUdh_3rfv.pdb \\
      --reference-sites atudh_active_site.txt \\
      --target-dir uronate_dehydrogenase_alphafold/ \\
      --output active_site_results/
        """
    )

    parser.add_argument('--reference', '-r', required=True,
                       help='Reference AtUdh PDB file')
    parser.add_argument('--reference-sites', '-s', required=True,
                       help='Reference active site definition file')
    parser.add_argument('--target-dir', '-d', required=True,
                       help='Directory containing target UDH PDB files')
    parser.add_argument('--output', '-o', default='active_site_results',
                       help='Output directory (default: active_site_results)')
    parser.add_argument('--cutoff', type=float, default=4.0,
                       help='Distance cutoff in Angstroms (default: 4.0)')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}")
        sys.exit(1)

    if not os.path.exists(args.reference_sites):
        print(f"Error: Reference sites file not found: {args.reference_sites}")
        sys.exit(1)

    if not os.path.isdir(args.target_dir):
        print(f"Error: Target directory not found: {args.target_dir}")
        sys.exit(1)

    # Initialize extractor/aligner
    extractor = ActiveSiteExtractorAligner(
        reference_pdb=args.reference,
        reference_sites_file=args.reference_sites
    )

    # Process targets
    target_ids = [struct_id for struct_id, _, _ in TARGET_STRUCTURES]

    results = extractor.process_multiple_targets(
        target_dir=args.target_dir,
        target_ids=target_ids,
        output_dir=args.output
    )

    print(f"\n{'='*70}")
    print(f"Processing complete!")
    print(f"Processed {len(results)}/{len(TARGET_STRUCTURES)} structures")
    print(f"Results saved to: {args.output}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
