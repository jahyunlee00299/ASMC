#!/usr/bin/env python3
"""
Find active sites in UDH structures using AtUdh as reference

This script aligns target UDH structures to AtUdh reference structure
and identifies corresponding active site residues based on spatial proximity.

Usage:
    python find_udh_active_sites.py --reference AtUdh_3rfv.pdb \
                                     --reference-sites atudh_active_site.txt \
                                     --target other_udh.pdb \
                                     --output results/

Author: ASMC
Date: 2025-11-20
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import subprocess

try:
    from Bio.PDB import PDBParser, Superimposer, PDBIO, Selection
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    import numpy as np
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


class UDHActiveSiteFinder:
    """Find active sites in UDH structures using reference structure alignment"""

    def __init__(self, reference_pdb: str, reference_sites_file: str,
                 distance_cutoff: float = 4.0):
        """
        Initialize the active site finder

        Args:
            reference_pdb: Path to reference AtUdh PDB file
            reference_sites_file: Path to file containing reference active site residues
            distance_cutoff: Distance cutoff (Angstroms) for identifying corresponding residues
        """
        self.reference_pdb = reference_pdb
        self.reference_sites_file = reference_sites_file
        self.distance_cutoff = distance_cutoff
        self.parser = PDBParser(QUIET=True)

        # Load reference structure
        self.ref_structure = self.parser.get_structure("reference", reference_pdb)

        # Load reference active site information
        self.ref_chain, self.ref_residues = self.load_active_sites(reference_sites_file)

    def load_active_sites(self, sites_file: str) -> Tuple[str, List[int]]:
        """
        Load active site residue information from file

        Args:
            sites_file: Path to active site definition file

        Returns:
            Tuple of (chain_id, list of residue numbers)
        """
        with open(sites_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 3:
                    # Format: pdb_file[TAB]chain[TAB]residue_numbers
                    chain = parts[1]
                    residues = [int(x.strip()) for x in parts[2].split(',')]
                    return chain, residues

        raise ValueError(f"No valid active site definition found in {sites_file}")

    def get_ca_atoms(self, structure, chain_id: str, residue_numbers: Optional[List[int]] = None):
        """
        Get CA atoms from specified residues

        Args:
            structure: BioPython structure object
            chain_id: Chain identifier
            residue_numbers: List of residue numbers (None = all residues)

        Returns:
            List of CA atoms
        """
        atoms = []
        try:
            chain = structure[0][chain_id]
            for residue in chain:
                # Skip hetero atoms (water, ligands, etc.)
                if residue.id[0] != ' ':
                    continue

                res_num = residue.id[1]

                # Filter by residue numbers if specified
                if residue_numbers is not None and res_num not in residue_numbers:
                    continue

                # Get CA atom
                if 'CA' in residue:
                    atoms.append(residue['CA'])

        except KeyError as e:
            print(f"Warning: Chain {chain_id} not found or incomplete")

        return atoms

    def align_structures(self, target_pdb: str, target_chain: str = None) -> Tuple[object, str, float]:
        """
        Align target structure to reference structure

        Args:
            target_pdb: Path to target PDB file
            target_chain: Target chain ID (auto-detect if None)

        Returns:
            Tuple of (aligned target structure, chain_id, rmsd)
        """
        # Load target structure
        target_structure = self.parser.get_structure("target", target_pdb)

        # Auto-detect chain if not specified
        if target_chain is None:
            chains = list(target_structure[0].get_chains())
            if len(chains) == 0:
                raise ValueError(f"No chains found in {target_pdb}")
            target_chain = chains[0].id
            print(f"Auto-detected chain: {target_chain}")

        # Get CA atoms for alignment (use all residues for global alignment)
        ref_ca_atoms = self.get_ca_atoms(self.ref_structure, self.ref_chain)
        target_ca_atoms = self.get_ca_atoms(target_structure, target_chain)

        if len(ref_ca_atoms) == 0 or len(target_ca_atoms) == 0:
            raise ValueError("No CA atoms found for alignment")

        # For alignment, use minimum number of atoms
        n_atoms = min(len(ref_ca_atoms), len(target_ca_atoms))
        ref_ca_atoms = ref_ca_atoms[:n_atoms]
        target_ca_atoms = target_ca_atoms[:n_atoms]

        print(f"Aligning {n_atoms} CA atoms...")

        # Perform superposition
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_ca_atoms, target_ca_atoms)

        # Apply transformation to entire target structure
        super_imposer.apply(target_structure.get_atoms())

        rmsd = super_imposer.rms
        print(f"Alignment RMSD: {rmsd:.2f} Å")

        return target_structure, target_chain, rmsd

    def find_corresponding_residues(self, target_structure, target_chain: str) -> List[Tuple[int, str, float]]:
        """
        Find residues in target structure corresponding to reference active sites

        Args:
            target_structure: Aligned target structure
            target_chain: Target chain ID

        Returns:
            List of tuples (residue_number, residue_name, min_distance)
        """
        # Get reference active site CA atoms
        ref_active_atoms = self.get_ca_atoms(self.ref_structure, self.ref_chain, self.ref_residues)

        if len(ref_active_atoms) == 0:
            raise ValueError("No active site atoms found in reference")

        # Get all target residues
        target_chain_obj = target_structure[0][target_chain]

        corresponding_residues = []

        for ref_atom in ref_active_atoms:
            ref_coord = ref_atom.get_coord()
            ref_res = ref_atom.get_parent()

            min_dist = float('inf')
            closest_residue = None

            # Find closest residue in target
            for target_res in target_chain_obj:
                # Skip hetero residues
                if target_res.id[0] != ' ':
                    continue

                if 'CA' not in target_res:
                    continue

                target_coord = target_res['CA'].get_coord()
                dist = np.linalg.norm(ref_coord - target_coord)

                if dist < min_dist:
                    min_dist = dist
                    closest_residue = target_res

            if closest_residue and min_dist <= self.distance_cutoff:
                res_num = closest_residue.id[1]
                res_name = closest_residue.get_resname()
                corresponding_residues.append((res_num, res_name, min_dist))
                print(f"  Reference {ref_res.get_resname()}{ref_res.id[1]} -> "
                      f"Target {res_name}{res_num} (distance: {min_dist:.2f} Å)")

        return corresponding_residues

    def process_target(self, target_pdb: str, target_chain: str = None,
                      output_dir: str = "results") -> Dict:
        """
        Process a target UDH structure

        Args:
            target_pdb: Path to target PDB file
            target_chain: Target chain ID (auto-detect if None)
            output_dir: Output directory for results

        Returns:
            Dictionary with results
        """
        print(f"\n{'='*60}")
        print(f"Processing: {target_pdb}")
        print(f"{'='*60}")

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Align structures
        aligned_structure, chain_id, rmsd = self.align_structures(target_pdb, target_chain)

        # Find corresponding residues
        print(f"\nFinding active sites (cutoff: {self.distance_cutoff} Å)...")
        corresponding = self.find_corresponding_residues(aligned_structure, chain_id)

        # Save results
        target_name = Path(target_pdb).stem
        output_file = os.path.join(output_dir, f"{target_name}_active_sites.txt")

        with open(output_file, 'w') as f:
            f.write(f"# Active sites found in {target_name}\n")
            f.write(f"# Reference: {os.path.basename(self.reference_pdb)}\n")
            f.write(f"# Alignment RMSD: {rmsd:.2f} Å\n")
            f.write(f"# Distance cutoff: {self.distance_cutoff} Å\n")
            f.write(f"# Chain: {chain_id}\n")
            f.write(f"#\n")
            f.write(f"# Format: PDB_file[TAB]Chain[TAB]Residue_numbers\n")
            f.write(f"#\n")

            residue_numbers = [str(r[0]) for r in corresponding]
            f.write(f"{os.path.basename(target_pdb)}\t{chain_id}\t{','.join(residue_numbers)}\n")

            f.write(f"\n# Detailed residue information:\n")
            f.write(f"# ResNum\tResName\tDistance(Å)\n")
            for res_num, res_name, dist in corresponding:
                f.write(f"# {res_num}\t{res_name}\t{dist:.2f}\n")

        print(f"\nResults saved to: {output_file}")
        print(f"Found {len(corresponding)} active site residues")

        return {
            'target': target_name,
            'chain': chain_id,
            'rmsd': rmsd,
            'active_sites': corresponding,
            'output_file': output_file
        }

    def process_multiple_targets(self, target_dir: str, output_dir: str = "results",
                                pattern: str = "*.pdb") -> List[Dict]:
        """
        Process multiple target structures

        Args:
            target_dir: Directory containing target PDB files
            output_dir: Output directory for results
            pattern: File pattern to match (default: *.pdb)

        Returns:
            List of result dictionaries
        """
        target_path = Path(target_dir)
        pdb_files = sorted(target_path.glob(pattern))

        if not pdb_files:
            print(f"No PDB files found in {target_dir} matching pattern {pattern}")
            return []

        print(f"Found {len(pdb_files)} PDB files to process")

        results = []
        for pdb_file in pdb_files:
            try:
                result = self.process_target(str(pdb_file), output_dir=output_dir)
                results.append(result)
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")
                continue

        # Create summary file
        self.create_summary(results, output_dir)

        return results

    def create_summary(self, results: List[Dict], output_dir: str):
        """Create summary file comparing all structures"""
        summary_file = os.path.join(output_dir, "active_sites_comparison.tsv")

        with open(summary_file, 'w') as f:
            # Header
            f.write("Structure\tChain\tRMSD(Å)\tN_Sites\tResidue_Positions\tResidue_Types\n")

            # Reference
            ref_name = Path(self.reference_pdb).stem
            ref_res_str = ','.join(map(str, self.ref_residues))
            f.write(f"{ref_name}\t{self.ref_chain}\t0.00\t{len(self.ref_residues)}\t{ref_res_str}\t-\n")

            # Targets
            for result in results:
                target = result['target']
                chain = result['chain']
                rmsd = result['rmsd']
                sites = result['active_sites']
                n_sites = len(sites)

                res_positions = ','.join(str(s[0]) for s in sites)
                res_types = ','.join(s[1] for s in sites)

                f.write(f"{target}\t{chain}\t{rmsd:.2f}\t{n_sites}\t{res_positions}\t{res_types}\n")

        print(f"\nSummary saved to: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Find active sites in UDH structures using AtUdh as reference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single target structure
  python find_udh_active_sites.py --reference AtUdh_3rfv.pdb \\
                                   --reference-sites atudh_active_site.txt \\
                                   --target other_udh.pdb

  # Process multiple structures in a directory
  python find_udh_active_sites.py --reference AtUdh_3rfv.pdb \\
                                   --reference-sites atudh_active_site.txt \\
                                   --target-dir other_udh_structures/ \\
                                   --output results/

  # Use custom distance cutoff
  python find_udh_active_sites.py --reference AtUdh_3rfv.pdb \\
                                   --reference-sites atudh_active_site.txt \\
                                   --target other_udh.pdb \\
                                   --cutoff 3.5
        """
    )

    parser.add_argument('--reference', '-r', required=True,
                       help='Reference AtUdh PDB file')
    parser.add_argument('--reference-sites', '-s', required=True,
                       help='Reference active site definition file')
    parser.add_argument('--target', '-t',
                       help='Target UDH PDB file (single structure)')
    parser.add_argument('--target-dir', '-d',
                       help='Directory containing target PDB files')
    parser.add_argument('--target-chain', '-c',
                       help='Target chain ID (auto-detect if not specified)')
    parser.add_argument('--output', '-o', default='results',
                       help='Output directory (default: results)')
    parser.add_argument('--cutoff', type=float, default=4.0,
                       help='Distance cutoff in Angstroms (default: 4.0)')
    parser.add_argument('--pattern', default='*.pdb',
                       help='File pattern for target directory (default: *.pdb)')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.reference):
        print(f"Error: Reference file not found: {args.reference}")
        sys.exit(1)

    if not os.path.exists(args.reference_sites):
        print(f"Error: Reference sites file not found: {args.reference_sites}")
        sys.exit(1)

    if not args.target and not args.target_dir:
        print("Error: Either --target or --target-dir must be specified")
        parser.print_help()
        sys.exit(1)

    # Initialize finder
    finder = UDHActiveSiteFinder(
        reference_pdb=args.reference,
        reference_sites_file=args.reference_sites,
        distance_cutoff=args.cutoff
    )

    # Process target(s)
    if args.target:
        # Single target
        if not os.path.exists(args.target):
            print(f"Error: Target file not found: {args.target}")
            sys.exit(1)

        finder.process_target(
            target_pdb=args.target,
            target_chain=args.target_chain,
            output_dir=args.output
        )
    else:
        # Multiple targets
        if not os.path.isdir(args.target_dir):
            print(f"Error: Target directory not found: {args.target_dir}")
            sys.exit(1)

        finder.process_multiple_targets(
            target_dir=args.target_dir,
            output_dir=args.output,
            pattern=args.pattern
        )

    print("\n" + "="*60)
    print("Processing complete!")
    print("="*60)


if __name__ == "__main__":
    main()
