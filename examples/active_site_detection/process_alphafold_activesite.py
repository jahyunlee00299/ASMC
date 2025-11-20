#!/usr/bin/env python3
"""
Process AlphaFold UDH structures to identify active sites

This script:
1. Processes all AlphaFold models (both .pdb and .cif files)
2. Aligns each model to the 3RFV reference structure
3. Identifies active site residues based on spatial alignment
4. Selects the best model for each UDH variant
5. Generates comprehensive report

Usage:
    python process_alphafold_activesite.py

Author: ASMC
Date: 2025-11-20
"""

import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import json

# Fix Windows console encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

try:
    from Bio.PDB import PDBParser, MMCIFParser, Superimposer
    import numpy as np
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


class AlphaFoldActiveSitePredictor:
    """Predict active sites in AlphaFold structures using reference alignment"""

    def __init__(self, reference_pdb: str, reference_sites_file: str,
                 distance_cutoff: float = 4.0):
        """
        Initialize predictor

        Args:
            reference_pdb: Path to reference structure (3RFV)
            reference_sites_file: Active site definition file
            distance_cutoff: Distance cutoff for identifying corresponding residues
        """
        self.reference_pdb = reference_pdb
        self.distance_cutoff = distance_cutoff

        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)

        # Load reference structure
        self.ref_structure = self.pdb_parser.get_structure("reference", reference_pdb)

        # Load reference active sites
        self.ref_chain, self.ref_residues = self.load_active_sites(reference_sites_file)

        print(f"Loaded reference: {reference_pdb}")
        print(f"Reference chain: {self.ref_chain}")
        print(f"Reference active sites: {len(self.ref_residues)} residues")
        print(f"Active site positions: {self.ref_residues}")
        print()

    def load_active_sites(self, sites_file: str) -> Tuple[str, List[int]]:
        """Load active site residue information"""
        with open(sites_file, 'r', encoding='utf-8') as f:
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

    def get_ca_atoms(self, structure, chain_id: str,
                    residue_numbers: Optional[List[int]] = None):
        """Get CA atoms from specified residues"""
        atoms = []
        try:
            # Try to get the specified chain
            chains = list(structure[0].get_chains())

            # If specified chain not found, use first chain
            chain = None
            for c in chains:
                if c.id == chain_id:
                    chain = c
                    break

            if chain is None and len(chains) > 0:
                chain = chains[0]
                # print(f"  Warning: Chain {chain_id} not found, using chain {chain.id}")

            if chain is None:
                return atoms

            for residue in chain:
                # Skip hetero atoms
                if residue.id[0] != ' ':
                    continue

                res_num = residue.id[1]

                # Filter by residue numbers if specified
                if residue_numbers is not None and res_num not in residue_numbers:
                    continue

                # Get CA atom
                if 'CA' in residue:
                    atoms.append(residue['CA'])

        except Exception as e:
            print(f"  Warning: Error getting CA atoms: {e}")

        return atoms

    def align_structure(self, target_structure, target_name: str) -> Tuple[object, str, float]:
        """Align target structure to reference"""
        # Get first chain
        chains = list(target_structure[0].get_chains())
        if len(chains) == 0:
            raise ValueError(f"No chains found in {target_name}")

        target_chain = chains[0].id

        # Get CA atoms for alignment
        ref_ca_atoms = self.get_ca_atoms(self.ref_structure, self.ref_chain)
        target_ca_atoms = self.get_ca_atoms(target_structure, target_chain)

        if len(ref_ca_atoms) == 0 or len(target_ca_atoms) == 0:
            raise ValueError(f"No CA atoms found for alignment")

        # Use minimum number of atoms
        n_atoms = min(len(ref_ca_atoms), len(target_ca_atoms))
        ref_ca_atoms = ref_ca_atoms[:n_atoms]
        target_ca_atoms = target_ca_atoms[:n_atoms]

        # Perform superposition
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_ca_atoms, target_ca_atoms)
        super_imposer.apply(target_structure.get_atoms())

        rmsd = super_imposer.rms

        return target_structure, target_chain, rmsd

    def find_corresponding_residues(self, target_structure, target_chain: str,
                                   target_name: str) -> List[Tuple[int, str, float, int]]:
        """
        Find corresponding active site residues

        Returns:
            List of (target_res_num, res_name, distance, ref_res_num)
        """
        # Get reference active site CA atoms
        ref_active_atoms = self.get_ca_atoms(self.ref_structure, self.ref_chain,
                                             self.ref_residues)

        if len(ref_active_atoms) == 0:
            raise ValueError("No active site atoms found in reference")

        # Get target chain
        target_chain_obj = None
        for chain in target_structure[0].get_chains():
            if chain.id == target_chain:
                target_chain_obj = chain
                break

        if target_chain_obj is None:
            raise ValueError(f"Chain {target_chain} not found in target")

        corresponding_residues = []

        for ref_atom in ref_active_atoms:
            ref_coord = ref_atom.get_coord()
            ref_res = ref_atom.get_parent()
            ref_res_num = ref_res.id[1]

            min_dist = float('inf')
            closest_residue = None

            # Find closest residue in target
            for target_res in target_chain_obj:
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
                corresponding_residues.append((res_num, res_name, min_dist, ref_res_num))

        return corresponding_residues

    def process_structure_file(self, file_path: str) -> Optional[Dict]:
        """Process a single structure file (PDB or CIF)"""
        file_path = Path(file_path)
        file_name = file_path.name

        try:
            # Load structure based on file type
            if file_path.suffix.lower() == '.pdb':
                structure = self.pdb_parser.get_structure(file_name, str(file_path))
            elif file_path.suffix.lower() == '.cif':
                structure = self.cif_parser.get_structure(file_name, str(file_path))
            else:
                print(f"  Skipping unsupported file type: {file_name}")
                return None

            # Align structure
            aligned_structure, chain_id, rmsd = self.align_structure(structure, file_name)

            # Find corresponding residues
            corresponding = self.find_corresponding_residues(aligned_structure, chain_id,
                                                            file_name)

            result = {
                'file': file_name,
                'file_path': str(file_path),
                'chain': chain_id,
                'rmsd': rmsd,
                'n_active_sites': len(corresponding),
                'active_sites': corresponding,
                'coverage': len(corresponding) / len(self.ref_residues) * 100
            }

            return result

        except Exception as e:
            print(f"  Error processing {file_name}: {e}")
            return None

    def process_alphafold_directory(self, alphafold_dir: str, output_dir: str):
        """Process all AlphaFold structures"""
        alphafold_path = Path(alphafold_dir)
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        # Find all subdirectories (UDH variants)
        variant_dirs = [d for d in alphafold_path.iterdir() if d.is_dir()]

        print(f"Found {len(variant_dirs)} UDH variants to process")
        print("="*70)
        print()

        all_results = {}

        for variant_dir in sorted(variant_dirs):
            variant_name = variant_dir.name
            print(f"Processing {variant_name}...")
            print("-"*70)

            # Find all structure files in this variant
            pdb_files = list(variant_dir.glob("*.pdb"))
            cif_files = list(variant_dir.glob("*.cif"))
            all_files = sorted(pdb_files + cif_files)

            print(f"  Found {len(all_files)} structure files ({len(pdb_files)} PDB, {len(cif_files)} CIF)")

            variant_results = []

            for struct_file in all_files:
                print(f"    Processing {struct_file.name}...", end=' ')
                result = self.process_structure_file(struct_file)

                if result:
                    variant_results.append(result)
                    print(f"RMSD={result['rmsd']:.2f}Å, Sites={result['n_active_sites']}/{len(self.ref_residues)} ({result['coverage']:.1f}%)")
                else:
                    print("FAILED")

            # Select best model (lowest RMSD)
            if variant_results:
                best_model = min(variant_results, key=lambda x: x['rmsd'])
                all_results[variant_name] = {
                    'all_models': variant_results,
                    'best_model': best_model,
                    'n_models': len(variant_results)
                }

                print(f"\n  Best model: {best_model['file']}")
                print(f"    RMSD: {best_model['rmsd']:.2f} Å")
                print(f"    Active sites found: {best_model['n_active_sites']}/{len(self.ref_residues)} ({best_model['coverage']:.1f}%)")

                # Save individual variant results
                variant_output = output_path / f"{variant_name}_active_sites.txt"
                self.save_variant_results(variant_name, best_model, variant_output)

            print()

        # Save summary
        self.save_summary(all_results, output_path)

        # Save detailed JSON
        json_output = output_path / "all_results.json"
        with open(json_output, 'w', encoding='utf-8') as f:
            # Convert to serializable format
            serializable_results = {}
            for variant, data in all_results.items():
                serializable_results[variant] = {
                    'n_models': data['n_models'],
                    'best_model': {
                        'file': data['best_model']['file'],
                        'rmsd': float(data['best_model']['rmsd']),
                        'n_active_sites': data['best_model']['n_active_sites'],
                        'coverage': float(data['best_model']['coverage']),
                        'active_sites': [
                            {
                                'target_position': int(s[0]),
                                'residue_type': s[1],
                                'distance': float(s[2]),
                                'reference_position': int(s[3])
                            }
                            for s in data['best_model']['active_sites']
                        ]
                    }
                }
            json.dump(serializable_results, f, indent=2)

        print(f"\nDetailed JSON results saved to: {json_output}")

        return all_results

    def save_variant_results(self, variant_name: str, result: Dict, output_file: Path):
        """Save results for a single variant"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f"# Active sites predicted for {variant_name}\n")
            f.write(f"# Reference: {os.path.basename(self.reference_pdb)}\n")
            f.write(f"# Best model: {result['file']}\n")
            f.write(f"# Alignment RMSD: {result['rmsd']:.2f} Å\n")
            f.write(f"# Distance cutoff: {self.distance_cutoff} Å\n")
            f.write(f"# Chain: {result['chain']}\n")
            f.write(f"# Coverage: {result['n_active_sites']}/{len(self.ref_residues)} ({result['coverage']:.1f}%)\n")
            f.write(f"#\n")
            f.write(f"# Format: Target_Position[TAB]Residue[TAB]Distance[TAB]Reference_Position\n")
            f.write(f"#\n\n")

            # Write residue numbers
            residue_numbers = [str(s[0]) for s in result['active_sites']]
            f.write(f"{result['file']}\t{result['chain']}\t{','.join(residue_numbers)}\n\n")

            # Detailed information
            f.write(f"# Detailed mapping:\n")
            f.write(f"# Target_Pos\tResidue\tDistance(Å)\tRef_Pos\tRef_Residue\n")

            for target_pos, res_name, dist, ref_pos in result['active_sites']:
                # Get reference residue name
                ref_res_name = "?"
                for residue in self.ref_structure[0][self.ref_chain]:
                    if residue.id[0] == ' ' and residue.id[1] == ref_pos:
                        ref_res_name = residue.get_resname()
                        break

                f.write(f"{target_pos}\t{res_name}\t{dist:.2f}\t{ref_pos}\t{ref_res_name}\n")

    def save_summary(self, all_results: Dict, output_path: Path):
        """Save summary comparison file"""
        summary_file = output_path / "alphafold_activesite_summary.tsv"

        with open(summary_file, 'w', encoding='utf-8') as f:
            # Header
            f.write("Variant\tBest_Model\tRMSD(Å)\tN_Sites\tCoverage(%)\tResidue_Positions\n")

            # Reference
            ref_name = Path(self.reference_pdb).stem
            ref_res_str = ','.join(map(str, self.ref_residues))
            f.write(f"Reference_3RFV\t{ref_name}\t0.00\t{len(self.ref_residues)}\t100.0\t{ref_res_str}\n")

            # All variants
            for variant_name in sorted(all_results.keys()):
                data = all_results[variant_name]
                best = data['best_model']

                res_positions = ','.join(str(s[0]) for s in best['active_sites'])

                f.write(f"{variant_name}\t{best['file']}\t{best['rmsd']:.2f}\t"
                       f"{best['n_active_sites']}\t{best['coverage']:.1f}\t{res_positions}\n")

        print(f"Summary saved to: {summary_file}")


def main():
    """Main execution"""
    # Paths
    base_dir = Path(__file__).parent
    reference_pdb = base_dir / "pdb3rfv_chainA.pdb"
    reference_sites = base_dir / "atudh_active_site.txt"
    alphafold_dir = base_dir / "Alphafold server"
    output_dir = base_dir / "alphafold_activesite_results"

    # Validate inputs
    if not reference_pdb.exists():
        print(f"Error: Reference PDB not found: {reference_pdb}")
        sys.exit(1)

    if not reference_sites.exists():
        print(f"Error: Reference sites file not found: {reference_sites}")
        sys.exit(1)

    if not alphafold_dir.exists():
        print(f"Error: AlphaFold directory not found: {alphafold_dir}")
        sys.exit(1)

    print("="*70)
    print("AlphaFold UDH Active Site Prediction")
    print("="*70)
    print()

    # Initialize predictor
    predictor = AlphaFoldActiveSitePredictor(
        reference_pdb=str(reference_pdb),
        reference_sites_file=str(reference_sites),
        distance_cutoff=4.0
    )

    # Process all structures
    results = predictor.process_alphafold_directory(
        str(alphafold_dir),
        str(output_dir)
    )

    print()
    print("="*70)
    print(f"Processing complete! Results saved to: {output_dir}")
    print(f"Total variants processed: {len(results)}")
    print("="*70)


if __name__ == "__main__":
    main()
