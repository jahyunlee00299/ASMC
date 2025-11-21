#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Active Site Extractor for Protein Structures

This module provides a unified interface for extracting active site residues
from protein structures using structural alignment and sequence alignment methods.

Author: ASMC
Date: 2025-11-21
"""

import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import warnings

# Fix Windows console encoding for Unicode characters
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

try:
    from Bio import SeqIO, pairwise2
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.PDB import PDBParser, Superimposer, PDBIO, Selection
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    import numpy as np
except ImportError as e:
    print(f"Error: Required BioPython package not found. Install with: pip install biopython")
    print(f"Details: {e}")
    sys.exit(1)


class ActiveSiteExtractor:
    """
    Unified class for extracting active site residues from protein structures.

    This class combines two approaches:
    1. Structural alignment-based extraction (for PDB structures)
    2. Sequence alignment-based extraction (for FASTA sequences)

    Features:
    - Extract active sites from single or multiple structures
    - Support for both PDB structures and FASTA sequences
    - Customizable distance cutoffs for structural matching
    - Automatic chain detection
    - Batch processing capabilities
    - Comprehensive result reporting
    """

    # Standard amino acid conversion
    AA_3TO1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    def __init__(self,
                 reference_pdb: Optional[str] = None,
                 reference_sites: Optional[Union[str, List[int]]] = None,
                 reference_chain: str = 'A',
                 distance_cutoff: float = 4.0,
                 verbose: bool = True):
        """
        Initialize the Active Site Extractor.

        Args:
            reference_pdb: Path to reference PDB file
            reference_sites: Either a file path or list of residue numbers defining active sites
            reference_chain: Chain ID in reference structure (default: 'A')
            distance_cutoff: Distance cutoff in Angstroms for structural matching (default: 4.0)
            verbose: Print detailed progress information
        """
        self.reference_pdb = reference_pdb
        self.reference_chain = reference_chain
        self.distance_cutoff = distance_cutoff
        self.verbose = verbose

        # BioPython PDB parser
        self.parser = PDBParser(QUIET=not verbose)

        # Reference structure data
        self.ref_structure = None
        self.ref_sequence = None
        self.ref_residue_numbers = None
        self.active_site_positions = None

        # Load reference if provided
        if reference_pdb:
            self._load_reference_structure()

        if reference_sites:
            self._load_active_sites(reference_sites)

    def _log(self, message: str):
        """Print message if verbose mode is enabled"""
        if self.verbose:
            print(message)

    def _load_reference_structure(self):
        """Load reference PDB structure"""
        if not os.path.exists(self.reference_pdb):
            raise FileNotFoundError(f"Reference PDB file not found: {self.reference_pdb}")

        self._log(f"Loading reference structure: {self.reference_pdb}")
        self.ref_structure = self.parser.get_structure("reference", self.reference_pdb)

        # Extract sequence
        self.ref_sequence, self.ref_residue_numbers = self.read_pdb_sequence(
            self.reference_pdb,
            chain_id=self.reference_chain
        )
        self._log(f"  Reference sequence length: {len(self.ref_sequence)}")

    def _load_active_sites(self, sites: Union[str, List[int]]):
        """
        Load active site residue positions.

        Args:
            sites: Either a file path (str) or list of residue numbers (List[int])
        """
        if isinstance(sites, list):
            # Direct list of residue numbers
            self.active_site_positions = sorted(sites)
            self._log(f"Loaded {len(self.active_site_positions)} active site positions")
        elif isinstance(sites, str):
            # File path
            if not os.path.exists(sites):
                raise FileNotFoundError(f"Active sites file not found: {sites}")

            self.active_site_positions = self._parse_sites_file(sites)
            self._log(f"Loaded {len(self.active_site_positions)} active site positions from {sites}")
        else:
            raise ValueError("sites must be either a file path (str) or list of residue numbers")

    def _parse_sites_file(self, sites_file: str) -> List[int]:
        """
        Parse active site definition file.

        File format:
        - Tab-delimited: pdb_file[TAB]chain[TAB]comma-separated residue numbers
        - Lines starting with '#' are comments

        Args:
            sites_file: Path to sites definition file

        Returns:
            List of residue numbers
        """
        with open(sites_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 3:
                    residues = [int(x.strip()) for x in parts[2].split(',')]
                    return sorted(residues)

        raise ValueError(f"No valid active site definition found in {sites_file}")

    @staticmethod
    def read_pdb_sequence(pdb_file: str, chain_id: str = 'A') -> Tuple[str, List[int]]:
        """
        Extract amino acid sequence from PDB file.

        Args:
            pdb_file: Path to PDB file
            chain_id: Chain identifier to extract

        Returns:
            Tuple of (sequence string, list of residue numbers)
        """
        residues = {}

        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain = line[21]
                    resnum = int(line[22:26].strip())
                    resname = line[17:20].strip()

                    if chain == chain_id and resnum not in residues:
                        if resname in ActiveSiteExtractor.AA_3TO1:
                            residues[resnum] = ActiveSiteExtractor.AA_3TO1[resname]

        # Sort by residue number
        sorted_residues = sorted(residues.items())
        sequence = ''.join([aa for _, aa in sorted_residues])
        residue_numbers = [num for num, _ in sorted_residues]

        return sequence, residue_numbers

    @staticmethod
    def align_sequences(seq1: str, seq2: str) -> Tuple[str, str, float]:
        """
        Perform global sequence alignment.

        Args:
            seq1: First sequence (reference)
            seq2: Second sequence (target)

        Returns:
            Tuple of (aligned_seq1, aligned_seq2, alignment_score)
        """
        alignments = pairwise2.align.globalxx(seq1, seq2)

        if alignments:
            best = alignments[0]
            return best.seqA, best.seqB, best.score
        else:
            return None, None, 0

    def extract_sites_from_alignment(self,
                                     aligned_ref: str,
                                     aligned_target: str,
                                     ref_residue_numbers: List[int],
                                     active_positions: List[int]) -> str:
        """
        Extract active site residues from aligned sequences.

        Args:
            aligned_ref: Aligned reference sequence
            aligned_target: Aligned target sequence
            ref_residue_numbers: Residue numbers in reference
            active_positions: Active site positions in reference numbering

        Returns:
            String of active site residues from target
        """
        # Map alignment positions to reference residue numbers
        ref_pos_to_resnum = {}
        ref_idx = 0

        for aln_idx, char in enumerate(aligned_ref):
            if char != '-':
                ref_pos_to_resnum[aln_idx] = ref_residue_numbers[ref_idx]
                ref_idx += 1

        # Extract active site residues from target
        active_site_residues = []
        for aln_idx, resnum in ref_pos_to_resnum.items():
            if resnum in active_positions:
                target_aa = aligned_target[aln_idx]
                active_site_residues.append(target_aa if target_aa != '-' else 'X')

        return ''.join(active_site_residues)

    def get_ca_atoms(self, structure, chain_id: str,
                     residue_numbers: Optional[List[int]] = None) -> List:
        """
        Get CA atoms from specified residues in a structure.

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

        except KeyError:
            warnings.warn(f"Chain {chain_id} not found or incomplete")

        return atoms

    def align_structures(self,
                        target_pdb: str,
                        target_chain: Optional[str] = None) -> Tuple[object, str, float]:
        """
        Align target structure to reference structure using superposition.

        Args:
            target_pdb: Path to target PDB file
            target_chain: Target chain ID (auto-detect if None)

        Returns:
            Tuple of (aligned target structure, chain_id, RMSD)
        """
        # Load target structure
        target_structure = self.parser.get_structure("target", target_pdb)

        # Auto-detect chain if not specified
        if target_chain is None:
            chains = list(target_structure[0].get_chains())
            if len(chains) == 0:
                raise ValueError(f"No chains found in {target_pdb}")
            target_chain = chains[0].id
            self._log(f"Auto-detected chain: {target_chain}")

        # Get CA atoms for alignment
        ref_ca_atoms = self.get_ca_atoms(self.ref_structure, self.reference_chain)
        target_ca_atoms = self.get_ca_atoms(target_structure, target_chain)

        if len(ref_ca_atoms) == 0 or len(target_ca_atoms) == 0:
            raise ValueError("No CA atoms found for alignment")

        # Use minimum number of atoms
        n_atoms = min(len(ref_ca_atoms), len(target_ca_atoms))
        ref_ca_atoms = ref_ca_atoms[:n_atoms]
        target_ca_atoms = target_ca_atoms[:n_atoms]

        self._log(f"Aligning {n_atoms} CA atoms...")

        # Perform superposition
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_ca_atoms, target_ca_atoms)
        super_imposer.apply(target_structure.get_atoms())

        rmsd = super_imposer.rms
        self._log(f"Alignment RMSD: {rmsd:.2f} Å")

        return target_structure, target_chain, rmsd

    def find_corresponding_residues(self,
                                   target_structure,
                                   target_chain: str) -> List[Tuple[int, str, float]]:
        """
        Find residues in target corresponding to reference active sites.

        Uses spatial proximity after structural alignment.

        Args:
            target_structure: Aligned target structure
            target_chain: Target chain ID

        Returns:
            List of tuples (residue_number, residue_name, min_distance)
        """
        # Get reference active site CA atoms
        ref_active_atoms = self.get_ca_atoms(
            self.ref_structure,
            self.reference_chain,
            self.active_site_positions
        )

        if len(ref_active_atoms) == 0:
            raise ValueError("No active site atoms found in reference")

        # Get target chain
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
                self._log(f"  Ref {ref_res.get_resname()}{ref_res.id[1]} -> "
                         f"Target {res_name}{res_num} (dist: {min_dist:.2f} Å)")

        return corresponding_residues

    def extract_from_structure(self,
                              target_pdb: str,
                              target_chain: Optional[str] = None,
                              output_file: Optional[str] = None) -> Dict:
        """
        Extract active sites from a PDB structure using structural alignment.

        Args:
            target_pdb: Path to target PDB file
            target_chain: Target chain ID (auto-detect if None)
            output_file: Optional output file for results

        Returns:
            Dictionary with extraction results
        """
        if not self.ref_structure:
            raise ValueError("Reference structure not loaded. Set reference_pdb in constructor.")

        if not self.active_site_positions:
            raise ValueError("Active site positions not defined. Set reference_sites in constructor.")

        self._log(f"\n{'='*60}")
        self._log(f"Processing structure: {target_pdb}")
        self._log(f"{'='*60}")

        # Align structures
        aligned_structure, chain_id, rmsd = self.align_structures(target_pdb, target_chain)

        # Find corresponding residues
        self._log(f"\nFinding active sites (cutoff: {self.distance_cutoff} Å)...")
        corresponding = self.find_corresponding_residues(aligned_structure, chain_id)

        # Prepare results
        target_name = Path(target_pdb).stem
        results = {
            'target': target_name,
            'chain': chain_id,
            'rmsd': rmsd,
            'n_sites': len(corresponding),
            'active_sites': corresponding,
            'residue_numbers': [r[0] for r in corresponding],
            'residue_types': [r[1] for r in corresponding],
            'distances': [r[2] for r in corresponding]
        }

        # Save to file if requested
        if output_file:
            self._write_structure_results(results, output_file, target_pdb)
            results['output_file'] = output_file

        self._log(f"\nFound {len(corresponding)} active site residues")

        return results

    def extract_from_sequences(self,
                              fasta_file: str,
                              output_file: Optional[str] = None,
                              progress_interval: int = 100) -> List[SeqRecord]:
        """
        Extract active sites from sequences using sequence alignment.

        Args:
            fasta_file: Path to FASTA file with target sequences
            output_file: Optional output FASTA file for active site sequences
            progress_interval: Show progress every N sequences

        Returns:
            List of SeqRecord objects with active site sequences
        """
        if not self.ref_sequence:
            raise ValueError("Reference sequence not loaded. Set reference_pdb in constructor.")

        if not self.active_site_positions:
            raise ValueError("Active site positions not defined. Set reference_sites in constructor.")

        self._log(f"\n{'='*60}")
        self._log(f"Processing sequences: {fasta_file}")
        self._log(f"{'='*60}")

        # Read sequences
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        self._log(f"Total sequences: {len(sequences)}")
        self._log(f"Active site positions: {len(self.active_site_positions)}")

        # Process sequences
        self._log("\nAligning sequences and extracting active sites...")
        active_site_sequences = []

        for i, record in enumerate(sequences):
            if (i + 1) % progress_interval == 0:
                self._log(f"  Progress: {i + 1}/{len(sequences)}")

            # Align to reference
            aligned_ref, aligned_target, score = self.align_sequences(
                self.ref_sequence,
                str(record.seq)
            )

            if aligned_ref is not None:
                # Extract active site residues
                active_site = self.extract_sites_from_alignment(
                    aligned_ref,
                    aligned_target,
                    self.ref_residue_numbers,
                    self.active_site_positions
                )

                if len(active_site) == len(self.active_site_positions):
                    # Create new record with active site sequence
                    new_record = SeqRecord(
                        Seq(active_site),
                        id=record.id,
                        description=f"Active site residues ({len(active_site)} positions)"
                    )
                    active_site_sequences.append(new_record)

        self._log(f"\nSuccessfully processed: {len(active_site_sequences)}/{len(sequences)}")

        # Save to file if requested
        if output_file:
            SeqIO.write(active_site_sequences, output_file, "fasta")
            self._log(f"Active site sequences saved: {output_file}")

        return active_site_sequences

    def batch_process_structures(self,
                                 structure_dir: str,
                                 output_dir: str,
                                 pattern: str = "*.pdb") -> List[Dict]:
        """
        Process multiple PDB structures in batch.

        Args:
            structure_dir: Directory containing PDB files
            output_dir: Output directory for results
            pattern: File pattern to match (default: *.pdb)

        Returns:
            List of result dictionaries
        """
        structure_path = Path(structure_dir)
        pdb_files = sorted(structure_path.glob(pattern))

        if not pdb_files:
            self._log(f"No PDB files found in {structure_dir} matching {pattern}")
            return []

        self._log(f"Found {len(pdb_files)} PDB files to process")

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        results = []
        for pdb_file in pdb_files:
            try:
                target_name = pdb_file.stem
                output_file = os.path.join(output_dir, f"{target_name}_active_sites.txt")

                result = self.extract_from_structure(
                    str(pdb_file),
                    output_file=output_file
                )
                results.append(result)

            except Exception as e:
                self._log(f"Error processing {pdb_file}: {e}")
                continue

        # Create summary
        if results:
            summary_file = os.path.join(output_dir, "active_sites_summary.tsv")
            self._write_batch_summary(results, summary_file)

        return results

    def _write_structure_results(self, results: Dict, output_file: str, target_pdb: str):
        """Write structure extraction results to file"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f"# Active sites found in {results['target']}\n")
            f.write(f"# Reference: {os.path.basename(self.reference_pdb)}\n")
            f.write(f"# Alignment RMSD: {results['rmsd']:.2f} Å\n")
            f.write(f"# Distance cutoff: {self.distance_cutoff} Å\n")
            f.write(f"# Chain: {results['chain']}\n")
            f.write(f"#\n")
            f.write(f"# Format: PDB_file[TAB]Chain[TAB]Residue_numbers\n")
            f.write(f"#\n")

            residue_numbers = ','.join(map(str, results['residue_numbers']))
            f.write(f"{os.path.basename(target_pdb)}\t{results['chain']}\t{residue_numbers}\n")

            f.write(f"\n# Detailed residue information:\n")
            f.write(f"# ResNum\tResName\tDistance(Å)\n")
            for site in results['active_sites']:
                f.write(f"# {site[0]}\t{site[1]}\t{site[2]:.2f}\n")

        self._log(f"Results saved to: {output_file}")

    def _write_batch_summary(self, results: List[Dict], summary_file: str):
        """Write batch processing summary"""
        with open(summary_file, 'w', encoding='utf-8') as f:
            # Header
            f.write("Structure\tChain\tRMSD(Å)\tN_Sites\tResidue_Positions\tResidue_Types\n")

            # Reference
            ref_name = Path(self.reference_pdb).stem
            ref_res_str = ','.join(map(str, self.active_site_positions))
            f.write(f"{ref_name}\t{self.reference_chain}\t0.00\t"
                   f"{len(self.active_site_positions)}\t{ref_res_str}\t-\n")

            # Targets
            for result in results:
                res_positions = ','.join(map(str, result['residue_numbers']))
                res_types = ','.join(result['residue_types'])

                f.write(f"{result['target']}\t{result['chain']}\t{result['rmsd']:.2f}\t"
                       f"{result['n_sites']}\t{res_positions}\t{res_types}\n")

        self._log(f"Summary saved to: {summary_file}")


def main():
    """Example usage of ActiveSiteExtractor"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract active site residues from protein structures and sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--reference-pdb', '-r', required=True,
                       help='Reference PDB structure')
    parser.add_argument('--reference-sites', '-s', required=True,
                       help='Active site definition (file or comma-separated numbers)')
    parser.add_argument('--reference-chain', '-c', default='A',
                       help='Reference chain ID (default: A)')
    parser.add_argument('--target-pdb', '-t',
                       help='Target PDB structure (single file)')
    parser.add_argument('--target-fasta', '-f',
                       help='Target sequences (FASTA file)')
    parser.add_argument('--target-dir', '-d',
                       help='Directory with multiple PDB files')
    parser.add_argument('--output', '-o', default='active_sites_output',
                       help='Output file or directory')
    parser.add_argument('--distance-cutoff', type=float, default=4.0,
                       help='Distance cutoff for structure matching (Å)')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress messages')

    args = parser.parse_args()

    # Parse sites (file or list)
    if os.path.exists(args.reference_sites):
        sites = args.reference_sites
    else:
        sites = [int(x.strip()) for x in args.reference_sites.split(',')]

    # Initialize extractor
    extractor = ActiveSiteExtractor(
        reference_pdb=args.reference_pdb,
        reference_sites=sites,
        reference_chain=args.reference_chain,
        distance_cutoff=args.distance_cutoff,
        verbose=not args.quiet
    )

    # Process based on input type
    if args.target_pdb:
        # Single PDB structure
        extractor.extract_from_structure(
            args.target_pdb,
            output_file=args.output
        )
    elif args.target_fasta:
        # FASTA sequences
        extractor.extract_from_sequences(
            args.target_fasta,
            output_file=args.output
        )
    elif args.target_dir:
        # Batch PDB structures
        extractor.batch_process_structures(
            args.target_dir,
            output_dir=args.output
        )
    else:
        print("Error: Specify --target-pdb, --target-fasta, or --target-dir")
        parser.print_help()
        return 1

    print("\nProcessing complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
