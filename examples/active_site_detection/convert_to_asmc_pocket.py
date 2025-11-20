#!/usr/bin/env python3
"""
Convert active site comparison results to ASMC pocket format

This script converts the TSV output from find_udh_active_sites.py
to the pocket format used by ASMC pipeline.

Usage:
    python convert_to_asmc_pocket.py results/active_sites_comparison.tsv udh_pocket.txt

Input format (TSV):
    Structure	Chain	RMSD(Å)	N_Sites	Residue_Positions	Residue_Types
    udh1	A	1.23	22	10,12,34,36,...	THR,ASP,SER,...

Output format (ASMC pocket):
    /path/to/udh1.pdb	A	10,12,34,36,...

Author: ASMC
Date: 2025-11-20
"""

import sys
import os
import argparse
from pathlib import Path


def convert_to_pocket_format(input_tsv: str, output_file: str,
                            pdb_dir: str = "pdb_structures",
                            skip_reference: bool = True):
    """
    Convert TSV comparison file to ASMC pocket format

    Args:
        input_tsv: Input TSV file from find_udh_active_sites.py
        output_file: Output pocket file path
        pdb_dir: Directory containing PDB files
        skip_reference: Skip the reference structure (first line)
    """
    if not os.path.exists(input_tsv):
        print(f"Error: Input file not found: {input_tsv}")
        sys.exit(1)

    # Read TSV file
    with open(input_tsv, 'r') as f:
        lines = f.readlines()

    if len(lines) < 2:
        print("Error: Input file is empty or has no data")
        sys.exit(1)

    # Parse header
    header = lines[0].strip().split('\t')
    print(f"Header: {header}")

    # Find column indices
    try:
        struct_idx = header.index('Structure')
        chain_idx = header.index('Chain')
        residues_idx = header.index('Residue_Positions')
    except ValueError as e:
        print(f"Error: Required column not found in header: {e}")
        sys.exit(1)

    # Write output
    with open(output_file, 'w') as out:
        out.write("# ASMC Pocket Definition File\n")
        out.write("# Generated from UDH active site finder results\n")
        out.write("# Format: PDB_path[TAB]Chain[TAB]Residue_numbers\n")
        out.write("#\n")

        # Process data lines
        for i, line in enumerate(lines[1:], start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) <= max(struct_idx, chain_idx, residues_idx):
                print(f"Warning: Skipping malformed line {i+1}")
                continue

            structure = parts[struct_idx]
            chain = parts[chain_idx]
            residues = parts[residues_idx]

            # Skip reference if requested
            if skip_reference and i == 1:
                out.write(f"# Reference (skipped): {structure}\n")
                continue

            # Construct PDB path
            # Assume structure name is the PDB filename without extension
            pdb_filename = f"{structure}.pdb"
            pdb_path = os.path.join(pdb_dir, pdb_filename)

            # Write in ASMC pocket format
            out.write(f"{pdb_path}\t{chain}\t{residues}\n")

    print(f"\nConversion complete!")
    print(f"Input:  {input_tsv}")
    print(f"Output: {output_file}")
    print(f"Converted {len(lines)-1-int(skip_reference)} structures")


def main():
    parser = argparse.ArgumentParser(
        description="Convert active site comparison TSV to ASMC pocket format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion
  python convert_to_asmc_pocket.py results/active_sites_comparison.tsv udh_pocket.txt

  # Specify PDB directory
  python convert_to_asmc_pocket.py results/active_sites_comparison.tsv udh_pocket.txt \\
                                    --pdb-dir /path/to/pdb/files

  # Include reference structure
  python convert_to_asmc_pocket.py results/active_sites_comparison.tsv udh_pocket.txt \\
                                    --include-reference
        """
    )

    parser.add_argument('input', help='Input TSV file from find_udh_active_sites.py')
    parser.add_argument('output', help='Output pocket file for ASMC')
    parser.add_argument('--pdb-dir', default='pdb_structures',
                       help='Directory containing PDB files (default: pdb_structures)')
    parser.add_argument('--include-reference', action='store_true',
                       help='Include reference structure in output (default: skip)')

    args = parser.parse_args()

    convert_to_pocket_format(
        input_tsv=args.input,
        output_file=args.output,
        pdb_dir=args.pdb_dir,
        skip_reference=not args.include_reference
    )


if __name__ == "__main__":
    main()
