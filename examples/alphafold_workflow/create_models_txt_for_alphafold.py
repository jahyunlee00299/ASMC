#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create models.txt file for AlphaFold structures

This script creates a models.txt file that ASMC can use directly
with pre-existing AlphaFold structures, bypassing the MODELLER alignment step.

Usage:
    python create_models_txt_for_alphafold.py \\
        --alphafold-dir path/to/alphafold/structures/ \\
        --reference path/to/reference.pdb \\
        --output models.txt

Author: ASMC
Date: 2025-11-20
"""

import argparse
import sys
from pathlib import Path

# Fix Windows encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

def create_models_txt(alphafold_dir, reference_pdb, output_file, pattern="*.pdb"):
    """
    Create models.txt file mapping AlphaFold structures to reference

    Args:
        alphafold_dir: Directory containing AlphaFold PDB files
        reference_pdb: Path to reference structure
        output_file: Output models.txt path
        pattern: File pattern to match (default: *.pdb)

    Returns:
        Number of structures found
    """
    alphafold_path = Path(alphafold_dir)
    reference_path = Path(reference_pdb).absolute()
    output_path = Path(output_file)

    if not alphafold_path.exists():
        print(f"Error: AlphaFold directory not found: {alphafold_dir}")
        return 0

    if not reference_path.exists():
        print(f"Error: Reference PDB not found: {reference_pdb}")
        return 0

    # Find all PDB files
    if alphafold_path.is_file():
        # Single file provided
        pdb_files = [alphafold_path]
    else:
        # Directory provided
        pdb_files = sorted(alphafold_path.glob(pattern))
        # Also check subdirectories
        pdb_files += sorted(alphafold_path.glob(f"**/{pattern}"))

    if not pdb_files:
        print(f"Error: No PDB files found matching pattern '{pattern}' in {alphafold_dir}")
        return 0

    # Remove duplicates
    pdb_files = list(set(pdb_files))
    pdb_files.sort()

    print(f"Found {len(pdb_files)} structure(s)")
    print(f"Reference: {reference_path}")
    print(f"Output: {output_path}")
    print()

    # Create models.txt
    with open(output_path, 'w') as f:
        for pdb_file in pdb_files:
            pdb_abs = pdb_file.absolute()
            # Format: model_path reference_path
            f.write(f"{pdb_abs} {reference_path}\n")
            print(f"  Added: {pdb_file.name}")

    print(f"\n✓ Created {output_path} with {len(pdb_files)} entries")
    return len(pdb_files)


def main():
    parser = argparse.ArgumentParser(
        description="Create models.txt for AlphaFold structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single directory of AlphaFold structures
  python create_models_txt_for_alphafold.py \\
      --alphafold-dir alphafold_structures/ \\
      --reference reference_3rfv.pdb \\
      --output models.txt

  # Specific file pattern
  python create_models_txt_for_alphafold.py \\
      --alphafold-dir alphafold_structures/ \\
      --reference reference_3rfv.pdb \\
      --pattern "*_relaxed_*.pdb" \\
      --output models.txt

Then run ASMC with:
  python -m asmc.run_asmc run \\
      -r references.txt \\
      -p pocket.txt \\
      -m models.txt \\
      -o output_dir \\
      --end pocket
        """
    )

    parser.add_argument('--alphafold-dir', '-d', required=True,
                       help='Directory containing AlphaFold PDB structures')
    parser.add_argument('--reference', '-r', required=True,
                       help='Reference PDB structure file')
    parser.add_argument('--output', '-o', default='models.txt',
                       help='Output models.txt file (default: models.txt)')
    parser.add_argument('--pattern', '-p', default='*.pdb',
                       help='File pattern to match (default: *.pdb)')

    args = parser.parse_args()

    count = create_models_txt(
        args.alphafold_dir,
        args.reference,
        args.output,
        args.pattern
    )

    if count == 0:
        sys.exit(1)

    print("\n" + "="*70)
    print("Next steps:")
    print("="*70)
    print(f"1. Run ASMC pocket detection:")
    print(f"   python -m asmc.run_asmc run \\")
    print(f"       -r <references.txt> \\")
    print(f"       -p <pocket.txt> \\")
    print(f"       -m {args.output} \\")
    print(f"       -o <output_dir> \\")
    print(f"       --end pocket")
    print()
    print("2. Or run full pipeline (pocket + clustering):")
    print(f"   python -m asmc.run_asmc run \\")
    print(f"       -r <references.txt> \\")
    print(f"       -p <pocket.txt> \\")
    print(f"       -m {args.output} \\")
    print(f"       -o <output_dir>")
    print("="*70)


if __name__ == "__main__":
    main()
