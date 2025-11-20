#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run ASMC with AlphaFold structures (no MODELLER required)

This example shows how to use ASMC with pre-existing AlphaFold structures,
completely bypassing the MODELLER alignment and modeling steps.

Usage:
    python run_asmc_with_alphafold.py

Author: ASMC
Date: 2025-11-20
"""

import subprocess
import sys
import os
from pathlib import Path

# Fix Windows encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')


def run_asmc_with_alphafold_structures():
    """
    Run ASMC using AlphaFold structures from the pdb_structures directory
    """

    print("="*70)
    print("🧬 ASMC with AlphaFold Structures (MODELLER-free)")
    print("="*70)
    print()

    # Paths
    base_dir = Path(__file__).parent
    alphafold_dir = base_dir / "pdb_structures" / "Alphafold_server"
    reference_pdb = base_dir / "pdb_structures" / "pdb3rfv_chainA.pdb"

    # Check if directories exist
    if not alphafold_dir.exists():
        print(f"❌ Error: AlphaFold directory not found: {alphafold_dir}")
        print(f"   Please ensure your AlphaFold structures are in this directory")
        return 1

    if not reference_pdb.exists():
        print(f"❌ Error: Reference PDB not found: {reference_pdb}")
        return 1

    print(f"📁 AlphaFold structures: {alphafold_dir}")
    print(f"📁 Reference structure: {reference_pdb}")
    print()

    # Step 1: Create models.txt
    print("Step 1: Creating models.txt file...")
    print("-" * 70)

    models_txt = base_dir / "alphafold_models.txt"

    cmd_create = [
        sys.executable,
        str(base_dir / "create_models_txt_for_alphafold.py"),
        "--alphafold-dir", str(alphafold_dir),
        "--reference", str(reference_pdb),
        "--output", str(models_txt),
        "--pattern", "*.pdb"
    ]

    result = subprocess.run(cmd_create, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"❌ Error creating models.txt:")
        print(result.stderr)
        return 1

    print("✓ models.txt created successfully")
    print()

    # Step 2: Create reference and pocket files
    print("Step 2: Setting up reference and pocket files...")
    print("-" * 70)

    references_txt = base_dir / "references_alphafold.txt"
    pocket_txt = base_dir / "pocket_alphafold.txt"

    # Create references.txt
    with open(references_txt, 'w') as f:
        f.write(f"{reference_pdb.absolute()}\n")
    print(f"✓ Created {references_txt}")

    # Create pocket.txt (using active site info from atudh_active_site.txt)
    active_site_file = base_dir / "pdb_structures" / "atudh_active_site.txt"
    if active_site_file.exists():
        with open(active_site_file, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.startswith('#') and line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        chain = parts[1]
                        residues = parts[2]

                        # Create pocket.txt in ASMC format
                        # Format: pdb_name,chain,residue_numbers
                        with open(pocket_txt, 'w') as pf:
                            pf.write(f"{reference_pdb.stem},{chain},{residues}\n")

                        print(f"✓ Created {pocket_txt}")
                        print(f"  Chain: {chain}")
                        print(f"  Active sites: {residues[:50]}..." if len(residues) > 50 else f"  Active sites: {residues}")
                        break
    else:
        print(f"⚠️  Warning: {active_site_file} not found")
        print(f"   Creating basic pocket file...")
        with open(pocket_txt, 'w') as f:
            f.write(f"{reference_pdb.stem},A,10,12,34,36,58\n")
        print(f"✓ Created basic {pocket_txt}")

    print()

    # Step 3: Run ASMC
    print("Step 3: Running ASMC (pocket detection only)...")
    print("-" * 70)
    print()

    output_dir = base_dir / "asmc_alphafold_output"

    cmd_asmc = [
        sys.executable, "-m", "asmc.run_asmc", "run",
        "-r", str(references_txt),
        "-p", str(pocket_txt),
        "-m", str(models_txt),
        "-o", str(output_dir),
        "--end", "pocket",
        "-t", "4"
    ]

    print("Running command:")
    print(" ".join(cmd_asmc))
    print()

    result = subprocess.run(cmd_asmc)

    if result.returncode == 0:
        print()
        print("="*70)
        print("✓ SUCCESS! ASMC completed")
        print("="*70)
        print(f"Output directory: {output_dir}")
        print()
        print("Generated files:")
        if output_dir.exists():
            for f in output_dir.iterdir():
                print(f"  - {f.name}")
    else:
        print()
        print("="*70)
        print("❌ ASMC failed")
        print("="*70)
        return result.returncode

    return 0


if __name__ == "__main__":
    try:
        exit_code = run_asmc_with_alphafold_structures()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
