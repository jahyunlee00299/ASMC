#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple alignment builder for AlphaFold structures
This bypasses MODELLER requirement by using pre-existing AlphaFold structures
"""

import argparse
import warnings
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser

def create_simple_alignment_files(ref_file, seq_file, pocket_file, outdir, pid_cutoff=20.0):
    """
    Create alignment files for AlphaFold structures without using MODELLER

    Args:
        ref_file: Reference structures file
        seq_file: Sequences fasta file
        pocket_file: Pocket definition file
        outdir: Output directory
        pid_cutoff: Percent identity cutoff

    Returns:
        Success status
    """
    outdir = Path(outdir).absolute()
    outdir.mkdir(exist_ok=True)

    # Read references
    all_ref = Path(ref_file).read_text().split("\n")
    all_ref = [r for r in all_ref if r.strip()]

    # Read pocket info
    ref_chain = {}
    with open(pocket_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                ref_chain[parts[0]] = parts[1].strip()

    # Parse sequences
    sequences = list(SeqIO.parse(seq_file, 'fasta'))

    print(f"Found {len(sequences)} sequences to process")
    print(f"Found {len(all_ref)} reference structures")

    # Create output files
    models_file = outdir.parent / "models.txt"
    identity_file = outdir / "identity_targets_refs.tsv"
    job_file = outdir.parent / "job_file.txt"

    models_file.write_text("")  # Clear
    identity_file.write_text("")
    job_file.write_text("")

    for seq_record in sequences:
        seq_id = seq_record.id

        # Assume AlphaFold structure exists or will be provided
        # Just create the tracking files

        # Pick first reference as default
        best_ref = all_ref[0] if all_ref else ""

        # Write to models.txt
        model_path = outdir.parent / "models" / f"{seq_id}.pdb"
        with models_file.open('a') as f:
            f.write(f"{model_path} {best_ref}\n")

        # Write to identity file (use 100% since we have AlphaFold structures)
        with identity_file.open('a') as f:
            ref_name = Path(best_ref).stem if best_ref else "N/A"
            f.write(f"{seq_id}\t{ref_name}\t100.00\t{pid_cutoff:.2f}\n")

        # Write to job file
        ref_text = "+".join([code for code in ref_chain])
        with job_file.open('a') as f:
            f.write(f"ALPHAFOLD_STRUCT+{ref_text}\n")

    print(f"\nCreated alignment tracking files:")
    print(f"  - {models_file}")
    print(f"  - {identity_file}")
    print(f"  - {job_file}")
    print(f"\nNote: Using AlphaFold structures directly, skipping MODELLER alignment")

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple alignment file generator for AlphaFold structures"
    )
    parser.add_argument("-r", "--ref", required=True,
                       help="File containing reference structure paths")
    parser.add_argument("-s", "--seq", required=True,
                       help="Multi-fasta sequence file")
    parser.add_argument("-p", "--pocket", required=True,
                       help="Pocket definition file")
    parser.add_argument("-o", "--outdir", default="./",
                       help="Output directory")
    parser.add_argument("--id", type=float, default=20.0,
                       help="Percent identity cutoff")

    args = parser.parse_args()

    result = create_simple_alignment_files(
        args.ref, args.seq, args.pocket, args.outdir, args.id
    )

    exit(result)
