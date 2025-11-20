#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
UDH 서열에서 Active Site 잔기 추출
참조 구조와 정렬 후 active site 위치의 잔기들만 추출
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import re

def read_pdb_sequence(pdb_file):
    """PDB 파일에서 서열 읽기"""
    sequence = []
    residues = {}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()

                if chain == 'A' and resnum not in residues:
                    # 3-letter to 1-letter amino acid code
                    aa_map = {
                        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                    }

                    if resname in aa_map:
                        residues[resnum] = aa_map[resname]

    # Sort by residue number
    sorted_residues = sorted(residues.items())
    sequence = ''.join([aa for _, aa in sorted_residues])
    residue_numbers = [num for num, _ in sorted_residues]

    return sequence, residue_numbers

def align_to_reference(ref_seq, target_seq):
    """참조 서열과 타겟 서열을 정렬"""
    # Global alignment
    alignments = pairwise2.align.globalxx(ref_seq, target_seq)

    if alignments:
        best_alignment = alignments[0]
        return best_alignment.seqA, best_alignment.seqB, best_alignment.score
    else:
        return None, None, 0

def extract_active_site_residues(aligned_ref, aligned_target, ref_residue_numbers, active_site_positions):
    """정렬된 서열에서 active site 위치의 잔기 추출"""
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
        if resnum in active_site_positions:
            target_aa = aligned_target[aln_idx]
            active_site_residues.append(target_aa if target_aa != '-' else 'X')

    return ''.join(active_site_residues)

def main():
    print("=" * 70)
    print("Extracting Active Site Residues from UDH Sequences")
    print("=" * 70)
    print()

    # Input files
    pdb_file = Path("test_data/AtUdh_pdb3rfv_chainA.pdb")
    fasta_file = Path("test_data/UDHs_filtered_std2.5.fasta")

    # Active site positions (from pocket file)
    active_site_positions = [137, 138, 139, 140, 141, 143, 165, 166, 167,
                             189, 190, 191, 213, 214, 215, 237, 238, 257, 258, 259]

    print(f"Reading reference structure: {pdb_file}")
    ref_seq, ref_residue_numbers = read_pdb_sequence(pdb_file)
    print(f"  Reference sequence length: {len(ref_seq)}")
    print(f"  Active site positions: {len(active_site_positions)}")
    print()

    # Read target sequences
    print(f"Reading target sequences: {fasta_file}")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"  Total sequences: {len(sequences)}")
    print()

    # Process sequences
    print("Aligning sequences and extracting active sites...")
    active_site_sequences = []
    successful = 0

    for i, record in enumerate(sequences):
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i + 1}/{len(sequences)}")

        # Align to reference
        aligned_ref, aligned_target, score = align_to_reference(ref_seq, str(record.seq))

        if aligned_ref is not None:
            # Extract active site residues
            active_site = extract_active_site_residues(
                aligned_ref, aligned_target, ref_residue_numbers, active_site_positions
            )

            if len(active_site) == len(active_site_positions):
                # Create new record with active site sequence
                new_record = SeqRecord(
                    Seq(active_site),
                    id=record.id,
                    description=f"Active site residues ({len(active_site)} positions)"
                )
                active_site_sequences.append(new_record)
                successful += 1

    print(f"\n  Successfully processed: {successful}/{len(sequences)}")
    print()

    # Save active site sequences
    output_file = Path("udh_active_sites.fasta")
    SeqIO.write(active_site_sequences, output_file, "fasta")
    print(f"Active site sequences saved: {output_file}")
    print(f"  Total sequences: {len(active_site_sequences)}")
    print(f"  Active site length: {len(active_site_positions)} residues")
    print()
    print("=" * 70)
    print("Extraction Complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
