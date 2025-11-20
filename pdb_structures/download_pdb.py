#!/usr/bin/env python3
"""Download PDB structure using BioPython"""

from Bio.PDB import PDBList
import sys

def download_pdb(pdb_id, output_dir='.'):
    """Download PDB file using BioPython PDBList"""
    pdbl = PDBList()
    try:
        filename = pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format='pdb')
        print(f"Downloaded {pdb_id} to {filename}")
        return filename
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None

if __name__ == "__main__":
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "3RFV"
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    download_pdb(pdb_id, output_dir)
