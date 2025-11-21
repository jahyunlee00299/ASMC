#!/usr/bin/env python3
"""
Download multiple UDH structures from PDB

This script downloads known UDH/UGDH structures from the PDB database.
"""

import os
import sys
from urllib import request
import time

# Known UDH structures in PDB
UDH_STRUCTURES = {
    '3RFV': 'Arabidopsis thaliana UDP-glucose dehydrogenase',
    '2Q3E': 'Human UDP-glucose dehydrogenase',
    '3PTZ': 'E. coli UDP-glucose dehydrogenase',
    '2B69': 'Streptococcus pyogenes UGDH',
    '3TF5': 'Klebsiella pneumoniae UGDH',
    '4RKQ': 'Burkholderia cepacia UGDH',
    '3UXH': 'Pseudomonas aeruginosa UGDH',
    '4YK5': 'Bacillus anthracis UGDH',
    '3MZ7': 'Candida albicans UGDH',
    '5U3M': 'Acinetobacter baumannii UGDH',
    '6B3H': 'Cryptococcus neoformans UGDH',
}

def download_pdb(pdb_id, output_dir='.'):
    """Download PDB file using urllib"""
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    output_file = os.path.join(output_dir, f'{pdb_id}.pdb')

    try:
        print(f"Downloading {pdb_id}...", end=' ', flush=True)
        with request.urlopen(url, timeout=30) as response:
            content = response.read()

            if len(content) < 100:  # Too small to be a valid PDB
                print(f"FAILED (invalid content)")
                return False

            with open(output_file, 'wb') as f:
                f.write(content)

            file_size = len(content) / 1024  # KB
            print(f"OK ({file_size:.1f} KB)")
            return True

    except Exception as e:
        print(f"FAILED ({str(e)[:50]})")
        return False

def main():
    output_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(output_dir, exist_ok=True)

    print(f"Downloading {len(UDH_STRUCTURES)} UDH structures to {output_dir}/\n")

    successful = []
    failed = []

    for i, (pdb_id, description) in enumerate(UDH_STRUCTURES.items(), 1):
        print(f"[{i}/{len(UDH_STRUCTURES)}] {pdb_id}: {description}")

        if download_pdb(pdb_id, output_dir):
            successful.append(pdb_id)
        else:
            failed.append(pdb_id)

        # Be nice to the server
        if i < len(UDH_STRUCTURES):
            time.sleep(0.5)

    print(f"\n{'='*60}")
    print(f"Download complete!")
    print(f"Successful: {len(successful)}/{len(UDH_STRUCTURES)}")
    if failed:
        print(f"Failed: {', '.join(failed)}")
    print(f"{'='*60}")

    return len(successful)

if __name__ == "__main__":
    sys.exit(0 if main() > 0 else 1)
