#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ASMC 실제 실행 및 결과 시각화"""

import subprocess
import sys
from pathlib import Path
import json

def create_test_data():
    """테스트 데이터 생성"""
    print("Creating test data...")
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)

    # Simple PDB file
    pdb_content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       1.950  -0.800   1.200  1.00 20.00           C
ATOM      6  N   VAL A   2       3.300   1.550   0.000  1.00 20.00           N
ATOM      7  CA  VAL A   2       3.950   2.850   0.000  1.00 20.00           C
ATOM      8  C   VAL A   2       5.450   2.750   0.000  1.00 20.00           C
ATOM      9  O   VAL A   2       6.100   1.700   0.000  1.00 20.00           O
ATOM     10  CB  VAL A   2       3.500   3.650   1.200  1.00 20.00           C
TER
END
"""

    # Create PDB files
    for i in range(1, 4):
        pdb_file = test_dir / f"protein{i}.pdb"
        pdb_file.write_text(pdb_content, encoding='utf-8')

    # Create references.txt
    refs_content = "\n".join([
        str((test_dir / f"protein{i}.pdb").absolute())
        for i in range(1, 3)
    ])
    (test_dir / "references.txt").write_text(refs_content, encoding='utf-8')

    # Create sequences.fasta
    fasta_content = """>Protein_1
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Protein_2
MKHLWFFLLLVAAPRWVLSAAGACGQEARPEAVGQHWEALGRFWDYLRWVQTLSEQVQEE
LLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Protein_3
MKTLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAG
"""
    (test_dir / "sequences.fasta").write_text(fasta_content, encoding='utf-8')

    # Create models.txt
    models_content = "\n".join([
        f"{(test_dir / f'protein{i}.pdb').absolute()}\tprotein1"
        for i in range(1, 4)
    ])
    (test_dir / "models.txt").write_text(models_content, encoding='utf-8')

    # Create pocket.txt
    pocket_content = "\n".join([
        f"{(test_dir / f'protein{i}.pdb').absolute()}\tA\t1,2"
        for i in range(1, 3)
    ])
    (test_dir / "pocket.txt").write_text(pocket_content, encoding='utf-8')

    print(f"Test data created in {test_dir.absolute()}")
    return test_dir

def run_identity_calculation():
    """서열 유사도 계산"""
    print("\n" + "="*60)
    print("Running sequence identity calculation...")
    print("="*60)

    test_dir = Path("test_data")
    cmd = [
        sys.executable, "-m", "asmc.run_asmc", "identity",
        "-s", str(test_dir / "sequences.fasta"),
        "-r", str(test_dir / "references.txt"),
        "-o", "identity_results.txt"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print("Identity calculation completed!")
        if Path("identity_results.txt").exists():
            content = Path("identity_results.txt").read_text()
            print("\nResults:")
            print(content)
            return content
    else:
        print(f"Error: {result.stderr}")
    return None

def run_clustering():
    """클러스터링 실행"""
    print("\n" + "="*60)
    print("Running clustering...")
    print("="*60)

    test_dir = Path("test_data")
    output_dir = Path("output_clustering")

    cmd = [
        sys.executable, "-m", "asmc.run_asmc", "run",
        "-m", str(test_dir / "models.txt"),
        "-r", str(test_dir / "references.txt"),
        "-p", str(test_dir / "pocket.txt"),
        "-o", str(output_dir),
        "--end", "clustering",
        "-e", "0.3",
        "--min-samples", "2"
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print("Clustering completed!")
        # Check for output files
        if output_dir.exists():
            print(f"\nOutput files in {output_dir}:")
            for file in output_dir.rglob("*"):
                if file.is_file():
                    print(f"  - {file}")
                    if file.suffix in ['.txt', '.tsv', '.csv']:
                        content = file.read_text()[:500]
                        print(f"    Content preview: {content[:100]}...")
    else:
        print(f"Error: {result.stderr}")
        print(f"Output: {result.stdout}")

    return output_dir

def main():
    print("="*60)
    print("ASMC Execution and Visualization")
    print("="*60)

    # Create test data
    test_dir = create_test_data()

    # Run identity calculation
    identity_results = run_identity_calculation()

    # Run clustering
    output_dir = run_clustering()

    print("\n" + "="*60)
    print("Execution completed!")
    print("="*60)

if __name__ == "__main__":
    main()