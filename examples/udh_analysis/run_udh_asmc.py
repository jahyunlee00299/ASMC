#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
UDH (Uronate Dehydrogenase) ASMC 분석
AtUdh PDB 구조 기반 Active Site Motif Clustering
"""

import subprocess
import sys
from pathlib import Path
from datetime import datetime

def run_asmc_analysis():
    """ASMC를 사용한 UDH 클러스터링 분석"""

    print("=" * 70)
    print("UDH ASMC (Active Site Motif Clustering) Analysis")
    print("=" * 70)
    print()

    # 입력 파일 경로
    pdb_file = Path("test_data/AtUdh_pdb3rfv_chainA.pdb")
    fasta_file = Path("test_data/UDHs_filtered_std2.5.fasta")

    # 파일 확인
    if not pdb_file.exists():
        print(f"Error: PDB file not found: {pdb_file}")
        return 1
    if not fasta_file.exists():
        print(f"Error: FASTA file not found: {fasta_file}")
        return 1

    print(f"PDB structure: {pdb_file}")
    print(f"Sequences: {fasta_file}")
    print()

    # UDH active site 정의
    # PDB 3RFV는 UDP-glucose와 NAD+ 결합 부위를 가지고 있음
    # 문헌 기반으로 촉매 잔기와 기질 결합 잔기 선택
    # Typical dehydrogenase active site includes:
    # - NAD+ binding residues
    # - Catalytic residues (usually Ser, Lys, His, Tyr)
    # - Substrate binding pocket residues

    # 3RFV 구조에서 중요한 active site 잔기들 (예시)
    # 실제로는 ligand 주변 5-6Å 이내 잔기를 찾거나 문헌 참고
    active_site_residues = [
        "143",  # Catalytic residue (예: Ser)
        "144",
        "145",
        "165",  # NAD binding
        "166",
        "167",
        "189",  # Substrate binding
        "190",
        "191",
        "213",
        "214",
        "215",
        "257",  # Active site pocket
        "258",
        "259",
    ]

    # Active site를 pocket 파일로 저장
    pocket_file = Path("udh_active_site.txt")
    with open(pocket_file, 'w') as f:
        for res in active_site_residues:
            f.write(f"A {res}\n")  # Chain A, residue number

    print(f"Active site definition saved: {pocket_file}")
    print(f"Number of active site residues: {len(active_site_residues)}")
    print()

    # 출력 디렉토리 생성
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = Path(f"udh_asmc_output_{timestamp}")
    output_dir.mkdir(exist_ok=True)

    print(f"Output directory: {output_dir}")
    print()

    # ASMC 명령어 구성
    # asmc 명령어 형식:
    # python -m asmc.run_asmc -p PDB -s FASTA -o OUTPUT [-r RESIDUES] [options]

    cmd = [
        sys.executable,
        "-m", "asmc.run_asmc",
        "-p", str(pdb_file),
        "-s", str(fasta_file),
        "-o", str(output_dir),
        "-r", str(pocket_file),  # active site residues
        "-c", "dbscan",  # clustering method
        "-e", "0.5",  # epsilon for DBSCAN
        "-m", "3",  # min samples for DBSCAN
        "--logo",  # generate sequence logos
    ]

    print("Running ASMC...")
    print("Command:", " ".join(cmd))
    print("=" * 70)
    print()

    # ASMC 실행
    try:
        result = subprocess.run(
            cmd,
            capture_output=False,  # 실시간 출력
            text=True,
            check=True
        )

        print()
        print("=" * 70)
        print("ASMC Analysis Complete!")
        print("=" * 70)
        print(f"Results saved in: {output_dir}")

        # 생성된 파일 확인
        output_files = list(output_dir.glob("*"))
        if output_files:
            print("\nGenerated files:")
            for f in sorted(output_files):
                print(f"  - {f.name}")

        return 0

    except subprocess.CalledProcessError as e:
        print()
        print("=" * 70)
        print("Error: ASMC execution failed")
        print("=" * 70)
        print(f"Return code: {e.returncode}")
        return e.returncode

    except Exception as e:
        print()
        print("=" * 70)
        print("Error occurred during ASMC execution")
        print("=" * 70)
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    exit_code = run_asmc_analysis()
    sys.exit(exit_code)
