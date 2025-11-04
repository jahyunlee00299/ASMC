#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ASMC 실제 실행 스크립트
작성일: 2025년 11월 4일

이 스크립트는 ASMC를 실제로 실행하는 예제입니다.
테스트 데이터를 자동으로 생성하고 실행할 수 있습니다.
"""

import subprocess
import sys
import os
from pathlib import Path
import time

def ensure_test_data():
    """테스트 데이터가 없으면 생성"""
    if not Path("test_data").exists():
        print("📁 테스트 데이터 디렉토리 생성 중...")
        create_full_test_files()
    return True

def create_full_test_files():
    """실제 실행 가능한 테스트 파일들 생성"""
    print("\n" + "=" * 60)
    print("테스트 데이터 생성 중...")
    print("=" * 60)

    # test_data 디렉토리 생성
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)

    # 1. 간단한 PDB 파일 생성 (최소한의 구조)
    pdb_content = """HEADER    TEST PROTEIN                            01-NOV-24   TEST
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
ATOM     11  CG1 VAL A   2       4.100   5.050   1.200  1.00 20.00           C
ATOM     12  CG2 VAL A   2       3.850   2.950   2.500  1.00 20.00           C
TER      13      VAL A   2
END
"""

    # 테스트 PDB 파일들 생성
    for i in range(1, 3):
        pdb_file = test_dir / f"test_protein{i}.pdb"
        pdb_file.write_text(pdb_content, encoding='utf-8')
        print(f"✅ {pdb_file} 생성됨")

    # 2. references.txt 생성 (실제 경로 사용)
    refs_content = f"{test_dir.absolute()}/test_protein1.pdb\n{test_dir.absolute()}/test_protein2.pdb"
    refs_file = test_dir / "references.txt"
    refs_file.write_text(refs_content, encoding='utf-8')
    print(f"✅ {refs_file} 생성됨")

    # 3. sequences.fasta 생성
    fasta_content = """>Test_Protein_1
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Test_Protein_2
MKHLWFFLLLVAAPRWVLSAAGACGQEARPEAVGQHWEALGRFWDYLRWVQTLSEQVQEE
LLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Test_Protein_3
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEELLSSQVTQELRALM
"""
    fasta_file = test_dir / "sequences.fasta"
    fasta_file.write_text(fasta_content, encoding='utf-8')
    print(f"✅ {fasta_file} 생성됨")

    # 4. models.txt 생성
    models_content = f"{test_dir.absolute()}/test_protein1.pdb\ttest_protein1\n{test_dir.absolute()}/test_protein2.pdb\ttest_protein1"
    models_file = test_dir / "models.txt"
    models_file.write_text(models_content, encoding='utf-8')
    print(f"✅ {models_file} 생성됨")

    # 5. pocket.txt 생성 (선택사항)
    pocket_content = f"{test_dir.absolute()}/test_protein1.pdb\tA\t1,2\n{test_dir.absolute()}/test_protein2.pdb\tA\t1,2"
    pocket_file = test_dir / "pocket.txt"
    pocket_file.write_text(pocket_content, encoding='utf-8')
    print(f"✅ {pocket_file} 생성됨")

    print("\n📌 테스트 데이터 생성 완료!")
    print(f"   위치: {test_dir.absolute()}")
    return test_dir

def run_command(cmd, description="명령어 실행 중"):
    """명령어를 실제로 실행하고 결과 표시"""
    print(f"\n🔧 {description}...")
    print(f"명령어: {' '.join(cmd)}")
    print("-" * 40)

    try:
        # 실시간 출력을 위해 Popen 사용
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )

        # 실시간으로 출력 표시
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())

        # 남은 에러 메시지 확인
        stderr = process.stderr.read()
        if stderr:
            print(f"⚠️ 경고/에러:\n{stderr}")

        rc = process.poll()
        if rc == 0:
            print("✅ 실행 완료!")
        else:
            print(f"❌ 실행 중 오류 발생 (종료 코드: {rc})")

        return rc == 0

    except FileNotFoundError:
        print("❌ Python을 찾을 수 없습니다.")
        return False
    except Exception as e:
        print(f"❌ 실행 중 오류: {e}")
        return False

def run_asmc_help():
    """ASMC 도움말 표시"""
    cmd = [sys.executable, "-m", "asmc.run_asmc", "--help"]
    run_command(cmd, "ASMC 도움말 표시")

def run_asmc_test():
    """간단한 테스트 실행 (포켓 검출까지만)"""
    print("\n" + "=" * 60)
    print("간단한 테스트 실행 (포켓 검출)")
    print("=" * 60)

    ensure_test_data()
    test_dir = Path("test_data")
    output_dir = Path("output_test")

    cmd = [
        sys.executable,
        "-m", "asmc.run_asmc",
        "run",
        "-r", str(test_dir / "references.txt"),
        "-p", str(test_dir / "pocket.txt"),
        "-s", str(test_dir / "sequences.fasta"),
        "-o", str(output_dir),
        "-t", "2",
        "--end", "pocket",  # 포켓 검출까지만
        "--id", "20"        # 낮은 임계값
    ]

    return run_command(cmd, "포켓 검출 테스트")

def run_identity_test():
    """서열 유사도 계산 테스트"""
    print("\n" + "=" * 60)
    print("서열 유사도 계산 테스트")
    print("=" * 60)

    ensure_test_data()
    test_dir = Path("test_data")

    cmd = [
        sys.executable,
        "-m", "asmc.run_asmc",
        "identity",
        "-s", str(test_dir / "sequences.fasta"),
        "-r", str(test_dir / "references.txt"),
        "-o", "identity_results.txt"
    ]

    success = run_command(cmd, "서열 유사도 계산")

    if success and Path("identity_results.txt").exists():
        print("\n📄 결과 파일 내용:")
        print("-" * 40)
        content = Path("identity_results.txt").read_text(encoding='utf-8')
        print(content[:500])  # 처음 500자만 표시
        if len(content) > 500:
            print("... (생략) ...")

    return success

def run_clustering_test():
    """클러스터링 테스트 (모델 기반)"""
    print("\n" + "=" * 60)
    print("클러스터링 테스트")
    print("=" * 60)

    ensure_test_data()
    test_dir = Path("test_data")
    output_dir = Path("output_clustering")

    cmd = [
        sys.executable,
        "-m", "asmc.run_asmc",
        "run",
        "-m", str(test_dir / "models.txt"),
        "-r", str(test_dir / "references.txt"),
        "-o", str(output_dir),
        "--end", "clustering",
        "-e", "0.5",
        "--min-samples", "1"  # 작은 테스트이므로 1로 설정
    ]

    return run_command(cmd, "클러스터링 실행")

def run_full_pipeline():
    """전체 파이프라인 실행 (시간이 오래 걸릴 수 있음)"""
    print("\n" + "=" * 60)
    print("전체 파이프라인 실행")
    print("⚠️ 주의: 이 작업은 시간이 오래 걸릴 수 있습니다!")
    print("=" * 60)

    response = input("\n정말 실행하시겠습니까? (y/n): ").strip().lower()
    if response != 'y':
        print("취소되었습니다.")
        return False

    ensure_test_data()
    test_dir = Path("test_data")
    output_dir = Path("output_full")

    cmd = [
        sys.executable,
        "-m", "asmc.run_asmc",
        "run",
        "-s", str(test_dir / "sequences.fasta"),
        "-r", str(test_dir / "references.txt"),
        "-o", str(output_dir),
        "-t", "4",
        "--id", "20"
    ]

    return run_command(cmd, "전체 파이프라인 실행")

def check_outputs():
    """생성된 출력 파일들 확인"""
    print("\n" + "=" * 60)
    print("생성된 출력 파일 확인")
    print("=" * 60)

    output_dirs = ["output_test", "output_clustering", "output_full"]

    for dir_name in output_dirs:
        dir_path = Path(dir_name)
        if dir_path.exists():
            print(f"\n📁 {dir_name}/")
            for file_path in dir_path.rglob("*"):
                if file_path.is_file():
                    size = file_path.stat().st_size
                    print(f"  - {file_path.relative_to(dir_path)} ({size:,} bytes)")
        else:
            print(f"\n❌ {dir_name}/ 디렉토리가 없습니다.")

    if Path("identity_results.txt").exists():
        print(f"\n📄 identity_results.txt ({Path('identity_results.txt').stat().st_size:,} bytes)")

def clean_outputs():
    """출력 파일들 정리"""
    print("\n정리할 항목:")
    print("1. 출력 디렉토리들 (output_*)")
    print("2. 결과 파일 (identity_results.txt)")
    print("3. 테스트 데이터 (test_data/)")
    print("4. 모두 정리")

    choice = input("\n선택 (1-4): ").strip()

    if choice in ["1", "4"]:
        for dir_name in ["output_test", "output_clustering", "output_full"]:
            if Path(dir_name).exists():
                import shutil
                shutil.rmtree(dir_name)
                print(f"✅ {dir_name}/ 삭제됨")

    if choice in ["2", "4"]:
        if Path("identity_results.txt").exists():
            Path("identity_results.txt").unlink()
            print("✅ identity_results.txt 삭제됨")

    if choice in ["3", "4"]:
        if Path("test_data").exists():
            import shutil
            shutil.rmtree("test_data")
            print("✅ test_data/ 삭제됨")

def main():
    """메인 실행 함수"""
    while True:
        print("\n" + "=" * 60)
        print("ASMC 실제 실행 스크립트")
        print("=" * 60)

        print("\n실행 옵션:")
        print("1. 📋 ASMC 도움말 보기")
        print("2. 🧪 간단한 테스트 (포켓 검출)")
        print("3. 🔍 서열 유사도 계산")
        print("4. 📊 클러스터링 테스트")
        print("5. 🚀 전체 파이프라인 실행 (주의: 오래 걸림)")
        print("6. 📁 출력 파일 확인")
        print("7. 🗑️ 출력 파일 정리")
        print("8. 🔄 테스트 데이터 재생성")
        print("0. 종료")
        print()

        try:
            choice = input("선택하세요 (0-8): ").strip()

            if choice == "1":
                run_asmc_help()
            elif choice == "2":
                run_asmc_test()
            elif choice == "3":
                run_identity_test()
            elif choice == "4":
                run_clustering_test()
            elif choice == "5":
                run_full_pipeline()
            elif choice == "6":
                check_outputs()
            elif choice == "7":
                clean_outputs()
            elif choice == "8":
                create_full_test_files()
            elif choice == "0":
                print("\n프로그램을 종료합니다.")
                break
            else:
                print("❌ 잘못된 선택입니다.")

            input("\n계속하려면 Enter를 누르세요...")

        except KeyboardInterrupt:
            print("\n\n프로그램이 중단되었습니다.")
            break
        except Exception as e:
            print(f"\n❌ 에러 발생: {e}")
            input("\n계속하려면 Enter를 누르세요...")

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("💡 ASMC 실제 실행 스크립트에 오신 것을 환영합니다!")
    print("=" * 60)
    print("\n이 스크립트는 ASMC를 실제로 실행합니다.")
    print("테스트 데이터가 없으면 자동으로 생성됩니다.")

    # Python 버전 확인
    print(f"\n🐍 Python 버전: {sys.version}")

    # ASMC 설치 확인
    try:
        import asmc
        print("✅ ASMC 모듈이 설치되어 있습니다.")
    except ImportError:
        print("❌ ASMC 모듈을 찾을 수 없습니다.")
        print("   'pip install -e .' 명령으로 설치해주세요.")
        sys.exit(1)

    main()

    print("\n" + "=" * 60)
    print("프로그램이 종료되었습니다. 감사합니다!")
    print("=" * 60)