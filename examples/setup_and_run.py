#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ASMC 자동 설치 및 실행 스크립트
PyCharm 또는 다른 Python 환경에서 직접 실행 가능
"""

import sys
import subprocess
import os
from pathlib import Path


def check_python_version():
    """Python 버전 확인"""
    version = sys.version_info
    print(f"🐍 Python 버전: {sys.version}")

    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print("❌ Python 3.8 이상이 필요합니다!")
        return False
    return True


def install_dependencies():
    """필요한 의존성 설치"""
    print("\n📦 의존성 패키지 설치 중...")

    # Python 버전에 따라 다른 requirements 파일 사용
    if sys.version_info.minor >= 14:
        requirements_file = "requirements.txt"
    else:
        requirements_file = "requirements_py38.txt"

    if not Path(requirements_file).exists():
        # requirements_py38.txt가 없으면 기본 패키지 설치
        packages = [
            "biopython>=1.81",
            "numpy",
            "scikit-learn",
            "scipy",
            "pyyaml",
            "pillow",
            "plotnineseqsuite",
            "pytest",
            "xlsxwriter"
        ]

        for package in packages:
            print(f"  설치 중: {package}")
            subprocess.run([sys.executable, "-m", "pip", "install", package],
                         capture_output=True, text=True)
    else:
        # requirements 파일 사용
        print(f"  {requirements_file} 사용")
        result = subprocess.run(
            [sys.executable, "-m", "pip", "install", "-r", requirements_file],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"⚠️ 일부 패키지 설치 실패: {result.stderr}")


def install_asmc():
    """ASMC 모듈 설치 (개발 모드)"""
    print("\n🔧 ASMC 설치 중...")

    # 현재 디렉토리에서 설치
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "-e", "."],
        capture_output=True, text=True
    )

    if result.returncode == 0:
        print("✅ ASMC 설치 완료!")
        return True
    else:
        print(f"❌ ASMC 설치 실패: {result.stderr}")
        return False


def verify_installation():
    """설치 확인"""
    try:
        import asmc
        print("\n✅ ASMC 모듈이 정상적으로 임포트됩니다.")
        return True
    except ImportError as e:
        print(f"\n❌ ASMC 모듈을 찾을 수 없습니다: {e}")
        return False


def run_example_script():
    """실행_예제.py 실행"""
    print("\n" + "=" * 60)
    print("🚀 실행_예제.py를 시작합니다...")
    print("=" * 60)

    example_script = Path("실행_예제.py")

    if not example_script.exists():
        print("❌ 실행_예제.py 파일을 찾을 수 없습니다!")
        return

    # 실행_예제.py 실행
    subprocess.run([sys.executable, str(example_script)])


def main():
    """메인 함수"""
    print("\n" + "=" * 60)
    print("ASMC 자동 설치 및 실행 스크립트")
    print("=" * 60)

    # 1. Python 버전 확인
    if not check_python_version():
        return

    # 2. ASMC 설치 확인
    if not verify_installation():
        print("\nASMC가 설치되지 않았습니다. 설치를 시작합니다...")

        # 3. 의존성 설치
        install_dependencies()

        # 4. ASMC 설치
        if not install_asmc():
            print("\n설치에 실패했습니다. 다음 명령을 수동으로 실행해보세요:")
            print("  pip install -e .")
            return

        # 5. 재확인
        if not verify_installation():
            print("\n설치 후에도 모듈을 찾을 수 없습니다.")
            print("Python 환경을 확인하세요.")
            return

    # 6. 실행_예제.py 실행
    print("\n실행 옵션:")
    print("1. 실행_예제.py 실행")
    print("2. ASMC 도움말 보기")
    print("3. 종료")

    choice = input("\n선택 (1-3): ").strip()

    if choice == "1":
        run_example_script()
    elif choice == "2":
        subprocess.run([sys.executable, "-m", "asmc.run_asmc", "--help"])
    else:
        print("종료합니다.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n사용자에 의해 중단되었습니다.")
    except Exception as e:
        print(f"\n오류 발생: {e}")
        import traceback
        traceback.print_exc()