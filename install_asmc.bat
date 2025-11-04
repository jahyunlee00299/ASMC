@echo off
chcp 65001 >nul
echo.
echo ============================================
echo ASMC 설치 스크립트
echo ============================================
echo.

REM Python 버전 확인
echo 🐍 Python 버전 확인 중...
python --version

echo.
echo ASMC를 현재 Python 환경에 설치합니다.
echo.

REM 기존 설치 제거 (있을 경우)
echo 1. 기존 설치 제거 중...
pip uninstall -y asmc 2>nul

REM 의존성 설치
echo.
echo 2. 의존성 패키지 설치 중...
echo    (시간이 걸릴 수 있습니다)
pip install biopython>=1.81
pip install numpy
pip install scikit-learn
pip install scipy
pip install pyyaml
pip install pillow
pip install plotnineseqsuite
pip install pytest
pip install xlsxwriter

REM ASMC 설치 (개발 모드)
echo.
echo 3. ASMC 설치 중 (개발 모드)...
pip install -e .

echo.
echo ============================================
echo 설치 확인 중...
echo ============================================
python -c "import asmc; print('✅ ASMC 모듈이 성공적으로 설치되었습니다!')" 2>nul
if errorlevel 1 (
    echo ❌ 설치 실패! 오류를 확인하세요.
) else (
    echo.
    echo 🎉 설치 완료!
    echo.
    echo 이제 다음 명령으로 ASMC를 사용할 수 있습니다:
    echo - python 실행_예제.py
    echo - python -m asmc.run_asmc --help
)

echo.
pause