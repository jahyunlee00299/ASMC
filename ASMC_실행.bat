@echo off
REM ============================================
REM ASMC 간편 실행 스크립트 (Windows)
REM ============================================
chcp 65001 >nul
echo.
echo ============================================
echo ASMC (Active Site Motif Clustering) 실행기
echo ============================================
echo.

:menu
echo 실행할 작업을 선택하세요:
echo.
echo 1. ASMC 도움말 보기
echo 2. 테스트 실행 (예제 파일 필요)
echo 3. 서열 기반 분석 실행
echo 4. 구조 모델 기반 분석 실행
echo 5. 서열 유사도 계산
echo 6. 테스트 파일 생성
echo 7. 프로젝트 테스트 실행
echo 0. 종료
echo.
set /p choice=선택 (0-7):

if "%choice%"=="1" goto help
if "%choice%"=="2" goto test_run
if "%choice%"=="3" goto seq_run
if "%choice%"=="4" goto model_run
if "%choice%"=="5" goto identity
if "%choice%"=="6" goto create_test
if "%choice%"=="7" goto run_tests
if "%choice%"=="0" goto end
echo 잘못된 선택입니다.
goto menu

:help
echo.
echo ASMC 도움말을 표시합니다...
echo.
python -m asmc.run_asmc --help
echo.
pause
goto menu

:test_run
echo.
echo 간단한 테스트 실행 (예제 파일 필요)
echo.
if not exist "test_refs.txt" (
    echo test_refs.txt 파일이 없습니다. 먼저 테스트 파일을 생성하세요.
    pause
    goto menu
)
python -m asmc.run_asmc run -s test_seqs.fasta -r test_refs.txt -o test_output/ -t 2 --id 30
echo.
echo 실행 완료! test_output/ 디렉토리를 확인하세요.
pause
goto menu

:seq_run
echo.
set /p seq_file=서열 파일 경로 입력 (FASTA 형식):
set /p ref_file=참조 구조 파일 목록 경로:
set /p out_dir=출력 디렉토리 (기본: output/):
if "%out_dir%"=="" set out_dir=output/
echo.
echo 실행 중...
python -m asmc.run_asmc run -s "%seq_file%" -r "%ref_file%" -o "%out_dir%" -t 4 --id 30
echo.
echo 실행 완료!
pause
goto menu

:model_run
echo.
set /p model_file=모델 파일 목록 경로:
set /p ref_file=참조 구조 파일 목록 경로:
set /p out_dir=출력 디렉토리 (기본: output/):
if "%out_dir%"=="" set out_dir=output/
echo.
echo 실행 중...
python -m asmc.run_asmc run -m "%model_file%" -r "%ref_file%" -o "%out_dir%"
echo.
echo 실행 완료!
pause
goto menu

:identity
echo.
set /p seq_file=서열 파일 경로 (FASTA 형식):
set /p ref_file=참조 구조 파일 목록 경로:
set /p out_file=출력 파일 (기본: identity.txt):
if "%out_file%"=="" set out_file=identity.txt
echo.
echo 서열 유사도 계산 중...
python -m asmc.run_asmc identity -s "%seq_file%" -r "%ref_file%" -o "%out_file%"
echo.
echo 완료! %out_file% 파일을 확인하세요.
pause
goto menu

:create_test
echo.
echo 테스트 파일을 생성합니다...
echo.

echo # 테스트용 참조 구조 파일 > test_refs.txt
echo # 실제 PDB 파일 경로로 교체하세요 >> test_refs.txt
echo /path/to/reference1.pdb >> test_refs.txt
echo /path/to/reference2.pdb >> test_refs.txt

(
echo ^>Test_Protein_1
echo MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
echo LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA
echo.
echo ^>Test_Protein_2
echo MKHLWFFLLLVAAPRWVLSAAGACGQEARPEAVGQHWEALGRFWDYLRWVQTLSEQVQEE
echo LLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA
) > test_seqs.fasta

echo # 테스트용 모델 파일 목록 > test_models.txt
echo # 실제 파일 경로로 교체하세요 >> test_models.txt
echo /path/to/model1.pdb	reference1 >> test_models.txt
echo /path/to/model2.pdb	reference2 >> test_models.txt

echo.
echo 테스트 파일 생성 완료:
echo - test_refs.txt
echo - test_seqs.fasta
echo - test_models.txt
echo.
echo 실제 데이터 경로로 수정 후 사용하세요!
echo.
pause
goto menu

:run_tests
echo.
echo 프로젝트 테스트를 실행합니다...
echo.
python -m pytest tests/ -v
echo.
pause
goto menu

:end
echo.
echo 프로그램을 종료합니다.
echo.
exit