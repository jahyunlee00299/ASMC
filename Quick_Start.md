# ASMC Quick Start Guide (빠른 시작 가이드)

## 🚀 5분 안에 시작하기

### 1단계: 설치 확인
```bash
# ASMC가 제대로 설치되었는지 확인
python -m asmc.run_asmc --help
```

### 2단계: 테스트 파일 생성
```bash
# Python 스크립트로 테스트 파일 생성 및 실행
python 실행_예제.py
# 대화형 메뉴에서 원하는 옵션 선택
```

또는 Windows에서:
```bash
# 배치 파일 실행
ASMC_실행.bat
# 옵션 6 선택 (테스트 파일 생성)
```

### 3단계: 첫 번째 실행
```bash
# 가장 간단한 실행 예제
python -m asmc.run_asmc run --help
```

## 📝 가장 자주 사용되는 명령어

### 1. 서열만 있을 때 (구조 예측 + 분석)
```bash
python -m asmc.run_asmc run \
    -s my_sequences.fasta \
    -r reference_pdbs.txt \
    -o results/
```

### 2. 구조가 이미 있을 때
```bash
python -m asmc.run_asmc run \
    -m my_models.txt \
    -r reference_pdbs.txt \
    -o results/
```

### 3. 서열 유사도만 빠르게 확인
```bash
python -m asmc.run_asmc identity \
    -s sequences.fasta \
    -r references.txt \
    -o similarity.txt
```

## 💡 실용적인 팁

1. **메모리 부족시**: 스레드 수 줄이기 `-t 2`
2. **빠른 테스트**: `--end pocket`으로 포켓 검출만
3. **클러스터링 조정**: `-e auto --min-samples auto`로 자동 설정
4. **시각화**: `--format svg`로 벡터 이미지 생성

## 📂 필수 파일 형식 예시

### references.txt
```
./structures/1ABC.pdb
./structures/2DEF.pdb
```

### sequences.fasta
```
>Protein1
MKVLWAALLVTFLAGCQA...
>Protein2
MKHLWFFLLLVAAPRWVL...
```

### models.txt (탭으로 구분)
```
./models/model1.pdb	1ABC
./models/model2.pdb	2DEF
```

## ⚠️ 주의사항

- PDB 파일 경로는 절대 경로 또는 상대 경로 사용 가능
- FASTA 파일은 표준 형식 준수
- 출력 디렉토리는 자동 생성됨
- Windows에서는 경로에 백슬래시(\) 대신 슬래시(/) 사용 권장

## 🆘 도움이 필요하면

1. 상세 사용법: `ASMC_사용법.md` 파일 참조
2. 대화형 실행: `python 실행_예제.py` 실행
3. Windows 사용자: `ASMC_실행.bat` 실행
4. 각 명령어 도움말: `python -m asmc.run_asmc --help` 옵션 사용

## 📚 관련 파일

- **실행_예제.py**: 대화형 메뉴로 다양한 ASMC 기능 테스트
- **run_asmc_test.py**: 기본적인 실행 테스트 스크립트
- **run_asmc_demo.py**: 시각화 포함 데모 스크립트
- **ASMC_실행.bat**: Windows용 배치 스크립트

---
작성일: 2025년 11월 4일