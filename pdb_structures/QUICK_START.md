# UDH Active Site Finder - Quick Start Guide

AtUdh를 참조로 다른 UDH들의 active site를 빠르게 찾는 방법

## 5분 안에 시작하기

### 1. AtUdh 구조 다운로드 (1분)

```bash
cd pdb_structures
wget https://files.rcsb.org/download/3RFV.pdb -O AtUdh_3rfv.pdb

# wget이 안되면:
curl https://files.rcsb.org/download/3RFV.pdb -o AtUdh_3rfv.pdb
```

또는 브라우저에서 다운로드:
- https://www.rcsb.org/structure/3RFV
- "Download Files" → "PDB Format" 클릭
- `AtUdh_3rfv.pdb`로 저장

### 2. 비교할 UDH 구조 준비 (1분)

```bash
# other_udh_structures 디렉토리에 PDB 파일들 복사
cp /path/to/your/udh_files/*.pdb other_udh_structures/

# 또는 직접 다운로드
cd other_udh_structures
# 예: Human UGDH (PDB: 2Q3E)
wget https://files.rcsb.org/download/2Q3E.pdb
# 예: E. coli UGDH (PDB: 3PTZ)
wget https://files.rcsb.org/download/3PTZ.pdb
cd ..
```

### 3. Active Site 찾기 실행 (2분)

```bash
# 모든 UDH 구조 분석
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir other_udh_structures/ \
    --output results/
```

### 4. 결과 확인 (1분)

```bash
# 요약 파일 확인
cat results/active_sites_comparison.tsv

# 개별 결과 확인
ls results/*_active_sites.txt
```

## 완료! 🎉

이제 `results/` 디렉토리에 다음 파일들이 생성됩니다:
- `active_sites_comparison.tsv` - 모든 구조의 active site 비교
- `[구조명]_active_sites.txt` - 각 구조의 상세 정보

## 다음 단계

### ASMC 파이프라인과 통합

```bash
# 1. ASMC pocket 형식으로 변환
python convert_to_asmc_pocket.py \
    results/active_sites_comparison.tsv \
    udh_pocket.txt

# 2. ASMC 실행 (예시)
python -m asmc.run_asmc run \
    -m models.txt \
    -r references.txt \
    -p udh_pocket.txt \
    -o asmc_results/
```

### 결과 시각화 (PyMOL)

```bash
# PyMOL 스크립트로 active site 시각화
pymol AtUdh_3rfv.pdb other_udh_structures/*.pdb
```

PyMOL에서:
```python
# 첫 번째 결과 파일에서 잔기 번호 복사
select active_site, resi 10+12+34+36+58+80+102+104+106+135+137+139+161+163+165+187+189+211+233+235+257+259+281
show sticks, active_site
color red, active_site
```

## 자주 하는 질문

**Q: AtUdh 구조를 다운로드할 수 없어요**
```bash
# 방법 1: Python으로 다운로드
python download_pdb.py 3RFV .

# 방법 2: 브라우저 사용
# https://www.rcsb.org/structure/3RFV
```

**Q: 내 UDH 구조 파일 형식이 맞는지 확인하려면?**
```bash
# PDB 파일 헤더 확인
head -20 your_udh.pdb

# ATOM 레코드가 있어야 함
grep "^ATOM" your_udh.pdb | head -5
```

**Q: Active site가 너무 적게 나와요**
```bash
# Distance cutoff를 늘려서 재시도
python find_udh_active_sites.py ... --cutoff 5.0
```

**Q: RMSD가 너무 높아요 (>5 Å)**
- 정말 같은 protein family인지 확인
- Sequence identity 확인
- Chain이 올바른지 확인

## 도움말

더 자세한 사용법은 다음 파일들을 참조하세요:
- `USAGE_EXAMPLE.md` - 상세한 사용 예제
- `README.md` - 프로젝트 개요
- `../ASMC_사용법.md` - ASMC 전체 워크플로우

## 문제 해결

스크립트가 실행되지 않으면:

```bash
# BioPython 설치 확인
python -c "import Bio; print('OK')"

# 설치되어 있지 않다면
pip install biopython

# 실행 권한 확인
chmod +x find_udh_active_sites.py

# 도움말 확인
python find_udh_active_sites.py --help
```
