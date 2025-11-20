# UDH Active Site Finder - 사용 예제

AtUdh 구조와 active site를 이용하여 다른 UDH들의 active site를 찾는 방법

## 준비사항

### 1. AtUdh 참조 구조 다운로드

```bash
# 방법 1: wget 사용
wget https://files.rcsb.org/download/3RFV.pdb -O AtUdh_3rfv.pdb

# 방법 2: curl 사용
curl https://files.rcsb.org/download/3RFV.pdb -o AtUdh_3rfv.pdb

# 방법 3: 웹 브라우저로 다운로드
# https://www.rcsb.org/structure/3RFV 에서 "Download Files" > "PDB Format" 클릭
```

### 2. BioPython 설치 확인

```bash
python -c "import Bio; print(Bio.__version__)"

# 설치되어 있지 않다면:
pip install biopython
```

### 3. 파일 구조 확인

```
pdb_structures/
├── AtUdh_3rfv.pdb                    # 참조 구조 (다운로드 필요)
├── atudh_active_site.txt             # Active site 정의 (이미 생성됨)
├── find_udh_active_sites.py          # 메인 스크립트
├── other_udh_structures/             # 비교할 다른 UDH 구조들 넣는 디렉토리
│   ├── udh_species1.pdb
│   ├── udh_species2.pdb
│   └── ...
└── results/                          # 결과 저장 디렉토리
```

## 사용 방법

### 예제 1: 단일 UDH 구조 분석

```bash
# 다른 UDH 구조 하나의 active site 찾기
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target other_udh_structures/udh_species1.pdb \
    --output results/
```

**출력:**
- `results/udh_species1_active_sites.txt`: 발견된 active site 잔기 정보

### 예제 2: 여러 UDH 구조 일괄 분석

```bash
# other_udh_structures/ 디렉토리의 모든 PDB 파일 처리
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir other_udh_structures/ \
    --output results/
```

**출력:**
- `results/[각_구조]_active_sites.txt`: 각 구조의 active site 정보
- `results/active_sites_comparison.tsv`: 모든 구조의 active site 비교표

### 예제 3: 커스텀 파라미터 사용

```bash
# Distance cutoff를 3.5 Å으로 줄여서 더 엄격하게 매칭
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir other_udh_structures/ \
    --output results/ \
    --cutoff 3.5
```

### 예제 4: 특정 Chain 지정

```bash
# Target 구조에 여러 chain이 있을 때 특정 chain 지정
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target other_udh_structures/udh_dimer.pdb \
    --target-chain B \
    --output results/
```

## 출력 파일 형식

### 개별 결과 파일 (예: udh_species1_active_sites.txt)

```
# Active sites found in udh_species1
# Reference: AtUdh_3rfv.pdb
# Alignment RMSD: 1.23 Å
# Distance cutoff: 4.0 Å
# Chain: A
#
# Format: PDB_file[TAB]Chain[TAB]Residue_numbers
#
udh_species1.pdb	A	12,14,36,38,60,82,104,106,108,137,139,141,163,165,167,189,191,213,235,237,259,261,283

# Detailed residue information:
# ResNum	ResName	Distance(Å)
# 12	THR	0.52
# 14	ASP	0.73
# 36	SER	1.12
# ...
```

### 비교 요약 파일 (active_sites_comparison.tsv)

```
Structure	Chain	RMSD(Å)	N_Sites	Residue_Positions	Residue_Types
AtUdh_3rfv	A	0.00	23	10,12,34,36,58,...	-
udh_species1	A	1.23	22	12,14,36,38,60,...	THR,ASP,SER,...
udh_species2	A	2.45	21	11,13,35,37,59,...	THR,ASP,ALA,...
```

## ASMC 워크플로우와 통합하기

찾은 active site 정보를 ASMC 파이프라인에서 사용하려면:

### 1단계: Active site 발견 (이 스크립트 사용)

```bash
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir other_udh_structures/ \
    --output results/
```

### 2단계: ASMC pocket 형식으로 변환

`results/active_sites_comparison.tsv`를 ASMC의 pocket.txt 형식으로 변환:

```bash
# Python을 사용하여 변환
python << EOF
import pandas as pd

# TSV 읽기
df = pd.read_csv('results/active_sites_comparison.tsv', sep='\t')

# ASMC pocket 형식으로 변환
with open('udh_pocket.txt', 'w') as f:
    for _, row in df.iterrows():
        pdb_name = row['Structure']
        chain = row['Chain']
        residues = row['Residue_Positions']
        # 파일 경로 추가 (실제 경로로 수정 필요)
        f.write(f"pdb_structures/{pdb_name}.pdb\t{chain}\t{residues}\n")

print("Created udh_pocket.txt for ASMC")
EOF
```

### 3단계: ASMC 실행

```bash
# 발견된 active site를 사용하여 ASMC 클러스터링
python -m asmc.run_asmc run \
    -m models.txt \
    -r udh_references.txt \
    -p udh_pocket.txt \
    -o udh_asmc_results/ \
    --chain A
```

## Active Site 정의 파일 수정하기

실제 PDB 구조를 확인한 후 `atudh_active_site.txt`를 정확한 잔기 번호로 수정하세요:

### 방법 1: P2RANK로 자동 검출

```bash
# P2RANK를 사용하여 pocket 자동 검출
# (ASMC에 통합되어 있음)
python -m asmc.run_asmc run \
    -r AtUdh_3rfv.pdb \
    -o temp_pocket/ \
    --end pocket
```

### 방법 2: PyMOL로 수동 확인

```python
# PyMOL에서 실행
load AtUdh_3rfv.pdb
select active_site, resi 10+12+34+36+58+80+102+104+106+135+137+139+161+163+165+187+189+211+233+235+257+259+281
show sticks, active_site
color red, active_site
```

### 방법 3: 문헌 참조

UDP-glucose dehydrogenase의 active site에 대한 문헌:
- Campbell et al. (1997) Plant Cell
- Sommer et al. (2004) J. Biol. Chem.
- 각 UDH의 특정 연구 논문

## 고급 사용법

### Distance Cutoff 조정

- **Default (4.0 Å)**: 일반적인 구조 정렬에 적합
- **Strict (3.0-3.5 Å)**: 더 정확한 매칭, 일부 residue 누락 가능
- **Loose (5.0-6.0 Å)**: 더 많은 residue 포함, false positive 가능

```bash
# 엄격한 매칭
python find_udh_active_sites.py ... --cutoff 3.0

# 느슨한 매칭
python find_udh_active_sites.py ... --cutoff 5.0
```

### 특정 패턴의 파일만 처리

```bash
# 특정 이름 패턴만
python find_udh_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir other_udh_structures/ \
    --pattern "*_model*.pdb" \
    --output results/
```

## 문제 해결

### 1. "No chains found" 오류

**원인**: PDB 파일 형식 문제 또는 빈 파일

**해결**:
```bash
# PDB 파일 확인
head -20 your_structure.pdb

# Chain 정보 확인
grep "^ATOM" your_structure.pdb | awk '{print $5}' | sort -u
```

### 2. "No active site atoms found" 오류

**원인**: atudh_active_site.txt의 잔기 번호가 실제 PDB와 맞지 않음

**해결**:
```bash
# AtUdh PDB에서 실제 잔기 번호 확인
grep "^ATOM.*CA.*A" AtUdh_3rfv.pdb | awk '{print $6}' | head -20
```

### 3. RMSD가 너무 높음 (>5 Å)

**원인**: 구조가 너무 다르거나 잘못된 구조

**확인**:
- 정말 같은 protein family인지 확인
- PDB 파일 품질 확인 (homology model vs. crystal structure)
- Chain이 올바른지 확인

### 4. Active site가 너무 적게 발견됨

**원인**: Distance cutoff가 너무 엄격하거나 정렬이 불완전

**해결**:
```bash
# Cutoff를 늘려서 재시도
python find_udh_active_sites.py ... --cutoff 5.0
```

## 추가 리소스

- **ASMC 문서**: `../ASMC_사용법.md`
- **PDB 데이터베이스**: https://www.rcsb.org/
- **P2RANK (pocket detection)**: https://github.com/rdk/p2rank
- **BioPython 문서**: https://biopython.org/wiki/Documentation

## 연락처

문제가 발생하면 ASMC GitHub repository에 이슈를 등록하세요.
