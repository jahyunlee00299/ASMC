# UDH Active Site 추출 및 Superposition 가이드

AtUdh를 참조로 5개 UDH AlphaFold 구조의 active site를 추출하고, **active site만을 사용하여 superposition**하는 방법

## 대상 구조

Rossmann fold RMSD 기준으로 선택된 5개 구조:

| 구조 ID      | Rossmann RMSD | 의미       |
|-------------|---------------|-----------|
| A0A1I2LZE5  | 15.78 Å       | 가장 유사  |
| A0A1I1AQL9  | 17.27 Å       | -         |
| A0A2Z4AB20  | 19.85 Å       | -         |
| A0A1N6JJV2  | 20.89 Å       | -         |
| A0A1U7CQA8  | 25.19 Å       | 가장 다름  |

## 준비사항

### 1. 디렉토리 구조 확인

```bash
pdb_structures/
├── AtUdh_3rfv.pdb                          # 참조 구조
├── atudh_active_site.txt                   # Active site 정의
├── extract_and_align_active_sites.py       # 메인 스크립트
└── uronate_dehydrogenase_alphafold/        # AlphaFold UDH 구조들
    ├── A0A1I2LZE5_xxxxx.pdb
    ├── A0A1I1AQL9_xxxxx.pdb
    ├── A0A2Z4AB20_xxxxx.pdb
    ├── A0A1N6JJV2_xxxxx.pdb
    └── A0A1U7CQA8_xxxxx.pdb
```

### 2. PDB 파일 복사

Windows에서 AlphaFold 구조들을 복사:

```bash
# uronate_dehydrogenase_alphafold 폴더의 PDB 파일들을 복사
# 파일명에 A0A1I2LZE5, A0A1I1AQL9 등의 ID가 포함되어 있어야 함
```

Linux에서:
```bash
# 예시
cp /path/to/alphafold/structures/*.pdb pdb_structures/uronate_dehydrogenase_alphafold/
```

## 실행 방법

### 기본 실행

```bash
cd pdb_structures

python extract_and_align_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir uronate_dehydrogenase_alphafold/ \
    --output active_site_results/
```

### 커스텀 Distance Cutoff

```bash
# Active site 찾을 때 더 엄격한 기준 사용 (3.0 Å)
python extract_and_align_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir uronate_dehydrogenase_alphafold/ \
    --output active_site_results/ \
    --cutoff 3.0
```

## 스크립트의 동작 방식

### Step 1: Initial Global Alignment
모든 CA 원자를 사용하여 초기 구조 정렬
```
→ Global RMSD 계산
```

### Step 2: Active Site Identification
1. 참조 구조의 active site 잔기 추출
2. 초기 정렬 기준으로 target에서 대응 잔기 탐색
3. Distance cutoff 이내의 잔기만 선택

```
Reference SER78 -> Target SER80 (1.48 Å) ✓
Reference ALA101 -> Target SER103 (1.24 Å) ✓
...
```

### Step 3: Active Site Superposition
**Active site CA 원자만을 사용**하여 재정렬
```
→ Active Site RMSD 계산
→ 전체 구조에 transformation 적용
```

### Step 4: Results
- Active site RMSD vs. Global RMSD 비교
- Improvement 계산
- Aligned 구조 저장

## 출력 결과

### 파일 구조
```
active_site_results/
├── active_site_alignment_summary.tsv       # 전체 요약
├── A0A1I2LZE5_active_site_alignment.txt   # 각 구조의 상세 결과
├── A0A1I2LZE5_aligned.pdb                 # Active site로 정렬된 구조
├── A0A1I1AQL9_active_site_alignment.txt
├── A0A1I1AQL9_aligned.pdb
├── ...
```

### Summary TSV 형식

```tsv
Structure_ID    Rossmann_RMSD(Å)    Active_Site_RMSD(Å)    Global_RMSD(Å)    Improvement(Å)
A0A1I2LZE5     15.78                2.34                   12.45             10.11
A0A1I1AQL9     17.27                3.12                   14.23             11.11
A0A2Z4AB20     19.85                4.56                   16.78             12.22
A0A1N6JJV2     20.89                5.23                   18.34             13.11
A0A1U7CQA8     25.19                7.89                   22.45             14.56
```

### 개별 결과 파일 예시

`A0A1I2LZE5_active_site_alignment.txt`:
```
# Active Site Alignment Results
# Target: A0A1I2LZE5
# Rossmann RMSD: 15.78 Å
# Note: Most similar
#
# Chain: A
# Active sites found: 21/23
# Global RMSD: 12.45 Å
# Active site RMSD: 2.34 Å
# Improvement: 10.11 Å
#
# Correspondences:
# Ref_Res    Target_Res    Distance(Å)
SER78       SER80         1.48
ALA101      SER103        1.24
GLU102      GLU104        1.51
...
```

## 결과 분석

### 예상되는 패턴

1. **Rossmann RMSD가 낮을수록**:
   - Active site RMSD도 낮을 것으로 예상
   - 더 많은 active site 잔기 매칭
   - 잔기 타입 보존도 높음

2. **Active Site vs. Global RMSD**:
   - Active site RMSD < Global RMSD (일반적)
   - Active site가 더 잘 보존됨
   - Improvement 값이 클수록 active site가 구조적으로 중요

3. **구조별 특성**:
   - **A0A1I2LZE5** (가장 유사): 높은 active site 복구율, 낮은 RMSD
   - **A0A1U7CQA8** (가장 다름): 낮은 복구율, 높은 RMSD

## 시각화 (PyMOL)

### 모든 정렬된 구조 로드

```python
# PyMOL 실행
load active_site_results/A0A1I2LZE5_aligned.pdb
load active_site_results/A0A1I1AQL9_aligned.pdb
load active_site_results/A0A2Z4AB20_aligned.pdb
load active_site_results/A0A1N6JJV2_aligned.pdb
load active_site_results/A0A1U7CQA8_aligned.pdb
load AtUdh_3rfv.pdb

# 모든 구조를 cartoon으로
as cartoon, all

# 각 구조를 다른 색으로
color green, AtUdh_3rfv
color cyan, A0A1I2LZE5_aligned
color blue, A0A1I1AQL9_aligned
color yellow, A0A2Z4AB20_aligned
color orange, A0A1N6JJV2_aligned
color red, A0A1U7CQA8_aligned
```

### Active Site만 표시

```python
# AtUdh의 active site 표시
select atudh_active, AtUdh_3rfv and resi 10+12+34+36+58+80+102+104+106+135+137+139+161+163+165+187+189+211+233+235+257+259+281
show sticks, atudh_active
color red, atudh_active

# 각 구조의 active site 표시 (개별 결과 파일에서 잔기 번호 확인)
# 예: A0A1I2LZE5의 경우
select a0a1i2lze5_active, A0A1I2LZE5_aligned and resi 12+14+36+38+...
show sticks, a0a1i2lze5_active
color cyan, a0a1i2lze5_active
```

### RMSD 측정

```python
# Active site 간 RMSD 계산
rms_cur atudh_active, a0a1i2lze5_active

# 특정 chain만
rms_cur AtUdh_3rfv and chain A and name CA, A0A1I2LZE5_aligned and chain A and name CA
```

## 추가 분석

### 1. 잔기 보존도 분석

각 결과 파일에서 잔기 타입 변화를 확인:

```bash
# 잔기 타입 추출
grep -v "^#" active_site_results/A0A1I2LZE5_active_site_alignment.txt | awk '{print $1, $2}'
```

보존된 잔기 vs. 변이된 잔기 분석

### 2. Distance 분포 분석

```python
import pandas as pd

# 각 구조의 distance 분포 확인
for struct_id in ['A0A1I2LZE5', 'A0A1I1AQL9', 'A0A2Z4AB20', 'A0A1N6JJV2', 'A0A1U7CQA8']:
    df = pd.read_csv(f'active_site_results/{struct_id}_active_site_alignment.txt',
                     sep='\t', comment='#', header=None, names=['Ref', 'Target', 'Distance'])
    print(f"{struct_id}:")
    print(f"  Mean distance: {df['Distance'].mean():.2f} Å")
    print(f"  Max distance: {df['Distance'].max():.2f} Å")
    print(f"  Min distance: {df['Distance'].min():.2f} Å")
```

### 3. RMSD vs. Rossmann RMSD 상관관계

```python
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('active_site_results/active_site_alignment_summary.tsv', sep='\t')

plt.figure(figsize=(10, 6))
plt.scatter(df['Rossmann_RMSD(Å)'], df['Active_Site_RMSD(Å)'], s=100)
plt.xlabel('Rossmann RMSD (Å)', fontsize=12)
plt.ylabel('Active Site RMSD (Å)', fontsize=12)
plt.title('Active Site RMSD vs. Rossmann RMSD', fontsize=14)

for i, row in df.iterrows():
    plt.annotate(row['Structure_ID'],
                (row['Rossmann_RMSD(Å)'], row['Active_Site_RMSD(Å)']),
                fontsize=9)

plt.grid(True, alpha=0.3)
plt.savefig('active_site_results/rmsd_correlation.png', dpi=300)
plt.show()
```

## ASMC 통합

### Active Site 기반 Clustering

```bash
# 1. Active site 정보를 ASMC pocket 형식으로 변환
python << EOF
import pandas as pd

df = pd.read_csv('active_site_results/active_site_alignment_summary.tsv', sep='\t')

with open('udh_alphafold_pocket.txt', 'w') as f:
    f.write("# UDH AlphaFold Active Sites\n")
    f.write("# Extracted using active site superposition\n\n")

    for _, row in df.iterrows():
        struct_id = row['Structure_ID']
        # 각 결과 파일에서 잔기 번호 추출
        # (실제 구현시 파일 파싱 필요)
        f.write(f"uronate_dehydrogenase_alphafold/{struct_id}_aligned.pdb\tA\t...\n")
EOF
```

### 2. ASMC 클러스터링 실행

```bash
python -m asmc.run_asmc run \
    -m udh_alphafold_models.txt \
    -r udh_references.txt \
    -p udh_alphafold_pocket.txt \
    -o udh_alphafold_asmc_results/
```

## 문제 해결

### 1. "No PDB file found for XXX"

**원인**: 파일명에 구조 ID가 포함되지 않음

**해결**:
```bash
# 파일명 확인
ls uronate_dehydrogenase_alphafold/

# 파일명에 ID가 있는지 확인
# 예: AF-A0A1I2LZE5-F1-model_v4.pdb (OK)
# 예: structure_001.pdb (NG)

# 필요시 파일명 변경
mv structure_001.pdb A0A1I2LZE5_model.pdb
```

### 2. "Too few active site residues for alignment"

**원인**: Distance cutoff가 너무 엄격

**해결**:
```bash
# Cutoff를 늘려서 재시도
python extract_and_align_active_sites.py ... --cutoff 5.0
```

### 3. Active Site RMSD > Global RMSD

**의미**: Active site가 전체 구조보다 덜 보존됨 (unusual)

**확인**:
- 참조 구조의 active site 정의가 정확한지
- Chain ID가 올바른지
- 구조가 정말 같은 family인지

## 참고 자료

- **스크립트**: `extract_and_align_active_sites.py`
- **참조 구조**: AtUdh (PDB: 3RFV)
- **BioPython 문서**: https://biopython.org/wiki/Superimpose
- **ASMC 문서**: `../ASMC_사용법.md`

## 요약

이 도구는:
✅ **Active site만 추출**하여 정렬
✅ **5개 UDH 구조 자동 처리**
✅ **Global vs. Active site RMSD 비교**
✅ **상세한 잔기별 대응 정보 제공**
✅ **PyMOL 시각화 가능**
✅ **ASMC 통합 가능**

Active site 기반 superposition을 통해 기능적으로 중요한 부위의 구조적 보존도를 정확히 분석할 수 있습니다!
