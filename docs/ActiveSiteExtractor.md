# Active Site Extractor

단백질 구조와 서열로부터 active site 잔기를 추출하는 통합 도구입니다.

## 개요

`ActiveSiteExtractor`는 단백질의 active site 잔기를 추출하기 위한 두 가지 방법을 제공합니다:

1. **구조 정렬 기반 추출** (PDB 구조)
   - 참조 구조와 타겟 구조를 3D 공간에서 정렬
   - 공간적 근접성을 기반으로 대응하는 잔기 식별
   - RMSD 및 거리 정보 제공

2. **서열 정렬 기반 추출** (FASTA 서열)
   - 참조 서열과 타겟 서열을 정렬
   - 정렬 결과를 바탕으로 active site 위치의 잔기 추출
   - 대량의 서열 배치 처리 가능

## 주요 기능

- 단일 또는 다중 구조/서열 처리
- PDB 구조 및 FASTA 서열 모두 지원
- 자동 chain 감지
- 커스터마이징 가능한 거리 cutoff
- 배치 처리 기능
- 상세한 결과 리포트 생성
- TSV 형식의 요약 파일 생성

## 설치 및 요구사항

### 필수 패키지

```bash
pip install biopython numpy
```

### Python 버전
- Python 3.6 이상

## 사용 방법

### 1. Python 모듈로 사용

#### 기본 초기화

```python
from asmc.active_site_extractor import ActiveSiteExtractor

# 참조 구조와 active site 정의로 초기화
extractor = ActiveSiteExtractor(
    reference_pdb='path/to/reference.pdb',
    reference_sites='path/to/sites.txt',  # 또는 [137, 138, 139, ...]
    reference_chain='A',
    distance_cutoff=4.0,
    verbose=True
)
```

#### PDB 구조에서 추출

```python
# 단일 PDB 구조 처리
result = extractor.extract_from_structure(
    target_pdb='path/to/target.pdb',
    target_chain='A',  # None이면 자동 감지
    output_file='output/results.txt'
)

# 결과 확인
print(f"RMSD: {result['rmsd']:.2f} Å")
print(f"Found {result['n_sites']} active sites")
print(f"Residues: {result['residue_numbers']}")
```

#### FASTA 서열에서 추출

```python
# FASTA 파일의 서열들 처리
active_site_seqs = extractor.extract_from_sequences(
    fasta_file='path/to/sequences.fasta',
    output_file='output/active_sites.fasta',
    progress_interval=100
)

print(f"Processed {len(active_site_seqs)} sequences")
```

#### 배치 처리

```python
# 디렉토리의 모든 PDB 파일 처리
results = extractor.batch_process_structures(
    structure_dir='path/to/pdb_directory/',
    output_dir='output/',
    pattern='*.pdb'
)

# 요약 파일이 output/active_sites_summary.tsv에 생성됩니다
```

### 2. 명령줄 도구로 사용

#### 단일 PDB 구조 처리

```bash
python -m asmc.active_site_extractor \
    --reference-pdb reference.pdb \
    --reference-sites sites.txt \
    --target-pdb target.pdb \
    --output results.txt
```

#### FASTA 서열 처리

```bash
python -m asmc.active_site_extractor \
    --reference-pdb reference.pdb \
    --reference-sites "137,138,139,140,165,166,189" \
    --target-fasta sequences.fasta \
    --output active_sites.fasta
```

#### 배치 처리

```bash
python -m asmc.active_site_extractor \
    --reference-pdb reference.pdb \
    --reference-sites sites.txt \
    --target-dir structures/ \
    --output results/ \
    --distance-cutoff 3.5
```

#### 명령줄 옵션

| 옵션 | 설명 | 기본값 |
|------|------|--------|
| `--reference-pdb`, `-r` | 참조 PDB 파일 (필수) | - |
| `--reference-sites`, `-s` | Active site 정의 파일 또는 쉼표로 구분된 번호 (필수) | - |
| `--reference-chain`, `-c` | 참조 구조의 chain ID | A |
| `--target-pdb`, `-t` | 타겟 PDB 파일 (단일 구조) | - |
| `--target-fasta`, `-f` | 타겟 FASTA 파일 (서열) | - |
| `--target-dir`, `-d` | 타겟 PDB 파일들이 있는 디렉토리 | - |
| `--output`, `-o` | 출력 파일 또는 디렉토리 | active_sites_output |
| `--distance-cutoff` | 구조 매칭을 위한 거리 cutoff (Å) | 4.0 |
| `--quiet` | 진행 메시지 숨기기 | False |

## Active Site 정의 파일 형식

Active site를 정의하는 파일은 탭으로 구분된 형식을 사용합니다:

```text
# Comments start with #
pdb_file.pdb[TAB]chain_id[TAB]residue_numbers

# Example:
AtUdh_3rfv.pdb	A	137,138,139,140,141,143,165,166,167,189,190,191
```

**형식 규칙:**
- `#`으로 시작하는 줄은 주석
- 각 필드는 탭(TAB)으로 구분
- 잔기 번호는 쉼표로 구분

## 출력 파일 형식

### 구조 추출 결과 파일

```text
# Active sites found in target_structure
# Reference: reference.pdb
# Alignment RMSD: 1.23 Å
# Distance cutoff: 4.0 Å
# Chain: A
#
# Format: PDB_file[TAB]Chain[TAB]Residue_numbers
#
target.pdb	A	142,143,144,145,170,171,195

# Detailed residue information:
# ResNum	ResName	Distance(Å)
# 142	GLY	0.85
# 143	SER	1.12
...
```

### 배치 처리 요약 파일 (TSV)

```text
Structure	Chain	RMSD(Å)	N_Sites	Residue_Positions	Residue_Types
reference	A	0.00	20	137,138,139,...	-
target1	A	1.23	19	142,143,144,...	GLY,SER,THR,...
target2	A	2.45	18	140,141,142,...	ALA,GLY,SER,...
```

### FASTA 출력

```fasta
>seq1 Active site residues (20 positions)
GSTHIVLQNDKPFYAEWCMR
>seq2 Active site residues (20 positions)
GSTHIVLQNDKPFYAEWCMR
```

## 사용 예제

### 예제 1: UDH Active Site 추출

```python
from asmc.active_site_extractor import ActiveSiteExtractor

# UDH active site 위치 정의
udh_active_sites = [137, 138, 139, 140, 141, 143, 165, 166, 167,
                    189, 190, 191, 213, 214, 215, 237, 238, 257, 258, 259]

# Extractor 초기화
extractor = ActiveSiteExtractor(
    reference_pdb='test_data/AtUdh_pdb3rfv_chainA.pdb',
    reference_sites=udh_active_sites,
    reference_chain='A',
    distance_cutoff=4.0
)

# FASTA 서열에서 active site 추출
active_sites = extractor.extract_from_sequences(
    fasta_file='test_data/UDHs_filtered_std2.5.fasta',
    output_file='udh_active_sites.fasta'
)

print(f"Extracted active sites from {len(active_sites)} sequences")
```

### 예제 2: Substrate Binding Site 추출

```python
# Substrate binding site 위치 정의
substrate_direct = [75, 111, 112, 113, 136, 165, 174, 258]
substrate_proximal = [76, 163, 164, 175]
substrate_sites = sorted(substrate_direct + substrate_proximal)

# Extractor 초기화
extractor = ActiveSiteExtractor(
    reference_pdb='test_data/AtUdh_pdb3rfv_chainA.pdb',
    reference_sites=substrate_sites,
    reference_chain='A'
)

# 서열에서 substrate binding site 추출
substrate_site_seqs = extractor.extract_from_sequences(
    fasta_file='test_data/UDHs_filtered_std2.5.fasta',
    output_file='udh_substrate_sites.fasta'
)
```

### 예제 3: 다중 구조 비교

```python
# 여러 PDB 구조의 active site 비교
extractor = ActiveSiteExtractor(
    reference_pdb='reference/AtUdh.pdb',
    reference_sites='reference/active_sites.txt',
    distance_cutoff=3.5
)

# 디렉토리의 모든 구조 처리
results = extractor.batch_process_structures(
    structure_dir='structures/',
    output_dir='comparison_results/',
    pattern='*.pdb'
)

# RMSD 기준으로 정렬
sorted_results = sorted(results, key=lambda x: x['rmsd'])

for result in sorted_results[:10]:  # Top 10
    print(f"{result['target']}: RMSD={result['rmsd']:.2f} Å, "
          f"Sites={result['n_sites']}")
```

### 예제 4: 구조와 서열 통합 분석

```python
# 1단계: 구조 정렬로 active site 식별
extractor = ActiveSiteExtractor(
    reference_pdb='reference.pdb',
    reference_sites='sites.txt'
)

# 대표 구조에서 active site 추출
representative_result = extractor.extract_from_structure(
    target_pdb='representative_homolog.pdb',
    output_file='homolog_sites.txt'
)

# 2단계: 식별된 위치를 사용하여 서열 분석
homolog_sites = representative_result['residue_numbers']

extractor2 = ActiveSiteExtractor(
    reference_pdb='representative_homolog.pdb',
    reference_sites=homolog_sites
)

# 전체 패밀리 서열에서 active site 추출
family_active_sites = extractor2.extract_from_sequences(
    fasta_file='protein_family.fasta',
    output_file='family_active_sites.fasta'
)
```

## API 레퍼런스

### ActiveSiteExtractor 클래스

#### `__init__(reference_pdb, reference_sites, reference_chain='A', distance_cutoff=4.0, verbose=True)`

ActiveSiteExtractor 객체를 초기화합니다.

**Parameters:**
- `reference_pdb` (str, optional): 참조 PDB 파일 경로
- `reference_sites` (str or List[int], optional): Active site 정의 파일 경로 또는 잔기 번호 리스트
- `reference_chain` (str): 참조 구조의 chain ID (기본값: 'A')
- `distance_cutoff` (float): 구조 매칭을 위한 거리 cutoff in Å (기본값: 4.0)
- `verbose` (bool): 진행 정보 출력 여부 (기본값: True)

#### `extract_from_structure(target_pdb, target_chain=None, output_file=None)`

PDB 구조에서 active site를 추출합니다 (구조 정렬 사용).

**Parameters:**
- `target_pdb` (str): 타겟 PDB 파일 경로
- `target_chain` (str, optional): 타겟 chain ID (None이면 자동 감지)
- `output_file` (str, optional): 결과 출력 파일 경로

**Returns:**
- `dict`: 추출 결과를 담은 딕셔너리
  - `target`: 타겟 이름
  - `chain`: Chain ID
  - `rmsd`: 정렬 RMSD
  - `n_sites`: 찾은 active site 개수
  - `active_sites`: (residue_number, residue_name, distance) 튜플 리스트
  - `residue_numbers`: 잔기 번호 리스트
  - `residue_types`: 잔기 타입 리스트
  - `distances`: 거리 리스트

#### `extract_from_sequences(fasta_file, output_file=None, progress_interval=100)`

FASTA 서열에서 active site를 추출합니다 (서열 정렬 사용).

**Parameters:**
- `fasta_file` (str): 타겟 서열 FASTA 파일 경로
- `output_file` (str, optional): 결과 FASTA 파일 경로
- `progress_interval` (int): 진행 상황 표시 간격 (기본값: 100)

**Returns:**
- `List[SeqRecord]`: Active site 서열을 담은 SeqRecord 객체 리스트

#### `batch_process_structures(structure_dir, output_dir, pattern='*.pdb')`

여러 PDB 구조를 배치로 처리합니다.

**Parameters:**
- `structure_dir` (str): PDB 파일이 있는 디렉토리
- `output_dir` (str): 결과 출력 디렉토리
- `pattern` (str): 파일 패턴 (기본값: '*.pdb')

**Returns:**
- `List[dict]`: 각 구조의 결과 딕셔너리 리스트

#### 정적 메서드

##### `read_pdb_sequence(pdb_file, chain_id='A')`

PDB 파일에서 아미노산 서열을 추출합니다.

**Returns:** `(sequence, residue_numbers)` 튜플

##### `align_sequences(seq1, seq2)`

두 서열을 전역 정렬합니다.

**Returns:** `(aligned_seq1, aligned_seq2, score)` 튜플

## 알고리즘 설명

### 구조 기반 방법

1. **구조 로딩**: 참조 및 타겟 PDB 구조 로드
2. **CA 원자 추출**: 각 구조에서 Cα 원자 좌표 추출
3. **전역 정렬**: Superimposer를 사용한 구조 정렬
4. **RMSD 계산**: 정렬 품질 평가
5. **근접 잔기 찾기**: 참조 active site 각 잔기에 대해 타겟에서 가장 가까운 잔기 식별
6. **Cutoff 적용**: 거리 cutoff 이하의 잔기만 선택

### 서열 기반 방법

1. **서열 로딩**: 참조 PDB에서 서열 추출, 타겟 FASTA 서열 로드
2. **전역 정렬**: Pairwise2를 사용한 서열 정렬
3. **위치 매핑**: 정렬 위치를 참조 잔기 번호에 매핑
4. **Active site 추출**: 참조의 active site 위치에 해당하는 타겟 잔기 추출
5. **갭 처리**: 정렬 갭 위치는 'X'로 표시

## 제한사항 및 고려사항

1. **구조 정렬**
   - 구조가 크게 다른 경우 정렬 실패 가능
   - RMSD가 높은 경우 결과 신뢰도 감소
   - Distance cutoff 조정 필요할 수 있음

2. **서열 정렬**
   - 서열 유사도가 낮으면 정렬 품질 저하
   - 삽입/결실이 많은 경우 정확도 감소
   - Active site 위치가 보존되어 있다고 가정

3. **성능**
   - 대량의 구조 처리 시 시간 소요
   - 메모리 사용량은 처리하는 파일 크기에 비례

4. **Chain 선택**
   - 다중 chain 구조의 경우 적절한 chain 선택 필요
   - 자동 감지는 첫 번째 chain 사용

## 문제 해결

### ImportError: BioPython required

```bash
pip install biopython
```

### FileNotFoundError: PDB file not found

- 파일 경로 확인
- 절대 경로 또는 올바른 상대 경로 사용

### ValueError: No chains found

- PDB 파일 형식 확인
- ATOM 레코드가 있는지 확인

### 낮은 RMSD에도 불구하고 active site를 찾지 못함

- Distance cutoff 증가 (예: 5.0 Å)
- 참조 active site 정의 확인
- 구조 정렬 품질 확인

### 서열 정렬에서 많은 갭 발생

- 참조 서열과 타겟 서열의 유사도 확인
- 더 적절한 참조 서열 선택 고려

## 관련 스크립트

이 통합 클래스는 다음 스크립트들의 기능을 포함합니다:

- `examples/udh_analysis/extract_udh_active_sites.py` - 서열 정렬 기반 추출
- `examples/active_site_detection/find_udh_active_sites.py` - 구조 정렬 기반 추출
- `examples/udh_analysis/extract_udh_substrate_sites.py` - Substrate binding site 추출

## 라이선스

이 프로젝트는 ASMC 프로젝트의 일부입니다.

## 참고 문헌

- BioPython: Cock et al., Bioinformatics (2009)
- Structural alignment: Kabsch algorithm
- Sequence alignment: Needleman-Wunsch algorithm

## ⚠️ PDB 번호 오프셋 주의사항

**중요 (2025-11-21 확인)**: PDB 3RFV (AtUDH 참조 구조)는 **residue 2번부터 시작**합니다.

- **오프셋 공식**: `PDB_num = alignment_pos + 2`
- **예시**: alignment pos 134 (Y) = PDB 136 (Y) ← 핵심 촉매 잔기
- **PyMOL에서 PDB 파일 사용 시**: PDB 번호 그대로 사용
- **Alignment 기반 매핑 시**: 오프셋 고려 필요

**코드에서 처리 시:**
```python
# alignment position → PDB position 변환
pdb_position = alignment_position + 2

# PDB position → alignment position 변환
alignment_position = pdb_position - 2
```

## 업데이트 이력

- 2025-11-21:
  - PDB 번호 오프셋 문서 추가 (PDB_num = alignment_pos + 2)
  - Y136 촉매 잔기 확인 (이전 Y134 오류 수정)
- 2025-11-21: 초기 릴리스
  - 구조 및 서열 기반 추출 통합
  - 배치 처리 기능 추가
  - 상세한 문서화

## 문의

문제 발생 시 GitHub Issues에 등록해주세요.
