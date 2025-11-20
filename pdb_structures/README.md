# UDH Active Site Finder & Superposition Tools

이 디렉토리는 AtUdh (Arabidopsis thaliana UDP-glucose dehydrogenase)의 구조와 active site 정보를 이용하여 다른 UDH 단백질들의 active site를 찾고, active site 기반 superposition을 수행하는 도구들을 포함합니다.

## 제공 도구

### 1. `find_udh_active_sites.py`
- 전체 구조 정렬을 통한 active site 탐색
- 여러 UDH 구조 일괄 처리
- ASMC pocket 형식 변환 지원

### 2. `extract_and_align_active_sites.py` ⭐ NEW
- **Active site만을 사용한 구조 정렬**
- 5개 특정 UDH AlphaFold 구조 분석
- Active site RMSD vs. Global RMSD 비교
- 정렬된 구조 PDB 파일 저장

## 필요한 파일

### 1. AtUdh 참조 구조
- **파일명**: `AtUdh_3rfv.pdb`
- **다운로드**: https://www.rcsb.org/structure/3RFV
- 또는 다음 명령어로 다운로드:
  ```bash
  wget https://files.rcsb.org/download/3RFV.pdb -O AtUdh_3rfv.pdb
  ```

### 2. Active Site 정의 파일
- **파일명**: `atudh_active_site.txt`
- **형식**: Chain과 residue 번호를 포함
- 예시가 이미 생성되어 있습니다 (문헌 기반)

### 3. 비교할 다른 UDH 구조들
- `other_udh_structures/` 디렉토리에 PDB 파일들을 추가하세요

## 사용 방법

### 1단계: 참조 구조 준비
```bash
# AtUdh 구조를 이 디렉토리에 복사하거나 다운로드
cp /path/to/your/AtUdh_3rfv.pdb .
```

### 2단계: Active Site 찾기
```bash
# 단일 UDH 구조의 active site 찾기
python find_udh_active_sites.py --reference AtUdh_3rfv.pdb \
                                 --reference-sites atudh_active_site.txt \
                                 --target other_udh.pdb \
                                 --output results/

# 여러 UDH 구조들 한번에 처리
python find_udh_active_sites.py --reference AtUdh_3rfv.pdb \
                                 --reference-sites atudh_active_site.txt \
                                 --target-dir other_udh_structures/ \
                                 --output results/
```

### 3단계: 결과 확인
결과는 다음 형식으로 저장됩니다:
- `results/alignment_summary.txt`: 정렬 요약
- `results/[target_name]_active_sites.txt`: 각 타겟의 active site residues
- `results/active_sites_comparison.tsv`: 모든 구조의 active site 비교

## 출력 파일 형식

### active_sites_comparison.tsv
```
Structure	Chain	Residue_Positions	Residue_Types
AtUdh_3rfv	A	10,34,135,137,139,...	T,D,S,N,Y,...
Target_UDH	A	12,36,137,139,141,...	T,D,S,N,Y,...
```

## AtUdh Active Site 정보

UDP-glucose dehydrogenase의 주요 active site regions:
1. **NAD+ binding site**: Rossmann fold에 위치
2. **Substrate binding site**: UDP-glucose가 결합하는 부위
3. **Catalytic residues**: 산화 반응을 촉매하는 잔기들

자세한 정보는 `atudh_active_site.txt` 파일을 참조하세요.

## Active Site Superposition (새로운 기능!)

### 사용법: 5개 AlphaFold UDH 구조 분석

```bash
# Active site만을 사용하여 구조 정렬
python extract_and_align_active_sites.py \
    --reference AtUdh_3rfv.pdb \
    --reference-sites atudh_active_site.txt \
    --target-dir uronate_dehydrogenase_alphafold/ \
    --output active_site_results/
```

### 대상 구조 (Rossmann RMSD 기준)
- A0A1I2LZE5 (15.78 Å) - 가장 유사
- A0A1I1AQL9 (17.27 Å)
- A0A2Z4AB20 (19.85 Å)
- A0A1N6JJV2 (20.89 Å)
- A0A1U7CQA8 (25.19 Å) - 가장 다름

### 출력 결과
- `active_site_alignment_summary.tsv`: 전체 RMSD 비교
- `[ID]_active_site_alignment.txt`: 개별 상세 결과
- `[ID]_aligned.pdb`: Active site로 정렬된 구조 (PyMOL로 시각화 가능)

자세한 사용법은 `ACTIVE_SITE_SUPERPOSITION_GUIDE.md`를 참조하세요.

## 문서

- **QUICK_START.md**: 5분 빠른 시작 가이드
- **USAGE_EXAMPLE.md**: 상세 사용 예제
- **ACTIVE_SITE_SUPERPOSITION_GUIDE.md**: Active site superposition 완전 가이드
- **TEST_RESULTS.md**: Dehydrogenase 테스트 결과
