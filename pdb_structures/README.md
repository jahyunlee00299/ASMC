# UDH Active Site Finder

이 디렉토리는 AtUdh (Arabidopsis thaliana UDP-glucose dehydrogenase)의 구조와 active site 정보를 이용하여 다른 UDH 단백질들의 active site를 찾기 위한 도구들을 포함합니다.

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
