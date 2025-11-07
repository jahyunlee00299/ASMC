# ASMC (Active Site Motif Clustering) 사용 설명서

## 📌 프로젝트 개요
ASMC는 동종 단백질 패밀리의 활성 부위(active site)를 구성하는 아미노산 다양성을 해독하기 위한 정밀 도구입니다.

## 🚀 빠른 시작

처음 사용하시는 분들은 다음 중 하나를 선택하여 빠르게 시작하세요:

1. **대화형 실행 (추천)**: `python 실행_예제.py` - 메뉴 방식으로 쉽게 실행
2. **간단한 가이드**: `Quick_Start.md` 파일 참조
3. **Windows 사용자**: `ASMC_실행.bat` 실행

## 🚀 설치 및 환경 설정

### 필수 요구사항
- Python 3.8 이상 (현재 Python 3.14 설치됨)
- 의존성 패키지들 (이미 설치 완료):
  - biopython 1.86
  - scikit-learn 1.7.2
  - pytest 8.4.2
  - plotnineseqsuite 1.2.0
  - numpy, pandas, matplotlib 등

### 설치 확인
```bash
# 패키지가 제대로 설치되었는지 확인
python -m asmc.run_asmc --help

# 테스트 실행 (선택사항)
python -m pytest tests/ -v
```

## 📖 기본 사용법

### 명령어 구조
```bash
python -m asmc.run_asmc <명령어> [옵션들]
```

### 사용 가능한 명령어
1. **run** - 메인 ASMC 워크플로우 실행
2. **identity** - 서열 유사도 계산
3. **extract** - 활성 부위 추출
4. **compare** - 활성 부위 비교
5. **unique** - 중복 제거
6. **to_xlsx** - 결과를 Excel 파일로 변환
7. **pymol** - PyMOL 시각화

## 💡 주요 워크플로우 실행 방법

### 1️⃣ 서열 데이터로부터 시작 (Homology Modeling 포함)
단백질 서열만 있고 구조가 없을 때 사용합니다.

```bash
# 프로젝트 루트 디렉토리에서 실행
python -m asmc.run_asmc run -s sequences.fasta -r udh_references.txt -o output_basic/ -t 6 --id 30
```

**사용 파일:**
- 서열 파일: `sequences.fasta` (프로젝트 루트)
- 참조 구조: `udh_references.txt` (test_data/AtUdh_pdb3rfv_chainA.pdb 포함)
- 출력: `output_basic/` 디렉토리에 생성

### 2️⃣ 이미 있는 3D 구조 모델 사용
구조 모델이 이미 있을 때 사용합니다.

```bash
# models.txt 파일 먼저 생성 필요 (아래 형식 참조)
python -m asmc.run_asmc run -m models.txt -r udh_references.txt -o output_models/ --chain A
```

**models.txt 파일 형식:**
```
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein1.pdb	AtUdh
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein2.pdb	AtUdh
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein3.pdb	AtUdh
```

### 3️⃣ Multiple Sequence Alignment(MSA)로 직접 클러스터링
MSA가 이미 준비되어 있을 때 사용합니다.

```bash
# MSA 파일이 있을 때 (예: MUSCLE, MAFFT 등으로 생성한 정렬 파일)
python -m asmc.run_asmc run -M your_msa_file.txt -o output_msa/ -e 0.3 --min-samples 5
```

### 4️⃣ 활성 부위 정렬 데이터로 서브그룹 생성
**기본 실행 (자동 파라미터):**
```bash
# UDH 활성 부위 클러스터링
python -m asmc.run_asmc run -a udh_active_sites.fasta -o udh_asmc_results/
```

**커스텀 파라미터 사용:**
```bash
# eps와 min-samples 수동 지정
python -m asmc.run_asmc run -a udh_active_sites.fasta -o udh_asmc_custom/ -e 0.65 --min-samples 3
```

**기질 결합 부위 클러스터링:**
```bash
# UDH 기질 결합 부위 분석
python -m asmc.run_asmc run -a udh_substrate_sites.fasta -o udh_substrate_results/
```

**사용 가능한 FASTA 파일:**
- `udh_active_sites.fasta` - UDH 활성 부위 서열 (107KB)
- `udh_substrate_sites.fasta` - UDH 기질 결합 부위 서열 (110KB)
- `test_data/UDHs_filtered_std2.5.fasta` - 필터링된 UDH 서열

## 📁 입력 파일 형식

### references.txt (참조 구조 파일 목록)
**예시: udh_references.txt**
```
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\AtUdh_pdb3rfv_chainA.pdb
```

**프로젝트에서 사용 가능한 PDB 파일:**
```
# UDH 참조 구조
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\AtUdh_pdb3rfv_chainA.pdb

# Tutorial 예시 파일들
C:\Users\Jahyun\PycharmProjects\ASMC\docs\tutorial\ADH4.pdb
C:\Users\Jahyun\PycharmProjects\ASMC\docs\tutorial\DH35.pdb
C:\Users\Jahyun\PycharmProjects\ASMC\docs\tutorial\DHP6.pdb
C:\Users\Jahyun\PycharmProjects\ASMC\docs\tutorial\MATA.pdb

# 테스트용 PDB 파일들
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein1.pdb
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein2.pdb
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein3.pdb
```

### sequences.fasta (타겟 서열 파일)
**프로젝트 루트의 sequences.fasta 내용:**
```
>Example_Protein_1
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Example_Protein_2
MKHLWFFLLLVAAPRWVLSAAGACGQEARPEAVGQHWEALGRFWDYLRWVQTLSEQVQEE
LLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA
```

**사용 가능한 FASTA 파일들:**
- `C:\Users\Jahyun\PycharmProjects\ASMC\sequences.fasta` - 예시 서열 2개
- `C:\Users\Jahyun\PycharmProjects\ASMC\test_data\sequences.fasta` - 테스트용 서열
- `C:\Users\Jahyun\PycharmProjects\ASMC\test_data\UDHs_filtered_std2.5.fasta` - UDH 필터링 서열
- `C:\Users\Jahyun\PycharmProjects\ASMC\docs\tutorial\sequences.fasta` - Tutorial 서열

### models.txt (모델 파일 목록)
**형식:** 각 줄에 `모델_경로[TAB]참조_이름`
```
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein1.pdb	AtUdh
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein2.pdb	AtUdh
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\protein3.pdb	AtUdh
```

### pocket.txt (선택사항 - 활성 부위 정의)
**예시: udh_pocket.txt (프로젝트에 이미 존재)**
```
C:\Users\Jahyun\PycharmProjects\ASMC\test_data\AtUdh_pdb3rfv_chainA.pdb	A	137,138,139,140,141,143,165,166,167,189,190,191,213,214,215,237,238,257,258,259
```

**형식:** `PDB_경로[TAB]체인[TAB]잔기번호(쉼표로 구분)`

## ⚙️ 주요 옵션 설명

### 기본 옵션
- **-o, --outdir**: 출력 디렉토리 (기본값: ./)
- **-t, --threads**: 사용할 CPU 스레드 수 (기본값: 6)
- **-l, --log**: 로그 파일 경로 (미지정시 stdout 출력)
- **--end**: 워크플로우 중단 지점
  - pocket: 포켓 검출 후 중단
  - modeling: 모델링 후 중단
  - alignment: 정렬 후 중단
  - clustering: 클러스터링 후 중단
  - logo: 시퀀스 로고 생성까지 완료 (기본값)

### 참조 구조 옵션
- **-r, --ref**: 참조 구조 파일 목록 (필수)
- **-p, --pocket**: 각 참조에 대한 활성 부위 정의 파일
- **--chain**: 포켓 검색할 체인 지정 (기본값: all)

### 타겟 옵션
- **-s, --seqs**: 타겟 서열 파일 (FASTA 형식)
- **-m, --models**: 타겟 모델 파일 목록
- **-M, --msa**: MSA 파일
- **-a, --active-sites**: 활성 부위 정렬 파일
- **--id**: 타겟-참조 간 최소 서열 유사도 (기본값: 30.0%)
- **-n, --nb-models**: MODELLER로 생성할 모델 수

### 클러스터링 옵션
- **-e, --eps**: 클러스터링 거리 임계값 [0,1] (기본값: auto)
- **--min-samples**: 코어 포인트 최소 샘플 수 (기본값: auto)
- **--test**: 다양한 eps 값 테스트 (0 또는 1, 기본값: 0)
- **-w, --weighted-pos**: 가중치를 둘 포켓 위치 (예: 1,6,12)

### 시퀀스 로고 옵션
- **--prefix**: 클러스터 ID 앞 접두사 (기본값: G)
- **--format**: 출력 로고 형식 (svg 또는 png, 기본값: png)
- **--resolution**: 이미지 해상도 (150 또는 300 dpi, 기본값: 300)
- **--units**: Y축 단위 (bits 또는 probability, 기본값: bits)

## 📊 출력 결과

### 출력 디렉토리 구조
```
output_dir/
├── pocket_detection/         # P2RANK 포켓 검출 결과
├── homology_modeling/        # MODELLER 모델링 결과
├── structural_alignment/     # 구조 정렬 결과
├── sequence_alignment/       # 서열 정렬 결과
├── clustering/              # 클러스터링 결과
│   ├── clusters.tsv        # 클러스터 할당 정보
│   └── dissimilarity.csv   # 거리 매트릭스
└── sequence_logos/          # 시퀀스 로고 이미지
    ├── cluster_1.png
    ├── cluster_2.png
    └── ...
```

## 🔍 추가 유틸리티 명령어

### 서열 유사도 계산
```bash
# 프로젝트 파일로 서열 유사도 확인
python -m asmc.run_asmc identity -s sequences.fasta -r udh_references.txt -o identity_results.txt

# UDH 필터링 데이터 사용
python -m asmc.run_asmc identity -s test_data/UDHs_filtered_std2.5.fasta -r udh_references.txt -o udh_identity.txt
```

### 활성 부위 추출
```bash
# 클러스터링 결과에서 특정 조건의 서열 추출
# 예: 위치 1에 Lysine(K)이 있는 서열만 추출
python -m asmc.run_asmc extract -i udh_clusters_20251104_145026.tsv -p 1 -a K -o extracted_K_position1.txt
```

### 결과를 Excel로 변환
```bash
# TSV 클러스터링 결과를 Excel 형식으로 변환
python -m asmc.run_asmc to_xlsx -i udh_clusters_20251104_145026.tsv -o udh_clusters.xlsx
```

## 🐛 문제 해결

### 일반적인 문제와 해결 방법

1. **명령어를 찾을 수 없음**
   ```bash
   # 직접 Python 모듈로 실행
   python -m asmc.run_asmc run [옵션들]
   ```

2. **메모리 부족**
   - 스레드 수를 줄여보세요: `-t 2`
   - 데이터를 작은 배치로 나누어 처리

3. **P2RANK 실행 실패**
   - pocket.txt 파일을 직접 제공하여 P2RANK 단계 건너뛰기
   - Java가 설치되어 있는지 확인

4. **클러스터링 결과가 만족스럽지 않음**
   - `--test 1` 옵션으로 다양한 eps 값 테스트
   - `--min-samples` 조정
   - `-w` 옵션으로 중요한 위치에 가중치 부여

## 💻 실제 사용 예시

### 예시 1: UDH 활성 부위 클러스터링 (가장 간단)
```bash
# 프로젝트에 이미 준비된 UDH 활성 부위 데이터 사용
python -m asmc.run_asmc run -a udh_active_sites.fasta -o udh_results/
```
**실행 결과:** `udh_results/` 디렉토리에 클러스터링 결과와 시퀀스 로고 생성

### 예시 2: 기본 서열 파일로 전체 파이프라인 실행
```bash
# 서열에서 시작하여 구조 모델링까지 수행
python -m asmc.run_asmc run -s sequences.fasta -r udh_references.txt -o output_full/ -t 6 --id 30
```
**실행 과정:**
1. 서열 유사도 계산
2. Homology modeling (MODELLER 필요)
3. 포켓 검출 (P2RANK 또는 pocket 파일 사용)
4. 구조 정렬
5. 클러스터링
6. 시퀀스 로고 생성

### 예시 3: 커스텀 파라미터로 UDH 분석
```bash
# eps 값과 최소 샘플 수를 수동으로 지정
python -m asmc.run_asmc run -a udh_substrate_sites.fasta -o udh_custom/ -e 0.25 --min-samples 2 --format svg --resolution 300
```

### 예시 4: 여러 eps 값 자동 테스트
```bash
# --test 1 옵션으로 여러 eps 값 시도
python -m asmc.run_asmc run -a udh_active_sites.fasta -o udh_test_params/ --test 1
```
**결과:** 다양한 eps 값에 대한 클러스터링 결과가 각 서브디렉토리에 저장됨

### 예시 5: 포켓 정의 파일 사용
```bash
# 미리 정의된 활성 부위(pocket) 정보 사용
python -m asmc.run_asmc run -m models.txt -r udh_references.txt -p udh_pocket.txt -o output_pocket/ --chain A
```

## 📝 참고사항

- 대용량 데이터셋의 경우 충분한 메모리와 디스크 공간 확보 필요
- 구조 모델링은 시간이 오래 걸릴 수 있음 (서열당 5-30분)
- P2RANK를 사용한 포켓 검출은 Java 설치 필요
- 시각화를 위해 PyMOL 설치 권장

## 🚀 빠른 실행 가이드 (복사 & 붙여넣기)

### 가장 간단한 실행 (추천)
```bash
python -m asmc.run_asmc run -a udh_active_sites.fasta -o quick_test/
```

### 기질 결합 부위 분석
```bash
python -m asmc.run_asmc run -a udh_substrate_sites.fasta -o substrate_analysis/
```

### 커스텀 파라미터 사용
```bash
python -m asmc.run_asmc run -a udh_active_sites.fasta -o custom_output/ -e 0.3 --min-samples 3
```

### 서열 파일로 전체 분석 (모델링 포함)
```bash
python -m asmc.run_asmc run -s sequences.fasta -r udh_references.txt -o full_pipeline/ -t 6 --id 30
```

### 서열 유사도만 확인
```bash
python -m asmc.run_asmc identity -s sequences.fasta -r udh_references.txt -o identity.txt
```

## 📁 프로젝트에서 사용 가능한 주요 파일

| 파일 유형 | 파일 경로 | 설명 |
|---------|----------|------|
| **활성 부위 FASTA** | `udh_active_sites.fasta` | UDH 활성 부위 서열 (107KB) |
| **기질 부위 FASTA** | `udh_substrate_sites.fasta` | UDH 기질 결합 부위 (110KB) |
| **예시 서열** | `sequences.fasta` | 테스트용 단백질 서열 2개 |
| **참조 PDB** | `udh_references.txt` | AtUdh PDB 구조 경로 |
| **포켓 정의** | `udh_pocket.txt` | UDH 활성 부위 잔기 정의 |
| **테스트 서열** | `test_data/UDHs_filtered_std2.5.fasta` | 필터링된 UDH 서열 |
| **Tutorial PDB** | `docs/tutorial/*.pdb` | 예시 PDB 파일들 |

## 📧 문제 발생시

테스트가 모두 통과했으므로 기본 기능은 정상 작동합니다.
추가 문제 발생시 다음을 확인하세요:
1. Python 버전 (3.8 이상, 현재 3.14 설치됨)
2. 모든 의존성 패키지 설치 여부 (`pip install -e .`)
3. 입력 파일 형식과 경로 (절대 경로 사용 권장)
4. 충분한 시스템 리소스 (메모리, 디스크 공간)
5. Windows 경로 사용 시 역슬래시(`\`) 이스케이프 확인

## 📚 관련 파일 및 스크립트

### 실행 도구
- **실행_예제.py**: 대화형 메뉴로 다양한 ASMC 기능을 쉽게 테스트
- **ASMC_실행.bat**: Windows 사용자를 위한 배치 스크립트
- **run_asmc_test.py**: 기본적인 기능 테스트 스크립트
- **run_asmc_demo.py**: 결과 시각화를 포함한 데모 스크립트

### 문서
- **Quick_Start.md**: 5분 안에 시작하기 위한 빠른 가이드
- **ASMC_사용법.md**: 이 파일, 상세한 사용 설명서
- **README.md**: 프로젝트 개요 및 설치 가이드 (영문)
- **PYTHON_3.14_COMPATIBILITY.md**: Python 3.14 호환성 정보

### 테스트
- **tests/test_asmc.py**: ASMC 주요 기능 단위 테스트
- **tests/test_utils.py**: 유틸리티 함수 단위 테스트

---
작성일: 2025년 11월 7일
ASMC 버전: 1.2.0
업데이트: 실제 프로젝트 파일 경로로 모든 예시 코드 변경