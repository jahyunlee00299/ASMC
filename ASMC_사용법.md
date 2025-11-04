# ASMC (Active Site Motif Clustering) 사용 설명서

## 📌 프로젝트 개요
ASMC는 동종 단백질 패밀리의 활성 부위(active site)를 구성하는 아미노산 다양성을 해독하기 위한 정밀 도구입니다.

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
python -m asmc.run_asmc run \
    -s sequences.fasta \          # 타겟 서열 파일 (FASTA 형식)
    -r references.txt \           # 참조 구조 파일 목록
    -o output_dir/ \              # 출력 디렉토리
    -t 6 \                        # CPU 스레드 수
    --id 30                       # 최소 서열 유사도 (%)
```

### 2️⃣ 이미 있는 3D 구조 모델 사용
구조 모델이 이미 있을 때 사용합니다.

```bash
python -m asmc.run_asmc run \
    -m models.txt \               # 모델 파일 경로 목록
    -r references.txt \           # 참조 구조 파일 목록
    -o output_dir/ \              # 출력 디렉토리
    --chain A,B                   # 분석할 체인 지정
```

### 3️⃣ Multiple Sequence Alignment(MSA)로 직접 클러스터링
MSA가 이미 준비되어 있을 때 사용합니다.

```bash
python -m asmc.run_asmc run \
    -M msa_file.txt \             # MSA 파일
    -o output_dir/ \              # 출력 디렉토리
    -e 0.3 \                      # 클러스터링 거리 임계값
    --min-samples 5               # 최소 샘플 수
```

### 4️⃣ 활성 부위 정렬 데이터로 서브그룹 생성
```bash
python -m asmc.run_asmc run \
    -a active_sites.fasta \       # 활성 부위 정렬 파일
    -o output_dir/                # 출력 디렉토리
```

## 📁 입력 파일 형식

### references.txt (참조 구조 파일 목록)
```
/path/to/protein1.pdb
/path/to/protein2.pdb
/path/to/protein3.pdb
```

### sequences.fasta (타겟 서열 파일)
```
>Protein1
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA
>Protein2
MKHLWFFLLLVAAPRWVLSAAGACGQEARPEAVGQHWEALGRFWDYLRWVQTLSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA
```

### models.txt (모델 파일 목록)
```
/path/to/model1.pdb	reference1
/path/to/model2.pdb	reference1
/path/to/model3.pdb	reference2
```

### pocket.txt (선택사항 - 활성 부위 정의)
```
reference1.pdb	A	10,15,20,25,30
reference2.pdb	B	12,17,22,27,32
```

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
python -m asmc.run_asmc identity \
    -s sequences.fasta \
    -r references.txt \
    -o identity_results.txt
```

### 활성 부위 추출
```bash
python -m asmc.run_asmc extract \
    -i input_file.tsv \
    -p 1 \                    # 위치 (1-3)
    -a K \                    # 아미노산 타입
    -o extracted.txt
```

### 결과를 Excel로 변환
```bash
python -m asmc.run_asmc to_xlsx \
    -i clustering_results.tsv \
    -o results.xlsx
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

### 예시 1: 간단한 테스트 실행
```bash
# 테스트 데이터 준비
echo "/path/to/test.pdb" > refs.txt
echo ">TestSeq
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQR" > seqs.fasta

# ASMC 실행
python -m asmc.run_asmc run \
    -s seqs.fasta \
    -r refs.txt \
    -o test_output/ \
    --end modeling
```

### 예시 2: 전체 파이프라인 실행
```bash
python -m asmc.run_asmc run \
    -s my_proteins.fasta \
    -r reference_structures.txt \
    -o full_analysis/ \
    -t 8 \
    --id 25 \
    -e auto \
    --min-samples auto \
    --format svg \
    --resolution 300
```

## 📝 참고사항

- 대용량 데이터셋의 경우 충분한 메모리와 디스크 공간 확보 필요
- 구조 모델링은 시간이 오래 걸릴 수 있음 (서열당 5-30분)
- P2RANK를 사용한 포켓 검출은 Java 설치 필요
- 시각화를 위해 PyMOL 설치 권장

## 📧 문제 발생시

테스트가 모두 통과했으므로 기본 기능은 정상 작동합니다.
추가 문제 발생시 다음을 확인하세요:
1. Python 버전 (3.8 이상)
2. 모든 의존성 패키지 설치 여부
3. 입력 파일 형식과 경로
4. 충분한 시스템 리소스

---
작성일: 2025년 11월 4일
ASMC 버전: 1.2.0