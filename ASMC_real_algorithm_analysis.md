# ASMC 실제 알고리즘 분석

## 🚨 중요: ASMC는 HMM을 사용하지 않습니다!

코드를 직접 확인한 결과, **ASMC는 표준 sklearn DBSCAN을 사용**합니다.

---

## 📝 실제 구현 확인 (asmc/asmc.py)

### 1. Dissimilarity 계산 (asmc.py:538-578)

```python
def pairwise_score(scoring_dict, seqA, seqB, weighted_pos):
    """BLOSUM62 기반 pairwise distance 계산"""
    score = 0
    for i, (posA, posB) in enumerate(zip(seqA, seqB)):
        if posA in ["-", "X"] or posB in ["-", "X"]:
            # Gap/Unknown에 페널티
            if i+1 in weighted_pos:
                score += 20 * 5  # 가중 위치
            else:
                score += 20
        else:
            # BLOSUM62 점수 사용
            if i+1 in weighted_pos:
                score += scoring_dict[posA][posB] * 5
            else:
                score += scoring_dict[posA][posB]

    return score
```

**특징:**
- BLOSUM62 substitution matrix 사용
- Position-specific weighting 지원
- Gap penalty: 20
- Weighted position penalty: 20 × 5 = 100

### 2. Distance Matrix 생성 (asmc.py:580-624)

```python
def dissimilarity(sequences, scoring_dict, weighted_pos):
    """All-vs-all distance matrix 계산"""
    data = []
    for key1 in sequences:
        row = []
        for key2 in sequences:
            if key1 == key2:
                score = 0.0
            else:
                score = pairwise_score(scoring_dict,
                                      sequences[key1],
                                      sequences[key2],
                                      weighted_pos)
            row.append(score)
        data.append(row)

    # MinMaxScaler로 0-1 정규화
    data = preprocessing.MinMaxScaler().fit_transform(
        data.reshape(-1,1)
    ).reshape(data.shape)

    return data
```

**특징:**
- O(n²) 복잡도의 all-vs-all 비교
- **MinMaxScaler로 0-1 범위 정규화** ← 중요!
- 정규화 후 값이 dissimilarity (거리)

### 3. DBSCAN 클러스터링 (asmc.py:626-653)

```python
def dbscan_clustering(data, threshold, min_samples, threads):
    """sklearn DBSCAN 사용"""
    dbscan = DBSCAN(eps=threshold,
                    metric="precomputed",  # 이미 계산된 거리 행렬 사용
                    n_jobs=threads,
                    min_samples=min_samples)

    labels = dbscan.fit_predict(X=data)

    return labels
```

**특징:**
- **sklearn.cluster.DBSCAN 그대로 사용**
- `metric="precomputed"`: 직접 계산한 distance matrix 사용
- 특별한 변형 없음

---

## ❌ HMM이 아닌 이유

### HMM(Hidden Markov Model)이란?

```
HMM은:
- 상태(state) 간 전이 확률 모델
- Profile HMM: 서열의 각 위치를 상태로 모델링
- 예: HMMER, SAM, pfam 등
- 용도: 서열 검색, 도메인 찾기, 상동성 검출
```

**HMM 특징:**
- Training: 다중 정렬로부터 확률 모델 학습
- Emission probability: 각 위치에서 아미노산 발생 확률
- Transition probability: 상태 간 전이 확률
- Viterbi algorithm: 최적 경로 찾기

### ASMC가 사용하는 방법

```
ASMC는:
- Pairwise comparison (쌍별 비교)
- BLOSUM62 substitution matrix
- DBSCAN density-based clustering
- 확률 모델 없음
- Training 과정 없음
```

---

## 🔍 ASMC의 특별한 점은?

ASMC는 일반 DBSCAN과 달리 **생물학적으로 의미 있는 distance**를 사용합니다:

### 1. BLOSUM62 기반 Distance

```python
# 일반 DBSCAN
distance = euclidean_distance(seq1, seq2)

# ASMC
distance = blosum62_based_dissimilarity(seq1, seq2)
```

**BLOSUM62의 의미:**
- 진화적으로 보존된 치환: 작은 페널티
  - I ↔ V (소수성 유사): 높은 점수 → 작은 거리
  - S ↔ T (극성 유사): 높은 점수 → 작은 거리
- 진화적으로 드문 치환: 큰 페널티
  - D ↔ K (전하 반대): 낮은 점수 → 큰 거리
  - W ↔ G (크기 차이): 낮은 점수 → 큰 거리

### 2. Position Weighting

```python
# 중요한 잔기에 5배 가중치
if i+1 in weighted_pos:
    score += blosum62[posA][posB] * 5
```

**의미:**
- 촉매 잔기, 기질 결합 잔기 등에 더 큰 가중치
- 기능적으로 중요한 위치의 차이를 강조

### 3. Gap Penalty

```python
if posA in ["-", "X"] or posB in ["-", "X"]:
    score += 20  # or 100 for weighted positions
```

**의미:**
- Alignment gap 또는 불확실한 잔기에 페널티
- 구조적 차이 반영

### 4. MinMaxScaler 정규화

```python
data = MinMaxScaler().fit_transform(data)
# → 모든 거리를 0-1 범위로 정규화
```

**의미:**
- DBSCAN의 epsilon 값이 일정한 의미를 가짐
- 0.25 = 최대 거리의 25%
- 서열 길이에 독립적

---

## 📊 ASMC vs 일반 클러스터링

| 특징 | 일반 K-means | 일반 DBSCAN | **ASMC** |
|------|--------------|-------------|----------|
| 알고리즘 | K-means | DBSCAN | **DBSCAN** |
| Distance | Euclidean | Euclidean | **BLOSUM62** |
| 생물학적 의미 | ❌ | ❌ | **✅** |
| Position weighting | ❌ | ❌ | **✅** |
| Active site focus | ❌ | ❌ | **✅** |
| Noise detection | ❌ | ✅ | **✅** |
| 클러스터 개수 | 사전 지정 | 자동 | **자동** |

---

## 🎯 왜 사용자가 HMM이라고 생각했을까?

### 가능한 이유:

1. **다른 논문과 혼동**
   - ASMC 외의 active site 분석 논문에서 HMM 사용
   - 예: Pfam, PROSITE 등은 Profile HMM 사용

2. **USalign의 특징**
   - ASMC는 구조 정렬에 USalign 사용
   - USalign은 TM-score 기반 (HMM 아님)

3. **용어 혼동**
   - "Modeling"이라는 용어 때문에 HMM으로 오해
   - 실제로는 Homology Modeling (MODELLER)

---

## 📖 ASMC의 전체 워크플로우

```
1. Homology Modeling (MODELLER)
   ├─ 서열 → 3D 구조 예측
   └─ 참조 구조 기반

2. Pocket Detection (P2RANK)
   ├─ Ligand binding pocket 자동 탐지
   └─ 또는 사용자 지정 pocket

3. Structural Alignment (USalign)
   ├─ Active site 구조 정렬
   └─ TM-score 계산

4. Sequence Extraction
   ├─ 정렬된 active site 서열 추출
   └─ 우리가 한 것: 직접 서열 정렬로 추출

5. Dissimilarity Matrix (BLOSUM62)
   ├─ Pairwise distance 계산
   └─ MinMaxScaler 정규화

6. DBSCAN Clustering
   ├─ Density-based clustering
   └─ Noise detection

7. Sequence Logo Generation
   └─ 각 클러스터의 consensus 시각화
```

**우리가 사용한 부분**: Step 4-7
- 구조 정보 없이 서열 정렬로 active site 추출
- ASMC의 dissimilarity + DBSCAN 사용

---

## 💡 결론

### ASMC는:
- ✅ **sklearn DBSCAN** 사용
- ✅ **BLOSUM62** 기반 생물학적 distance
- ✅ **Position weighting** 지원
- ✅ Active site에 특화
- ❌ **HMM 사용 안 함**

### 특별한 점:
ASMC의 혁신은 **새로운 클러스터링 알고리즘이 아니라**:
1. Active site만 선택적 분석
2. 생물학적으로 의미 있는 distance (BLOSUM62)
3. 구조 + 서열 정보 통합
4. 자동화된 워크플로우

단순 DBSCAN이지만, **적절한 distance metric**과 **적절한 feature selection** (active site)으로 생물학적으로 의미 있는 클러스터링을 달성!

---

## 📚 관련 알고리즘 비교

### Profile HMM (예: HMMER)
```
용도: 서열 검색, 도메인 탐지
방법:
- Multiple alignment → Profile
- Position-specific emission probability
- Viterbi algorithm
특징: 확률적, generative model
```

### ASMC (이 프로젝트)
```
용도: Active site 클러스터링
방법:
- Active site 추출
- BLOSUM62 dissimilarity
- DBSCAN clustering
특징: 거리 기반, discriminative
```

### 차이점:
- HMM: "이 서열이 이 패밀리에 속할 확률은?"
- ASMC: "이 서열들이 얼마나 유사한가?"

완전히 다른 목적과 방법!

---

## 🔬 혼동하기 쉬운 이유

### 공통점:
- ✅ 서열 분석
- ✅ 다중 서열 다룸
- ✅ 생물학적 의미 고려

### 차이점:
| | HMM | ASMC |
|---|-----|------|
| **모델 타입** | 확률적 | 거리 기반 |
| **Training** | 필요 (alignment) | 불필요 |
| **출력** | 확률, 점수 | 클러스터 |
| **주 용도** | 검색, 분류 | 클러스터링 |
| **알고리즘** | Viterbi, Forward | DBSCAN |

---

## 📌 요약

**질문**: "ASMC가 HMM 사용하나요?"

**답변**:
❌ **아니요!** ASMC는 **sklearn DBSCAN + BLOSUM62 dissimilarity**를 사용합니다.

**ASMC의 핵심**:
1. 생물학적 distance (BLOSUM62)
2. Active site focusing
3. Standard DBSCAN

단순하지만 효과적인 접근법입니다!
