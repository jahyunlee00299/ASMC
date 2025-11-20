# UDH 클러스터링 분석에 사용된 기법 설명

## 📊 Overview

이번 UDH 분석에서 3가지 클러스터링 기법을 사용했습니다:

1. **DBSCAN** (ASMC 사용) - 밀도 기반 클러스터링
2. **Hierarchical Clustering** - 계층적 클러스터링
3. **Identity-based Clustering** - 서열 유사도 기반 클러스터링

---

## 1️⃣ DBSCAN (Density-Based Spatial Clustering of Applications with Noise)

### 원리

**밀도 기반 클러스터링**: 데이터 포인트들이 밀집된 영역을 찾아 클러스터를 형성합니다.

```
핵심 개념:
- Core point: 반경 ε 내에 MinPts개 이상의 이웃이 있는 점
- Border point: 클러스터에 속하지만 core point는 아닌 점
- Noise point: 어떤 클러스터에도 속하지 않는 이상치
```

### 알고리즘 작동 방식

```python
for each point P:
    if P is visited:
        continue

    mark P as visited
    neighbors = get_neighbors(P, epsilon)

    if len(neighbors) < min_samples:
        mark P as NOISE
    else:
        C = new_cluster()
        expand_cluster(P, neighbors, C, epsilon, min_samples)
```

### 파라미터

1. **ε (epsilon)**: 이웃으로 간주할 최대 거리
   - 우리 분석: 0.25
   - 의미: dissimilarity 0.25 이하인 서열들끼리 이웃

2. **MinPts (min_samples)**: 클러스터 형성 최소 개수
   - 우리 분석: 2
   - 의미: 2개 이상 밀집된 곳만 클러스터

### 거리 계산 방식 (ASMC)

```python
# ASMC에서 사용하는 dissimilarity 계산
# BLOSUM62 scoring matrix 사용

dissimilarity(seq1, seq2) = 1 - (similarity_score / max_possible_score)

# 예시:
# seq1 = "SVSNHY"
# seq2 = "SVSNHY" → dissimilarity = 0 (완전 동일)
# seq2 = "AVSNHY" → dissimilarity = 0.1 (S→A 치환)
# seq2 = "XXXXXX" → dissimilarity = 1.0 (완전 다름)
```

**BLOSUM62**: 아미노산 치환 확률을 반영한 scoring matrix
- 보존적 치환 (예: I↔V): 높은 점수
- 비보존적 치환 (예: D↔K): 낮은 점수

### 장점

✅ **클러스터 개수를 미리 정하지 않아도 됨**
✅ **임의 모양의 클러스터 발견 가능**
✅ **이상치(noise) 자동 검출**
✅ **대용량 데이터에 효율적**

### 단점

❌ **밀도가 다른 클러스터 처리 어려움**
❌ **고차원 데이터에서 성능 저하**
❌ **ε와 MinPts 선택이 결과에 큰 영향**

### 우리 분석 결과

```
UDH 기질 결합 부위 (12 잔기)
- ε = 0.25
- MinPts = 2
- Silhouette score = 0.427 (매우 우수)

결과:
- Cluster 0: 849 sequences (71.9%) ← 표준 패턴
- Noise: 332 sequences (28.1%) ← 변이 패턴들
```

**해석**:
- 332개 변이가 "noise"로 분류된 이유: 서로 너무 다양해서 밀집된 영역을 형성하지 못함
- 각 변이 패턴이 1~22개씩만 있어서 MinPts=2 기준은 충족해도, ε=0.25 내에 충분한 이웃이 없음

---

## 2️⃣ Hierarchical Clustering (계층적 클러스터링)

### 원리

**상향식(Agglomerative)**: 각 데이터를 개별 클러스터로 시작해서 가장 가까운 것들을 순차적으로 합쳐나감

```
Step 1: 각 서열 = 개별 클러스터 (1181개 클러스터)
Step 2: 가장 가까운 2개 합침 (1180개 클러스터)
Step 3: 또 가장 가까운 2개 합침 (1179개 클러스터)
...
Step N: 모두 합쳐짐 (1개 클러스터)
```

### 연결 방법 (Linkage Methods)

우리가 사용한 방법: **Average Linkage**

```python
# Average linkage: 클러스터 간 평균 거리
distance(Cluster_A, Cluster_B) = mean of all pairwise distances

# 다른 방법들:
# - Single linkage: min distance (가장 가까운 점)
# - Complete linkage: max distance (가장 먼 점)
# - Ward: 분산 최소화
```

### 거리 계산 (Hamming Distance)

```python
def hamming_distance(seq1, seq2):
    """다른 위치의 개수"""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

# 예시:
seq1 = "SVSNHYIGSRMF"
seq2 = "XXSNHYIGSRXX"
#       ||        ||
# Hamming distance = 4 (위치 1,2,11,12가 다름)
```

### 클러스터 형성

```python
# Cut tree at specific distance
clusters = fcluster(linkage_matrix, t=3, criterion='distance')

# t=3: Hamming distance 3 이하인 서열들을 같은 클러스터로
```

### Dendrogram (수형도)

```
         ┌────────────────┐
         │                │
    ┌────┴────┐      ┌────┴────┐
    │         │      │         │
  ┌─┴─┐     ┌─┴─┐  ┌─┴─┐     ┌─┴─┐
  A   B     C   D  E   F     G   H

높이 = 클러스터 간 거리
가로 = 개별 서열들
```

### 장점

✅ **계층 구조 시각화 (dendrogram)**
✅ **클러스터 개수를 나중에 결정 가능**
✅ **거리 임계값만 지정하면 됨**
✅ **해석이 직관적**

### 단점

❌ **대용량 데이터에 느림 (O(n²) ~ O(n³))**
❌ **한번 합쳐진 클러스터는 분리 불가**
❌ **이상치에 민감**

### 우리 분석 결과

**변이체 재분석:**
```
332개 변이 → 130개 sub-clusters
- 가장 큰 sub-cluster: 22 sequences (XXSXHYIXSRMF)
- Distance threshold = 3 (Hamming)
```

**해석**:
- DBSCAN에서 noise로 분류된 332개를 계층적 클러스터링으로 재분석
- 유사한 변이 패턴들끼리 130개 그룹으로 세분화
- 각 그룹의 consensus를 구해 변이 패턴 특성 파악

---

## 3️⃣ Identity-based Clustering (서열 유사도 기반)

### 원리

**전체 서열 비교**: 모든 서열 쌍의 identity를 계산하여 거리 행렬 생성

```python
def calculate_identity(seq1, seq2):
    """같은 아미노산 비율"""
    min_len = min(len(seq1), len(seq2))
    identical = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    return identical / min_len

# 예시:
seq1 = "MVKVGVNG"  (8 aa)
seq2 = "MVKVGANG"  (8 aa)
#       ||||*|||
# Identity = 7/8 = 0.875 (87.5%)
```

### Distance Matrix

```python
# Identity → Distance 변환
distance = 1 - identity

# Identity matrix (500 x 500):
[[1.0  0.8  0.3  ...]  # seq1 vs all
 [0.8  1.0  0.4  ...]  # seq2 vs all
 [0.3  0.4  1.0  ...]  # seq3 vs all
 ...                 ]

# Distance matrix:
[[0.0  0.2  0.7  ...]
 [0.2  0.0  0.6  ...]
 [0.7  0.6  0.0  ...]
 ...                 ]
```

### 계층적 클러스터링 적용

```python
# scipy 사용
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# 1. Distance matrix를 condensed form으로 변환
condensed_dist = squareform(distance_matrix)

# 2. Hierarchical clustering
linkage_matrix = linkage(condensed_dist, method='average')

# 3. 클러스터 형성 (distance threshold = 0.5)
clusters = fcluster(linkage_matrix, 0.5, criterion='distance')
```

### 장점

✅ **전체 서열 정보 활용**
✅ **생물학적으로 의미 있는 distance**
✅ **정렬 없이 직접 비교 (단순)**

### 단점

❌ **길이가 다른 서열 비교 어려움**
❌ **삽입/결실 무시**
❌ **진화적 거리 반영 안 됨**

### 우리 분석 결과

**전체 서열 클러스터링 (271 aa):**
```
1181 sequences → sampled 500
Distance threshold = 0.5 (50% identity)

결과:
- 173 clusters
- Average identity = 9.9%
- 매우 다양한 서열들
```

**해석**:
- 전체 서열 수준에서는 UDH들이 매우 다양
- 평균 identity 9.9% = 진화적으로 먼 거리
- 173개 클러스터 = 다양한 UDH 아류들 존재

---

## 🔬 비교 분석

### 각 방법으로 얻은 클러스터 수

| 방법 | 잔기 수 | 클러스터 수 | 주요 클러스터 크기 | Score |
|------|---------|-------------|-------------------|-------|
| **Identity-based** | 271 aa (전체) | 173 | 다양 | - |
| **DBSCAN (ASMC)** | 20 aa (NAD/촉매) | 2 | 96.4% | 0.225 |
| **DBSCAN (ASMC)** | 12 aa (기질 결합) | 2 | 71.9% | **0.427** |
| **Hierarchical** | 12 aa (변이체만) | 130 | 다양 | - |

### 왜 다른 결과가 나왔나?

**1. 전체 서열 (173 clusters)**
- 보존되지 않은 영역 포함 → 높은 다양성
- 진화적 거리가 먼 서열들도 포함

**2. NAD/촉매 부위 (2 clusters, 96.4%)**
- 20개 잔기 중 대부분 100% 보존
- 기능적으로 필수적인 잔기들 → 변이 허용 안 됨

**3. 기질 결합 부위 (2 clusters, 71.9%)**
- 12개 잔기 중 일부는 변이 허용
- 표준 패턴 vs 변이 패턴으로 명확히 구분
- **가장 의미 있는 분류!**

**4. 변이체 재분석 (130 sub-clusters)**
- DBSCAN의 "noise"를 세분화
- 유사한 변이 메커니즘 파악

---

## 💡 왜 이 조합을 사용했나?

### 1단계: Identity-based (전체 서열)
**목적**: UDH 전체적인 다양성 파악
```python
calculate_identity_matrix(all_sequences)  # 전체 271 aa
hierarchical_clustering(threshold=0.5)     # 50% identity
```
→ 결과: 173 clusters (매우 다양)

### 2단계: DBSCAN (Active site)
**목적**: 기능적으로 중요한 잔기만 사용, 이상치 검출
```python
extract_active_sites(12_positions)         # 기질 결합 12 잔기
asmc_dbscan(epsilon=0.25, min_samples=2)  # BLOSUM62 distance
```
→ 결과: 2 clusters + 332 noise (표준 vs 변이)

### 3단계: Hierarchical (Noise 재분석)
**목적**: 변이 패턴 세분화
```python
hierarchical_clustering(noise_only)        # 332 sequences
cut_tree(distance=3)                       # Hamming distance
```
→ 결과: 130 sub-clusters (변이 메커니즘 파악)

---

## 📈 Silhouette Score (클러스터링 품질 평가)

### 정의

```python
silhouette_score = (b - a) / max(a, b)

where:
  a = 같은 클러스터 내 평균 거리 (intra-cluster)
  b = 가장 가까운 다른 클러스터까지 평균 거리 (inter-cluster)

Range: -1 to +1
- +1: 완벽한 클러스터링 (잘 분리됨)
-  0: 클러스터가 겹침
- -1: 잘못된 클러스터링
```

### 우리 결과

```
NAD/촉매 부위 (20 aa):
  Silhouette = 0.225
  → 클러스터가 다소 겹침, 분리가 명확하지 않음

기질 결합 부위 (12 aa):
  Silhouette = 0.427
  → 매우 우수한 분리!
  → 표준 vs 변이가 명확히 구분됨
```

**해석**: 기질 결합 부위 12개 잔기가 UDH 변이를 구분하는 **가장 좋은 마커**!

---

## 🎯 결론

### 최적의 클러스터링 전략

1. **DBSCAN**: 기능적 부위만 선택적으로 분석
   - 이상치 자동 검출
   - 높은 silhouette score

2. **Hierarchical**: 이상치들의 세부 패턴 분석
   - Dendrogram으로 관계 시각화
   - 변이 메커니즘 이해

3. **Identity-based**: 전체적인 다양성 파악
   - 간단하고 직관적
   - 진화적 거리 추정

### 생물학적 의미

- **전체 서열**: 다양성 (173 clusters)
- **NAD/촉매 부위**: 극도로 보존됨 (1 major cluster)
- **기질 결합 부위**: 두 가지 주요 타입 (2 clusters)
  - Type I (71.9%): 표준 `SVSNHYIGSRMF`
  - Type II (28.1%): 130가지 변이 패턴

→ **UDH의 기질 특이성은 12개 잔기로 결정됨!**

---

## 📚 참고 코드

### DBSCAN (ASMC 내부)
```python
from sklearn.cluster import DBSCAN

# Dissimilarity matrix 계산 (BLOSUM62)
dissim_matrix = compute_dissimilarities(sequences, blosum62)

# DBSCAN 클러스터링
clustering = DBSCAN(eps=0.25, min_samples=2, metric='precomputed')
labels = clustering.fit_predict(dissim_matrix)
```

### Hierarchical Clustering
```python
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist, squareform

# Hamming distance matrix
dist_matrix = compute_hamming_distances(sequences)

# Hierarchical clustering
linkage_matrix = linkage(pdist(dist_matrix), method='average')

# Cut tree
clusters = fcluster(linkage_matrix, t=3, criterion='distance')

# Visualize
dendrogram(linkage_matrix)
```

### Identity-based
```python
import numpy as np

# Calculate identity matrix
n = len(sequences)
identity_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i, n):
        identity = calculate_pairwise_identity(sequences[i], sequences[j])
        identity_matrix[i, j] = identity
        identity_matrix[j, i] = identity

# Convert to distance
distance_matrix = 1 - identity_matrix

# Apply hierarchical clustering
linkage_matrix = linkage(squareform(distance_matrix), method='average')
clusters = fcluster(linkage_matrix, t=0.5, criterion='distance')
```

---

## 🔧 파라미터 튜닝 팁

### DBSCAN
```python
# epsilon 선택 방법:
# 1. K-distance plot 확인
# 2. Median dissimilarity 참고 (우리: 0.333)
# 3. 여러 값 테스트 (0.1, 0.2, 0.25, 0.3, ...)

# min_samples 선택:
# - 작은 값 (2-3): 작은 클러스터도 검출
# - 큰 값 (5+): 큰 클러스터만 검출, 노이즈 많음
```

### Hierarchical
```python
# Distance threshold 선택:
# 1. Dendrogram에서 적절한 높이 선택
# 2. 원하는 클러스터 개수에 맞춰 조정
# 3. 생물학적 의미를 고려 (예: 3 aa 차이 = 유의미)
```

### Identity-based
```python
# Threshold 선택:
# - 0.3 (30%): 먼 관계
# - 0.5 (50%): 중간 관계
# - 0.7 (70%): 가까운 관계
# - 생물학적 기능 유사성 고려
```
