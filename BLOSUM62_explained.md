# BLOSUM62 Biological Distance 완전 해설

## 🧬 BLOSUM62란?

**BLOSUM** = **BLO**cks **SU**bstitution **M**atrix
**62** = 62% 서열 identity로 클러스터링된 데이터에서 유도

### 한 문장 요약:
> "진화적으로 관련된 단백질들을 분석하여, 각 아미노산 쌍이 얼마나 자주 서로 치환되는지를 점수화한 행렬"

---

## 📊 BLOSUM62 Matrix (일부)

```
     A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
A    4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0
R   -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3
N   -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3
D   -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3
C    0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1
Q   -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2
E   -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2
G    0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3
H   -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3
I   -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3
L   -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1
K   -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2
M   -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1
F   -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2
S    1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -1
T    0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0
W   -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3
Y   -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1
V    0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -1   0  -3  -1   4
```

### 점수 해석:

- **양수 (+)**: 진화적으로 자주 일어나는 치환 (보존적)
- **음수 (-)**: 진화적으로 드문 치환 (비보존적)
- **0**: 중립적 치환

---

## 🔬 BLOSUM62는 어떻게 만들어졌나?

### 1. 데이터 수집 (BLOCKS Database)

```
약 2,000개의 단백질 패밀리
약 500,000개의 정렬된 블록
진화적으로 관련된 단백질들
```

### 2. 서열 클러스터링 (62% identity)

```python
# 62% 이상 유사한 서열들을 하나로 묶음
if identity(seq1, seq2) >= 0.62:
    cluster_together(seq1, seq2)

# 이유: 너무 유사한 서열들이 과대 대표되는 것 방지
```

### 3. 아미노산 쌍 빈도 계산

```python
# 각 정렬 컬럼에서 아미노산 쌍 관찰
for column in alignment:
    for pair in all_pairs_in_column:
        count[pair] += 1

# 예:
# Column:  I I V I I V  →  I-I: 10쌍, I-V: 8쌍, V-V: 1쌍
```

### 4. Log-odds Score 계산

```python
# 관찰 빈도 vs 기대 빈도
score(a, b) = log2(observed_freq(a,b) / expected_freq(a,b))

# 예: I-V 치환
observed_freq(I, V) = 0.08  (8% 관찰됨)
expected_freq(I, V) = 0.05  (5% 기대됨)
score(I, V) = log2(0.08/0.05) = log2(1.6) ≈ 0.68

# 반올림하면 BLOSUM62[I][V] = 3 (양수!)
```

### 5. 정수로 스케일링

```python
# 2 bit 단위로 스케일링 (관례)
final_score = round(score * 2)
```

---

## 💡 점수의 생물학적 의미

### 높은 양수 점수 (+3 ~ +11)

```
I ↔ V  (+3)  : 소수성, 크기 유사 → 자주 치환됨
I ↔ L  (+2)  : 소수성, 분지 → 자주 치환됨
S ↔ T  (+1)  : 극성, 수산기 → 자주 치환됨
C ↔ C  (+9)  : Cysteine은 거의 보존 (이황화 결합)
W ↔ W  (+11) : Tryptophan은 극도로 보존 (크고 특수함)
```

**의미**: 구조/기능이 유사하여 서로 대체 가능

### 낮은 양수/0 점수 (0 ~ +2)

```
A ↔ S  (+1)  : 크기는 다르지만 가끔 치환
G ↔ A  (0)   : 작은 아미노산끼리 중립적
```

**의미**: 때때로 일어나지만 특별히 선호되지 않음

### 음수 점수 (-1 ~ -4)

```
D ↔ K  (-1)  : 전하 반대 (산성 vs 염기성)
I ↔ P  (-3)  : 구조 파괴 (Proline은 굴곡)
W ↔ G  (-2)  : 크기 극단적 차이
```

**의미**: 구조/기능 손상 가능, 진화적으로 회피

### 매우 낮은 음수 점수 (-3 ~ -4)

```
C ↔ R  (-3)  : Cys 이황화결합 vs Arg 염기성
W ↔ P  (-4)  : 거대+방향족 vs 작고 구조 파괴
D ↔ C  (-3)  : 산성+친수성 vs 중성+이황화결합
```

**의미**: 거의 일어나지 않는 치환, 치명적일 가능성

---

## 🎯 실제 예시: UDH 기질 결합 부위

### 표준 서열 (Cluster 0)
```
Position:  1   2   3   4   5   6   7   8   9  10  11  12
Sequence:  S   V   S   N   H   Y   I   G   S   R   M   F
```

### 변이 서열 1
```
Position:  1   2   3   4   5   6   7   8   9  10  11  12
Sequence:  X   X   S   X   H   Y   I   X   S   R   X   F
```

**Distance 계산 (BLOSUM62):**

```python
# Position 1: S → X
score += 20  # Gap penalty

# Position 2: V → X
score += 20  # Gap penalty

# Position 3: S → S
score += 4  # BLOSUM62[S][S] = 4 (동일)

# Position 4: N → X
score += 20  # Gap penalty

# Position 5: H → H
score += 8  # BLOSUM62[H][H] = 8 (동일)

# Position 6: Y → Y
score += 7  # BLOSUM62[Y][Y] = 7 (동일)

# ... 계속

Total raw score = 20+20+4+20+8+7+4+20+4+5+20+6 = 138
```

### 변이 서열 2 (보존적 치환)
```
Position:  1   2   3   4   5   6   7   8   9  10  11  12
Sequence:  S   I   S   N   H   Y   V   G   S   R   M   F
```

**Distance 계산:**

```python
# Position 1: S → S = +4
# Position 2: V → I = +3  (소수성 유사!)
# Position 3: S → S = +4
# Position 7: I → V = +3  (소수성 유사!)
# ... 나머지 동일

Total raw score = 4+3+4+6+8+7+3+6+4+5+5+6 = 61
```

**비교:**
- 변이 1 (Gap 많음): 138 (높은 거리)
- 변이 2 (보존적 치환): 61 (낮은 거리)

---

## 🔄 ASMC에서의 사용

### 1. Pairwise Score 계산 (asmc/asmc.py)

```python
def pairwise_score(scoring_dict, seqA, seqB, weighted_pos):
    score = 0
    for i, (posA, posB) in enumerate(zip(seqA, seqB)):
        if posA in ["-", "X"] or posB in ["-", "X"]:
            # Gap penalty
            if i+1 in weighted_pos:
                score += 20 * 5  # 가중 위치
            else:
                score += 20
        else:
            # BLOSUM62 점수 (dissimilarity이므로 양수 = 거리)
            if i+1 in weighted_pos:
                score += scoring_dict[posA][posB] * 5
            else:
                score += scoring_dict[posA][posB]

    return score
```

### 2. Dissimilarity로 변환

```python
# BLOSUM62는 similarity (높을수록 유사)
# ASMC는 dissimilarity 필요 (높을수록 다름)

# 방법 1: 점수를 음수로
dissimilarity = -1 * blosum_score

# 방법 2: MinMaxScaler로 정규화
dissimilarity = normalize(blosum_score)
```

ASMC 코드를 보면 **raw BLOSUM score를 그대로 사용**하고 있습니다!

```python
score += scoring_dict[posA][posB]  # BLOSUM62 값 그대로
```

이것이 작동하는 이유:
- BLOSUM62에서 **낮은 점수 = 비유사 = 높은 dissimilarity**
- 음수가 많을수록 다른 서열
- MinMaxScaler가 나중에 정규화

### 3. MinMaxScaler 정규화

```python
# 모든 거리를 0-1 범위로 정규화
data = MinMaxScaler().fit_transform(data)

# 예:
# Raw scores: [10, 50, 100, 150]
# Normalized: [0, 0.29, 0.64, 1.0]
```

**이점:**
- DBSCAN의 epsilon이 일정한 의미를 가짐
- 서열 길이에 독립적
- 0.25 = "최대 거리의 25%"

---

## 📊 다른 Substitution Matrix와 비교

### BLOSUM 시리즈

| Matrix | Identity | 용도 |
|--------|----------|------|
| **BLOSUM45** | 45% | 먼 진화적 거리 (원거리 상동성) |
| **BLOSUM62** | 62% | **범용** (가장 널리 사용) |
| **BLOSUM80** | 80% | 가까운 진화적 거리 (최근 분기) |

**규칙:**
- 숫자가 **낮을수록** → 먼 진화적 관계 → 보수적인 점수
- 숫자가 **높을수록** → 가까운 진화적 관계 → 엄격한 점수

### PAM 시리즈

| Matrix | Mutations | 용도 |
|--------|-----------|------|
| **PAM30** | 30 per 100 aa | 가까운 관계 |
| **PAM120** | 120 per 100 aa | 중간 |
| **PAM250** | 250 per 100 aa | 먼 관계 |

**차이점:**
- PAM: 진화 모델 기반 (extrapolation)
- BLOSUM: 실제 관찰 데이터 (empirical)

### BLOSUM62 vs PAM120 비교

```
I → V 치환:

BLOSUM62:  I-V = +3  (자주 관찰됨)
PAM120:    I-V = +1  (덜 허용적)

D → K 치환:

BLOSUM62:  D-K = -1  (드물지만 가능)
PAM120:    D-K = -5  (매우 불리)
```

---

## 🧮 왜 BLOSUM62가 기본값인가?

### 1. 실증적 성능

```
BLAST, FASTA, HMMER 등에서 기본값
수천 개의 벤치마크에서 검증됨
```

### 2. 중간 거리 최적화

```
너무 가까운 거리 (BLOSUM80): 너무 엄격
너무 먼 거리 (BLOSUM45): 너무 관대
BLOSUM62: 대부분의 경우에 적합
```

### 3. 생물학적 의미

```
62% identity = 단백질 패밀리 내 전형적인 거리
기능은 보존, 서열은 다양화된 상태
```

### 4. 역사적 이유

```
1992년 Henikoff & Henikoff가 제안
즉시 성능이 입증됨
표준으로 자리잡음
```

---

## 💻 실제 BLOSUM62 적용 예시

### 예제 1: 보존적 치환

```python
seq1 = "ILE"  # Isoleucine
seq2 = "VAL"  # Valine

# 둘 다 소수성, 분지형, 크기 유사

distance = 0
distance += BLOSUM62['I']['V']  # +3
distance += BLOSUM62['L']['A']  # +2
distance += BLOSUM62['E']['E']  # +5

total = 3 + 2 + 5 = 10
# 낮은 거리 = 유사한 서열
```

### 예제 2: 비보존적 치환

```python
seq1 = "ASP"  # Aspartate (산성)
seq2 = "LYS"  # Lysine (염기성)

distance = 0
distance += BLOSUM62['A']['L']  # -1
distance += BLOSUM62['S']['Y']  # -2
distance += BLOSUM62['P']['S']  # -1

total = -1 + -2 + -1 = -4
# 음수 점수가 누적되면 매우 다른 서열
```

### 예제 3: UDH 실제 데이터

```python
# Standard: SVSNHYIGSRMF
# Variant:  AVSNHYIGSRMF
#           |           (S→A 차이)

score = 0

# Position 1: S → A
score += BLOSUM62['S']['A']  # +1 (작고 극성 유사)

# Position 2-12: 모두 동일
score += BLOSUM62['V']['V']  # +4
score += BLOSUM62['S']['S']  # +4
# ... 나머지 +4, +6, +8, +7, +4, +6, +4, +5, +5, +6

total = 1 + 4 + 4 + 6 + 8 + 7 + 4 + 6 + 4 + 5 + 5 + 6 = 60
# 60 / 12 = 5.0 (평균 유사도 높음)
```

---

## 🎨 시각화: BLOSUM62 Heatmap

```
             Hydrophobic          Polar    Charged
             I L V M F W Y  S T N Q  C G P  K R H  D E
Hydrophobic
I (4)        █ █ █ █ ▓ ▒ ▒  ░ ░ ░ ░  ▓ ▒ ░  ░ ░ ░  ░ ░
L (4)        █ █ █ █ ▓ ▒ ▒  ░ ░ ░ ░  ▓ ▒ ░  ░ ░ ░  ░ ░
V (4)        █ █ █ █ ▓ ▒ ▒  ░ ░ ░ ░  ▓ ▒ ░  ░ ░ ░  ░ ░

Polar
S (4)        ░ ░ ░ ░ ▒ ░ ▒  █ █ █ ▓  ▓ ▓ ░  ▓ ░ ░  ▓ ▓
T (5)        ░ ░ ░ ░ ▒ ░ ▒  █ █ █ ▓  ▓ ▓ ░  ▓ ░ ░  ▓ ▓

Charged
K (5)        ░ ░ ░ ░ ░ ░ ░  ▓ ▓ ▓ █  ▒ ░ ░  █ █ █  ░ █
R (5)        ░ ░ ░ ░ ░ ░ ░  ▓ ▓ ▓ █  ▒ ░ ░  █ █ █  ░ █

Special
C (9)        ▓ ▓ ▓ ▓ ▒ ▒ ▒  ▓ ▓ ░ ░  █ ░ ░  ░ ░ ░  ░ ░
W (11)       ▒ ▒ ▒ ▒ █ █ █  ░ ░ ░ ░  ▒ ░ ░  ░ ░ ░  ░ ░

Legend:
█ = Highly favorable (+4 to +11)
▓ = Favorable (+2 to +3)
▒ = Neutral (0 to +1)
░ = Unfavorable (-1 to -2)
  = Highly unfavorable (-3 to -4)
```

**패턴:**
- 같은 그룹 내: 높은 점수 (대각선 블록)
- 다른 그룹 간: 낮은 점수 (off-diagonal)
- Cys, Trp: 매우 특수함

---

## 🔍 BLOSUM62의 한계

### 1. 문맥 무시
```
BLOSUM62는 각 위치를 독립적으로 취급
실제로는 인접 잔기가 영향을 줌
```

### 2. 구조 정보 없음
```
2차 구조 (helix, sheet, loop)를 고려하지 않음
같은 치환이라도 위치에 따라 다른 영향
```

### 3. 기능적 제약 무시
```
Active site vs 표면 vs 내부
모든 위치를 동등하게 취급
```

### 4. Position-specific 정보 부족
```
Profile HMM이나 PSSM이 더 정확
BLOSUM62는 일반적인 경향만 반영
```

---

## 💡 ASMC의 개선점

ASMC는 BLOSUM62의 한계를 보완:

### 1. Position Weighting
```python
if i+1 in weighted_pos:
    score *= 5  # 중요한 위치에 가중치
```

### 2. Active Site Focusing
```python
# 전체 서열 대신 active site만 사용
# 기능적으로 중요한 잔기만 비교
```

### 3. Gap Penalty
```python
if posA == 'X' or posB == 'X':
    score += 20  # 구조적 차이 반영
```

---

## 📚 정리

### BLOSUM62 = 생물학적 지식이 담긴 distance metric

**장점:**
✅ 진화적으로 검증된 치환 확률
✅ 생물학적 의미 반영
✅ 광범위한 검증과 사용
✅ 단순하고 효율적

**단점:**
❌ 문맥과 구조 무시
❌ Position-specific 정보 부족
❌ 모든 서열에 동일하게 적용

**ASMC에서의 역할:**
🎯 Active site 잔기의 생물학적 유사도 측정
🎯 단순 Hamming distance보다 의미 있는 클러스터링
🎯 DBSCAN의 입력으로 사용

---

## 🎓 결론

**BLOSUM62 biological distance**는:

1. **진화 데이터에서 학습**한 아미노산 치환 점수
2. **생물학적 유사성**을 수치화
3. **서열 비교의 기본 도구**
4. **ASMC의 핵심 요소** - 단순 DBSCAN을 생물학적으로 의미 있게 만듦

단순한 문자열 비교가 아니라 **생물학적 지식을 활용한 intelligent distance**입니다!
