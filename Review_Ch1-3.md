# Chapter 1–3 복습 시트 (Review Sheet)

> 한 번 학습한 뒤 개념을 한눈에 재정리하기 위한 요약본.
> 대상 교재: Strang, *Introduction to Linear Algebra* 6th Ed., Ch 1–3
> 시험이 영어라서, 중요한 용어는 한글 옆에 **영문을 괄호로** 병기.

---

## 0. 한 장으로 보는 큰 그림

모든 내용이 결국 **"$A\mathbf{x} = \mathbf{b}$를 기하적·대수적으로 이해하기"** 로 귀결된다.

```
[Ch1] 벡터와 행렬              ->  선형결합, 내적/길이/각도,
      (Vectors & Matrices)       열공간, A = CR, 행렬 곱
   |
[Ch2] Ax = b 풀기              ->  소거법, LU, PA = LU,
      (Solving Ax = b)           역행렬, 전치, 대칭, 유한 차분
   |
[Ch3] 왜 그게 되는가 (구조)    ->  벡터공간·부분공간,
      (Why it works)              4대 부분공간, 차원, 기본 정리
```

**네 가지 핵심 관점 (Four Viewpoints)** — 모든 예제에서 항상 떠올릴 것:

1. **행 관점 (row picture / inner product)**: $A\mathbf{x}$의 $i$번째 성분 = $A$의 $i$행 · $\mathbf{x}$
2. **열 관점 (column picture / linear combination)**: $A\mathbf{x} = x_1\mathbf{a}_1 + \cdots + x_n\mathbf{a}_n$ ← **가장 중요**
3. **기하 관점 (geometric picture)**: 직선·평면·초평면의 교집합, 생성 (span)
4. **부분공간 관점 (subspace picture)**: $\mathbf{b} \in C(A)$ 인가? $\mathbf{x} \in N(A)$ 인가?

---

## 공통 기호·용어 (Common Notation)

| 기호 | 뜻 | 어느 공간에 사는가 |
|:---|:---|:---|
| $A \in \mathbb{R}^{m\times n}$ | $m$행 $n$열 행렬 (matrix) | — |
| $C(A)$ | 열공간 (column space) | $\mathbb{R}^m$ |
| $C(A^T)$ | 행공간 (row space) | $\mathbb{R}^n$ |
| $N(A)$ | 영공간 (null space) | $\mathbb{R}^n$ |
| $N(A^T)$ | 좌영공간 (left null space) | $\mathbb{R}^m$ |
| $r$ | 랭크 (rank) = 독립 열의 수 = 독립 행의 수 | — |
| $\mathbf{x}_p$ | $A\mathbf{x}=\mathbf{b}$의 **특수해** (particular solution) | $\mathbb{R}^n$ |
| $\mathbf{x}_n$ | **영공간 해** (null-space part), $A\mathbf{x}_n=\mathbf{0}$ | $N(A)$ |
| $\mathbf{s}_j$ | **특수해** (special solution) — $N(A)$의 기저 벡터 | $N(A)$ |
| $R_0 = \mathrm{rref}(A)$ | 기약행사다리꼴 (reduced row echelon form) | — |
| $R$ | $R_0$에서 영행 제거 | — |
| 자유변수 (free variable) | 피벗이 없는 열에 대응하는 변수 | — |
| 피벗 (pivot) | 소거 후 $U$의 0이 아닌 대각 원소 | — |

> ⚠ 주의: "**특수해**"가 두 가지 뜻 — 한글은 같지만 영어는 다름.
> - $\mathbf{x}_p$ = **particular solution** ($A\mathbf{x}=\mathbf{b}$의 한 해)
> - $\mathbf{s}_j$ = **special solution** ($A\mathbf{x}=\mathbf{0}$의 기저 벡터)

---
---

# 📘 Part 1 — Chapter 1: 벡터와 행렬 (Introduction to Vectors)

**다루는 절**: 1.1 선형결합 · 1.2 내적·길이·각도 · 1.3 행렬과 열공간 · 1.4 행렬 곱셈 $AB$와 $CR$

## 1-A. 벡터의 기본 연산

- **선형결합 (linear combination)**: $c\mathbf{v} + d\mathbf{w}$ — 스칼라곱 + 벡터 덧셈
- $\mathbf{v}, \mathbf{w} \in \mathbb{R}^2$가 일차독립 → 모든 $c\mathbf{v}+d\mathbf{w}$는 **평면 전체를 채움**
- $\mathbf{v}, \mathbf{w} \in \mathbb{R}^3$가 독립 → 결합은 **3차원 속의 2D 평면**까지만 (3차원을 채우려면 독립 벡터 3개 필요)

### ⭐ 핵심 관념: "벡터는 위치가 없다" (vectors have no position)

선형대수의 벡터는 **"방향 + 크기"** 로만 정의됨. 시작점은 벡터의 정보가 아님.

- 공간 속 서로 다른 위치에 그려진 화살표라도, **같은 방향 + 같은 크기 → 같은 벡터**
- 수학적 편의상 원점에서 출발하도록 그리는 것일 뿐
- "끝점에서 수직으로 만나야 직교" (X) → "**방향 사이 각도가 90도**" (O)

**직교 벡터 집합이 왜 3D를 채우지 않고 평면인가?**
- 위치가 정보가 아니므로 여러 "위치의 같은 방향"은 모두 **하나의 벡터**
- $\mathbf{n}=(0,0,1)$에 직교 = $z$성분만 0이면 됨 ($w_1, w_2$ 자유)
- 자유도 2 → 2차원 평면

### 영벡터의 특별한 지위 (the zero vector)

- **영벡터는 모든 벡터와 "직교"** (관례적 정의)
- $\mathbf{0}\cdot\mathbf{v}=0$은 기하학적 "각도 90°"가 아니라 $|\mathbf{0}|=0$이라서 자동 0
- 엄밀히 영벡터는 방향이 없어 "각도" 자체가 정의 안 됨
- 그래도 "**내적 = 0 ⟺ 직교**" 로 정의를 통일 (이론을 깔끔하게 만들기 위한 약속)
- **영벡터를 포함한 집합은 항상 일차종속**: $1\cdot\mathbf{0} + 0\cdot\mathbf{v} = \mathbf{0}$

## 1-B. 내적·길이·각도 (dot product, length, angle)

$$\mathbf{v}\cdot\mathbf{w} = \sum v_i w_i = \mathbf{v}^T\mathbf{w},\qquad \|\mathbf{v}\| = \sqrt{\mathbf{v}\cdot\mathbf{v}}$$

$$\cos\theta = \frac{\mathbf{v}\cdot\mathbf{w}}{\|\mathbf{v}\|\,\|\mathbf{w}\|}$$

- **슈바르츠 부등식 (Schwarz inequality)**: $|\mathbf{v}\cdot\mathbf{w}| \le \|\mathbf{v}\|\,\|\mathbf{w}\|$ (등호는 $\mathbf{v}\parallel\mathbf{w}$일 때)
- **삼각 부등식 (triangle inequality)**: $\|\mathbf{v}+\mathbf{w}\| \le \|\mathbf{v}\| + \|\mathbf{w}\|$
- **직교 (orthogonal / perpendicular)** $\Leftrightarrow$ $\mathbf{v}\cdot\mathbf{w}=0$
  → **피타고라스 (Pythagoras)**: $\|\mathbf{v}\pm\mathbf{w}\|^2 = \|\mathbf{v}\|^2+\|\mathbf{w}\|^2$
- **단위벡터 (unit vector)**: $\|\mathbf{v}\|=1$. 정규화 (normalize) = $\mathbf{v}/\|\mathbf{v}\|$
- **$\mathbf{v}\cdot\mathbf{w}$의 부호**: $\theta<90°$면 $>0$, $\theta=90°$면 $=0$, $\theta>90°$면 $<0$

> 💡 **슈바르츠의 두 가지 쓰임**:
> 1. **이론적 의미**: $\cos\theta \in [-1, 1]$이라는 사실의 다른 표현
> 2. **검증 도구**: $\cos\theta$ 계산 후 $|\mathbf{v}\cdot\mathbf{w}| \le \|\mathbf{v}\|\|\mathbf{w}\|$ 확인 → 계산 실수 방지

## 1-C. 행렬-벡터 곱 $A\mathbf{x}$의 두 관점

$A \in \mathbb{R}^{m\times n}$, $\mathbf{x} \in \mathbb{R}^n$일 때:

1. **행 관점 (row view, 내적)**: $A\mathbf{x}$의 $i$번째 성분 = $A$의 $i$행 · $\mathbf{x}$
2. **열 관점 (column view, 선형결합)**:
$$A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n$$

### 크기 규칙과 "결과가 1열이 되는 이유"

$$(m \times n) \times (n \times p) = (m \times p)$$

- 가운데 $n$은 **맞아야** 함 (안 맞으면 곱셈 불가)
- 결과의 크기는 **바깥 두 숫자**
- $\mathbf{x}$가 $n \times 1$이면 결과도 $m \times 1$ — **여러 열이 하나로 "합쳐지는" 이유**: 각 열벡터에 가중치를 두어 더하는 선형결합이라서

## 1-D. 열공간, 생성, 독립, 랭크 (기초)

- **열공간 (column space) $C(A)$**: 모든 $A\mathbf{x}$ = 열들의 모든 선형결합
- **생성 (span)**: 주어진 벡터 집합의 모든 선형결합의 집합
- **독립 (independent)**: $A\mathbf{x}=\mathbf{0}$의 유일 해가 $\mathbf{x}=\mathbf{0}$일 때
- **랭크 (rank) $r$**: 독립 열의 수 ★ **행 랭크 = 열 랭크** 항상 성립 ★
- **기저 (basis)**: $C(A)$의 기저 = $A$의 처음 $r$개의 독립 열

### ⭐ 다섯 가지 동치 조건 (시험 필수 암기)

$n \times n$ 정방 행렬 $A$에 대해 **다음은 모두 동치**:

$$\boxed{\det(A) \ne 0 \iff A \text{ 가역} \iff \text{rank}(A) = n \iff \text{열 일차독립} \iff A\mathbf{x}=\mathbf{b} \text{ 해 유일}}$$

하나 확인되면 나머지 다 성립. 반대로 $\det = 0$이면 **전부** 부정됨.

### 랭크의 여러 의미 (한 개념, 네 가지 시각)

| 시각 | 의미 |
|:---|:---|
| **대수적** | 독립 열(=행)의 개수 |
| **기하학적** | 열공간(=행공간)의 차원 |
| **분해** | $A$를 $r$개의 랭크 1 행렬의 합으로 표현 가능 |
| **정보량** | 행렬이 담은 "진짜 독립적인 정보"의 수 |

**예시 ($3\times 3$ 행렬)**:
- $C(A) = \mathbb{R}^3$ (독립 열 3개)
- $C(A)$ = 평면 (독립 열 2개)
- $C(A)$ = 직선 (독립 열 1개)
- $C(A) = \{\mathbf{0}\}$ ($A = O$)

## 1-E. 랭크 1 행렬 (rank-1 matrix)

- 모든 열벡터가 **같은 직선 위**에 놓임
- 모든 행도 같은 직선 위 (∵ 행 랭크 = 열 랭크)
- $A = \mathbf{u}\mathbf{v}^T$로 쓸 수 있음 (외적, outer product)

### 랭크 1과 "정보 압축"

$A \in \mathbb{R}^{m \times n}$이고 랭크 1이면:
- 원본: $mn$개 숫자
- 외적 분해 $\mathbf{u}\mathbf{v}^T$: $m + n$개 숫자로 충분
- **SVD / 저랭크 근사 (low-rank approximation) 의 시작점**

### ⚠ "$A_6$ vs $\mathbf{a}_1$ 은 같은가?" — 절대 아님

랭크 1 행렬 $A_6 \in \mathbb{R}^{3 \times 3}$과 그 독립 열 $\mathbf{a}_1 \in \mathbb{R}^{3 \times 1}$은:

| | $A_6$ | $\mathbf{a}_1$ |
|:---|:---|:---|
| 크기 | $3 \times 3$ | $3 \times 1$ |
| 입력 차원 | 3D 벡터 받음 | **1D 스칼라** 받음 |
| 변환 | $\mathbb{R}^3 \to \mathbb{R}^3$ (차원 뭉갬) | $\mathbb{R}^1 \to \mathbb{R}^3$ |
| $A\mathbf{x}=\mathbf{b}$ 해 | 무한히 많음 (2D 평면) | 유일 (존재 시) |
| 영공간 차원 | 2 | 0 |

→ 열공간은 같아도 **행렬 자체는 완전히 다른 객체**

## 1-F. 행렬 곱 $AB$ — 세 가지 관점 (three views)

$A \in \mathbb{R}^{m\times n}$, $B \in \mathbb{R}^{n\times p}$일 때:

1. **내적 (inner product)**: $(AB)_{ij} = \text{row}_i(A)\cdot\text{col}_j(B)$
2. **열 결합 (column combination)**: $AB = [A\mathbf{b}_1 \mid A\mathbf{b}_2 \mid \cdots \mid A\mathbf{b}_p]$
3. **외적 합 (sum of outer products)**: $AB = \sum_{k=1}^{n} \mathbf{a}_k\mathbf{b}_k^*$ ← **랭크 1 행렬의 합 (sum of rank-1 matrices)**

**곱셈 비용 (cost)**: $mnp$번 (정방이면 $n^3$) — **세 관점 모두 동일**, 계산 순서만 다름

### 내적 vs 외적 — 헷갈리지 말 것

| | 내적 (inner) | 외적 (outer) |
|:---|:---|:---|
| 기호 | $\mathbf{u}^T\mathbf{v}$ (T가 **안쪽**) | $\mathbf{u}\mathbf{v}^T$ (T가 **바깥쪽**) |
| 크기 | $(1\times n)(n\times 1) = $ 스칼라 | $(m\times 1)(1\times n) = m\times n$ 행렬 |
| 차원 변화 | 축소 (압축) | 확장 |
| 결과의 성질 | 수 하나 | **항상 랭크 1 행렬** |

> 💡 암기: **내**적 = **내** 속으로 압축 / **외**적 = **외**부로 확장

### 외적 관점의 큰 그림 (시험·개념 둘 다 중요)

$$AB = \sum_{k=1}^{n} \mathbf{a}_k \mathbf{b}_k^* \quad \text{(}n\text{개의 랭크 1 행렬의 합)}$$

**핵심 통찰**: **"랭크 $r$ 행렬 = $r$개의 랭크 1 행렬의 합으로 표현 가능"**
→ $A = CR$ 분해도 이 관점으로 재해석: $A = \sum_{k=1}^{r} \mathbf{c}_k \mathbf{r}_k^*$

**법칙**:
- 결합법칙 (associative): $(AB)C = A(BC)$ ✓
- 분배법칙 (distributive): $A(B+C) = AB + AC$ ✓
- **교환법칙 없음 (NOT commutative)**: 일반적으로 $AB \ne BA$

## 1-G. $A = CR$ 분해

임의의 행렬에 대해:

$$A = CR$$

- $C$ = $A$의 **독립 열 (independent columns)** 을 모은 $m\times r$ 행렬
- $R = (I\ F)$ 형태, $r\times n$, **기약행사다리꼴 (rref)**
- $A$의 종속 열은 $C$의 독립 열들의 결합: $A = CR = C(I\ F) = (C\ CF)$
- $C$는 $A$와 **같은 열공간**, $R$은 $A$와 **같은 행공간**

### $C$ 만드는 기계적 절차 (왼쪽부터 스캔)

1. $\mathbf{a}_1 \ne \mathbf{0}$이면 → $C$에 넣음
2. $\mathbf{a}_2$가 이미 $C$의 열의 조합이 아니면 → $C$에 넣음
3. $\mathbf{a}_3$도 마찬가지로 체크
4. $C$의 열 개수 = $r$ (랭크)

### $R$ 만들기: "각 열을 $C$로 표현하는 계수"

- $R$의 $j$번째 열 = $A$의 $j$번째 열을 $C$의 열들로 표현한 계수
- **독립 열 자리** → 단위벡터 $(0,\ldots,1,\ldots,0)^T$
- **종속 열 자리** → 실제 조합 계수
- **영벡터 열 자리** → **영벡터** $(0,0,\ldots,0)^T$ (모든 계수가 0)

### ⚠ $CR$ 분해는 유일하지 않다 (non-uniqueness)

임의의 가역 행렬 $M$에 대해:
$$A = CR = (CM)(M^{-1}R)$$

특히 랭크 1에서는 임의의 스칼라 $k \ne 0$에 대해:
$$A = \mathbf{u}\mathbf{v}^T = (k\mathbf{u})\left(\tfrac{1}{k}\mathbf{v}^T\right)$$

**교재의 "표준 분해"는 $C$의 열을 $A$에서 그대로 뽑고, $R$이 $(I\ F)$ 형태가 되도록 하는 선택**. 시험에선 특별한 지시 없으면 이 방식으로 쓰는 게 안전.

### 영벡터 열이 있으면?

영벡터는 **종속 열로 취급** ($C$에 절대 안 넣음). $R$에서는 해당 위치가 **영벡터 계수**:

$$A = \begin{pmatrix}1 & 0 & 2\\ 3 & 0 & 6\\ 2 & 0 & 4\end{pmatrix} = \underbrace{\begin{pmatrix}1 & 2\\ 3 & 6\\ 2 & 4\end{pmatrix}}_{C}\underbrace{\begin{pmatrix}1 & \mathbf{0} & 0\\ 0 & \mathbf{0} & 1\end{pmatrix}}_{R}$$

## 1-H. Ch1 자주 헷갈리는 것 (pitfalls)

| ❌ 헷갈림 | ✅ 바로잡기 |
|:---|:---|
| 행 랭크 ≠ 열 랭크 | **항상 같다** $= r$ |
| $\mathbf{v}, \mathbf{w} \in \mathbb{R}^3$의 결합이 $\mathbb{R}^3$을 채운다 | 최대 **2D 평면**까지. 3차원은 독립 벡터 3개 필요 |
| $AB = BA$ | **일반적으로 거짓**. 결합·분배만 성립 |
| 랭크 1 행렬은 열만 같은 방향 | 행도 같은 방향 (행 랭크 = 열 랭크) |
| 벡터에 "위치"가 있다 (시작점이 다르면 다른 벡터) | 벡터는 **방향+크기만**. 시작점은 정보 아님 |
| "끝점이 만나야" 직교 | **방향 각도가 90°**이면 직교 (위치 무관) |
| 영벡터는 방향이 없어서 직교 개념 없다 | 관례적으로 **모든 벡터와 직교** ($\mathbf{0}\cdot\mathbf{v}=0$) |
| 영벡터 포함 집합도 독립 가능 | **무조건 종속** |
| 내적 $\mathbf{u}^T\mathbf{v}$ vs 외적 $\mathbf{u}\mathbf{v}^T$ 헷갈림 | T 위치 확인: 안쪽=내적(스칼라), 바깥쪽=외적(행렬) |
| 랭크 1 행렬 $A \in \mathbb{R}^{3\times 3}$과 그 열 $\mathbf{a}_1 \in \mathbb{R}^{3\times 1}$은 같다 | 크기·입력 차원·해의 형태 모두 다름 |
| $A = CR$ 분해가 유일하다 | $(CM)(M^{-1}R)$도 $A$. 표준형이 관례일 뿐 |
| $R$의 "영벡터 열" 자리에 뭘 넣지? | **영벡터** (모든 계수 0) |
| "$A\mathbf{x}$ 계산하면 왜 결과가 1열?" | $\mathbf{x}$가 1열이니까. 열들이 가중치로 합쳐짐 |

## 1-I. Ch1 체크리스트

- [ ] 내적·길이·각도 공식을 바로 쓸 수 있는가?
- [ ] 슈바르츠 / 삼각 부등식을 증명 없이 쓸 수 있는가?
- [ ] 슈바르츠로 $\cos\theta$ 계산 검증할 수 있는가?
- [ ] "벡터는 위치가 없다" 개념을 말로 설명할 수 있는가?
- [ ] 영벡터와 직교·종속성의 관계를 말할 수 있는가?
- [ ] $A\mathbf{x}$를 행/열 두 관점으로 계산할 수 있는가?
- [ ] 크기 규칙 $(m\times n)(n\times p)=(m\times p)$를 외웠는가?
- [ ] $AB$의 세 가지 관점 (내적·열결합·외적합)을 모두 쓸 수 있는가?
- [ ] 내적 $\mathbf{u}^T\mathbf{v}$와 외적 $\mathbf{u}\mathbf{v}^T$을 T 위치로 구분할 수 있는가?
- [ ] $A = CR$의 $C, R$을 직접 뽑을 수 있는가?
- [ ] 영벡터 열 있는 경우 $R$에서 어떻게 처리하는지 아는가?
- [ ] 랭크 1 행렬을 보고 $A = \mathbf{u}\mathbf{v}^T$ 형태로 쓸 수 있는가?
- [ ] **다섯 가지 동치 조건** ($\det \ne 0 \iff$ 가역 $\iff$ rank $= n \iff$ 독립 $\iff$ 해 유일)을 외웠는가?

---
---

# 📗 Part 2 — Chapter 2: 연립일차방정식 풀기 (Solving Linear Equations)

**다루는 절**: 2.1 소거법·후진대입 · 2.2 소거 행렬·역행렬 · 2.3 $A=LU$ · 2.4 치환·전치 · 2.5 도함수와 유한차분

> Ch2는 주로 **정방 행렬 $A \in \mathbb{R}^{n\times n}$** 에 집중. 일반 행렬은 Ch3에서.

## 2-A. $A\mathbf{x}=\mathbf{b}$ 해의 세 가지 경우 (정방 중심)

| 경우 | 조건 | 형태 |
|:---|:---|:---|
| **유일 해 (unique)** | $A$ 가역, rank $= n$ | $\mathbf{x} = A^{-1}\mathbf{b}$ |
| **해 없음 (no solution)** | $\mathbf{b} \notin C(A)$ | — |
| **무한히 많은 해 (infinite)** | rank $<n$, $\mathbf{b}\in C(A)$ | $\mathbf{x}_p + \alpha X$ 등 |

## 2-B. 소거법 (Gaussian elimination) + 후진 대입 (back-substitution)

**소거 단계**: 행 $i$에서 행 $j$의 $l_{ij}$배를 빼서 피벗 아래에 0 만들기
$$\text{row } i \leftarrow \text{row } i - l_{ij} \cdot \text{row } j$$

**흐름**: $(A \mid \mathbf{b}) \xrightarrow{\text{소거}} (U \mid \mathbf{c}) \xrightarrow{\text{후진 대입}} \mathbf{x}$

**후진 대입**: $U\mathbf{x}=\mathbf{c}$를 마지막 행부터 위로 풀어 올라감
- $x_n = c_n / u_{nn}$, $x_{n-1} = (c_{n-1} - u_{(n-1)n}x_n)/u_{(n-1)(n-1)}$, ...

**피벗 (pivot)**: $U$ 대각의 0이 아닌 원소. 0이면 **행 교환 (row exchange)** 필요.

### 왜 행 연산이 정당한가 — 허용되는 3가지 (elementary row operations)

행 = 방정식이고, 이 세 연산은 **해집합을 바꾸지 않음**:

1. 한 행에 **0이 아닌 상수 곱하기** (식 전체에 곱하는 것과 같음)
2. 한 행에 **다른 행의 배수를 더하거나 빼기** (소거법이 쓰는 연산)
3. **두 행 바꾸기** (순서만 바뀌므로 해 불변)

> 🔑 소거법은 (2)를 반복하고, 피벗이 0이면 (3)으로 피해가는 것.

## 2-C. 소거 행렬 $E_{ij}$과 $A = LU$

- **소거 행렬 (elimination matrix)**: 단위행렬에서 $(i,j)$ 위치만 $-l_{ij}$
- **한꺼번에**: $E_{n2}\cdots E_{32} E_{31} E_{21} A = U$, 즉 $EA = U$
- **역변환**: $A = E^{-1} U = LU$

$$\boxed{A = LU}$$

- $L$ = 하삼각 (lower triangular), 대각은 1, **승수 $l_{ij}$가 그 위치에 그대로** 들어감 ★
- $U$ = 상삼각 (upper triangular), 대각은 피벗
- **성립 조건**: 행 교환 불필요 (모든 좌상단 $k\times k$ 부분행렬이 가역)

### ⚠ 승수 $l_{ij}$ 부호 규칙 (시험 실수 포인트)

**정의 기준**: "$R_i \leftarrow R_i - l_{ij} R_j$" — 항상 **빼기** 기준.

| 실제 한 연산 | 승수 값 |
|:---|:---|
| $R_i - 2R_j$ | $l_{ij} = 2$ |
| $R_i + 2R_j$ | $l_{ij} = -2$ (더하면 승수는 음수!) |
| $R_i - (-3)R_j$ = $R_i + 3R_j$ | $l_{ij} = -3$ |

**실전 계산법**: $l_{ij} = \dfrac{\text{없앨 값}}{\text{피벗 값}}$ — 이 결과를 부호 그대로 $L$의 $(i,j)$에 꽂기.

### Singular 행렬의 $LU$ — 되긴 하나?

- **운 좋은 Singular** (종속성이 뒤쪽 열에 몰림): $LU$ 형태로 **분해 가능**. 단 $U$ 대각에 0 생김 → $Ax=b$는 못 풀어
- **운 나쁜 Singular** (앞쪽 열이 종속): 소거 중간에 2열 피벗 자체가 없음 → **표준 $LU$ 불가**
- 결론: "분해가 존재하는가"와 "분해로 문제를 풀 수 있는가"는 별개. Full-rank여야 $LU$가 **실용적으로** 쓸모 있음

## 2-D. $PA = LU$ — 행 교환 필요 시

- **치환 행렬 (permutation matrix) $P$**: $I$의 행을 재배치, $n!$가지
- 성질: $P^{-1} = P^T$, $P$의 열들이 직교
- **부분 피벗팅 (partial pivoting)**: 각 열에서 가장 큰 원소를 피벗으로 → 반올림 오차 최소화, $|L_{ij}|\le 1$

## 2-E. 역행렬 (inverse) 성질

- **정의**: $A^{-1}A = AA^{-1} = I$
- **존재 조건**: $A$가 $n$개의 독립 열을 가질 때 (rank $= n$)
- **유일성**: $BA = I$이고 $AC = I$이면 $B = C$
- **역순 규칙 (reverse order)**: $(AB)^{-1} = B^{-1}A^{-1}$, $(ABC)^{-1} = C^{-1}B^{-1}A^{-1}$
- **$2\times 2$ 공식**:
$$\begin{pmatrix}a&b\\c&d\end{pmatrix}^{-1} = \frac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}$$
여기서 $ad-bc$ = **행렬식 (determinant)** $\det A$
- **삼각 행렬의 역**: 대각 원소가 모두 0이 아니면 역 존재, 역도 삼각

## 2-F. 가우스-조르단 (Gauss–Jordan) — 명시적 $A^{-1}$ 계산

$$(A \mid I) \xrightarrow{\text{소거}} (I \mid A^{-1})$$

## 2-G. 가역성 동치 조건 ★ (Invertibility — equivalent conditions)

정방 $A$에 대해 아래는 모두 **동치 (equivalent)** — 하나 성립하면 전부 성립:

- $A$ **가역 (invertible)** = 정칙 (nonsingular)
- $\mathrm{rank}(A) = n$ (full rank)
- 열이 **선형독립 (linearly independent)** / 행이 LI
- $N(A) = \{\mathbf{0}\}$ (trivial null space)
- $A\mathbf{x}=\mathbf{0}$의 유일 해가 $\mathbf{x}=\mathbf{0}$ (unique trivial solution)
- 모든 $\mathbf{b}$에 대해 $A\mathbf{x}=\mathbf{b}$가 **유일 해 (unique solution)**
- $A$가 $n$개의 0이 아닌 피벗 (nonzero pivots)을 가짐
- $\det(A) \ne 0$
- $C(A) = \mathbb{R}^n$
- $A^T$ 가역

## 2-H. 계산 비용 (computational cost)

- $A \to U$ 소거: $\tfrac{1}{3}n^3$
- $\mathbf{b} \to \mathbf{c} \to \mathbf{x}$ (전진 + 후진 대입): $n^2$
- 행렬 곱: $n^3$
- **명시적으로 $A^{-1}$을 구한 뒤 $A^{-1}\mathbf{b}$ 곱하는 건 낭비** — 소거법으로 바로 풀 것

## 2-I. 전치 (transpose) — $A^T$

- **정의**: $(A^T)_{ij} = A_{ji}$. $A$의 열 = $A^T$의 행
- **전치의 진짜 의미**: $(A\mathbf{x})\cdot\mathbf{y} = \mathbf{x}\cdot(A^T\mathbf{y})$
- **규칙**:
  - $(A+B)^T = A^T + B^T$
  - $(AB)^T = B^T A^T$ — **역순 (reverse order)**
  - $(A^{-1})^T = (A^T)^{-1}$

## 2-J. 대칭 행렬 (symmetric matrix)

- **정의**: $S^T = S$
- $S^{-1}$도 대칭
- $A^T A$, $AA^T$는 **항상 대칭**
- $A^T A$가 가역 $\iff$ $A$의 열들이 선형독립
- **$S = LDL^T$ 분해**: $S$ 대칭일 때 $LU$보다 더 자연스러움 (대칭성 보존, 소거 비용 절반)

## 2-K. 유한 차분 행렬 $K, T, B$ (Ch 2.5)

| 이름 | 형태 | 성질 |
|:---|:---|:---|
| $K$ (고정-고정, fixed-fixed) | $(-1, 2, -1)$ **삼중대각 (tridiagonal)** | 대칭, 가역, **양의 정부호 (positive definite)**; $\tfrac{1}{h^2}K\mathbf{u}=\mathbf{f}$ |
| $T$ (자유-고정, free-fixed) | 왼쪽 위 $1$로 시작 | 가역, $T = LL^T$ (촐레스키, Cholesky) |
| $B$ (자유-자유, free-free) | 양쪽 모두 자유 | **특이 (singular)**, $B\mathbf{1}=\mathbf{0}$ |

- **양의 정부호 (positive definite)**: 대칭 + $n$개의 **양의 피벗** + 모든 $\mathbf{x}\ne\mathbf{0}$에 대해 $\mathbf{x}^T K\mathbf{x}>0$
- **양의 준정부호 (positive semi-definite)**: $\mathbf{x}^T K\mathbf{x}\ge 0$, 피벗 $\ge 0$

**유한 차분 공식**:
- 전진 차분 (forward): $\frac{y(x+h)-y(x)}{h}$, $O(h)$
- 후진 차분 (backward): $\frac{y(x)-y(x-h)}{h}$, $O(h)$
- **중심 차분 (centered)**: $\frac{y(x+h)-y(x-h)}{2h}$, $O(h^2)$ ★
- 이차 차분 (second): $\frac{y(x+h)-2y(x)+y(x-h)}{h^2}$, $O(h^2)$

## 2-L. 테일러 급수 (Taylor series) — 문제 풀 때 필수

### 일반 공식 (시험지 여백에 먼저 적어둘 것)

$$y(x_0+h) = y + hy' + \frac{h^2}{2}y'' + \frac{h^3}{6}y''' + \frac{h^4}{24}y^{(4)} + \cdots$$

$$y(x_0-h) = y - hy' + \frac{h^2}{2}y'' - \frac{h^3}{6}y''' + \frac{h^4}{24}y^{(4)} - \cdots$$

> 분모는 팩토리얼 $1!, 2!, 3!, 4!, 5! = 1, 2, 6, 24, 120$
> **홀수 차수만 부호 뒤집힘** → 더하면 홀수 소거(짝수만), 빼면 짝수 소거(홀수만)

### 자주 쓰는 3개 (암기)

$$\sin h = h - \frac{h^3}{6} + \frac{h^5}{120} - \cdots\quad(\text{홀수만, 부호 교대})$$

$$\cos h = 1 - \frac{h^2}{2} + \frac{h^4}{24} - \cdots\quad(\text{짝수만, 부호 교대})$$

$$e^h = 1 + h + \frac{h^2}{2} + \frac{h^3}{6} + \frac{h^4}{24} + \cdots\quad(\text{모든 차수, 전부 +})$$

### ⭐ 정확도 판정 — 실전 3단계

1. 주어진 함수의 **Taylor 전개 쓰기**
2. **(진짜값) − (근사값)** 계산 → 앞부분이 상쇄되고 오차 항만 남음
3. **남은 항의 최저 차수 = 정확도**

### 한 줄 규칙 ★

> **"근사식이 Taylor의 앞 $k$개 항까지 포함 → 버린 첫 항의 차수가 정확도"**

**예시**:
| 근사 | 포함된 항까지 | 버린 첫 항 | 정확도 |
|:---|:---|:---|:---:|
| $\sin h \approx h$ | $h$ | $h^3$ | $O(h^3)$ |
| $\cos h \approx 1 - h^2/2$ | $h^2$ | $h^4$ | $O(h^4)$ |
| $e^h \approx 1+h$ | $h$ | $h^2$ | $O(h^2)$ |
| $e^h \approx 1+h+h^2/2$ | $h^2$ | $h^3$ | $O(h^3)$ |

> ⚠ 자주 하는 실수: "오차 식에 $h$ 있으니 1차" — **틀림**. $h$가 **몇 제곱부터** 시작하느냐가 차수. $h$가 아예 없고 $h^2$부터 시작하면 2차 정확도.

### 유한 차분 공식 정확도 (위 규칙 적용)

- 전진 차분 (forward) $\dfrac{y(x+h)-y(x)}{h}$: $O(h)$
- 후진 차분 (backward) $\dfrac{y(x)-y(x-h)}{h}$: $O(h)$
- **중심 차분 (centered)** $\dfrac{y(x+h)-y(x-h)}{2h}$: $O(h^2)$ ★ — 뺄셈으로 짝수 차수 소거
- 이차 차분 (second) $\dfrac{y(x+h)-2y(x)+y(x-h)}{h^2}$: $O(h^2)$ — 덧셈으로 홀수 차수 소거

### 다항식 대입 검증법 (Q20 유형)

**"차분 공식이 $k$차 정확도"를 보이는 법**:
- 중심점을 $x = 0$으로 잡기 (간단화)
- $u = x^2, x^4, \ldots$ 대입해 공식값 = 참값 확인
- $u = x^k$까지 정확히 맞추면 **$k$차 정확도 이상** 보장
- 4차 정확도 증명 → $u = x^2, x^4$ 둘 다 통과 확인

> 💡 대칭적 공식(계수가 반대칭)이면 **홀수 다항식은 자동 통과**, 짝수만 테스트하면 됨.

## 2-L2. 테일러·차분 문제가 선형대수 책에 있는 이유

- **미분방정식 $-u'' = f$** → **차분 공식으로 이산화** → **행렬 $K\mathbf{u}=h^2\mathbf{f}$** 형태로 등장
- $K, T, B$는 **차분 공식을 행렬로 포장한 것**
- 차분 공식의 정확도 = 행렬 해의 정확도
- 앞에서 배운 도구 ($LDL^T$, 대칭성, tridiagonal 소거의 $O(n)$ 속도 등)이 여기서 **실전 활용**

## 2-M. Ch2 자주 헷갈리는 것

| ❌ 헷갈림 | ✅ 바로잡기 |
|:---|:---|
| $(AB)^{-1} = A^{-1}B^{-1}$ | **틀림**. $(AB)^{-1} = B^{-1}A^{-1}$ — **역순** |
| $(AB)^T = A^T B^T$ | **틀림**. $(AB)^T = B^T A^T$ — **역순** |
| 피벗이 0이면 $A$는 항상 특이 | 행 교환하면 가역일 수 있음 ($PA = LU$) |
| $A = LU$ 하면 $L$ 계산이 복잡 | **승수 $l_{ij}$가 그 자리에 그대로** — 별도 계산 불필요 |
| $L$은 대각이 피벗 | $L$의 대각은 **모두 1**. 피벗은 $U$의 대각 |
| $A^{-1}$ 구해서 $\mathbf{x}=A^{-1}\mathbf{b}$가 효율적 | **낭비**. 소거로 바로 풀 것 |
| "$R_i + 2R_j$ 했으니 승수는 $+2$" | 승수는 **빼기 기준**. $R_i + 2R_j$면 $l_{ij} = -2$ |
| "$3\times 3$에서 두 열이 서로 배수 아니면 독립" | **아님**. 세 열 조합으로 종속일 수 있음 (한 열 = 다른 두 열의 결합) |
| $A^T$ (대각선 뒤집기) $= A^{-1}$ | 일반적으로 **다름**. 직교행렬일 때만 $A^T = A^{-1}$ |
| $A^{-1}$의 원소가 $A$ 원소와 단순 대응 | **무관**. 전체 계산 필요 ($2\times2$는 공식, $n\times n$은 Gauss-Jordan) |
| Singular면 $LU$ 아예 불가 | 소거가 진행되면 $LU$ 형태는 나옴. 단 $U$에 0 대각 생겨 실용성 없음 |
| "오차에 $h$ 있으면 1차 정확도" | $h$의 **최저 차수**로 판단. $h$ 항 없고 $h^2$부터면 2차 |

## 2-N. Ch2 체크리스트

- [ ] 소거법과 후진 대입으로 $A\mathbf{x}=\mathbf{b}$를 손으로 풀 수 있는가?
- [ ] $A=LU$에서 $L, U$를 바로 뽑을 수 있는가? (승수가 어디로 가는지)
- [ ] 승수 부호 규칙 ($+$ 연산이면 음수 승수)?
- [ ] $PA=LU$에서 $P$를 추적하는 방법?
- [ ] 가역성 동치 조건을 5개 이상 말할 수 있는가?
- [ ] $(AB)^{-1}$, $(AB)^T$의 역순 규칙?
- [ ] $2\times 2$ 역행렬 공식?
- [ ] 가우스-조르단으로 $A^{-1}$ 구하기?
- [ ] $K, T, B$ 행렬의 성질 구분?
- [ ] 테일러 전개 3개 ($\sin, \cos, e^h$) 암기?
- [ ] 정확도 판정 3단계 (전개 → 빼기 → 최저 차수)?
- [ ] Q20 유형: $u = x^2, x^4$ 대입으로 4차 정확도 증명?

---
---

# 📕 Part 3 — Chapter 3: 벡터 공간과 부분공간 (Vector Spaces and Subspaces)

**다루는 절**: 3.1 벡터 공간·부분공간 · 3.2 영공간 계산 · 3.3 완전해 · 3.4 독립·기저·차원 · 3.5 4대 부분공간

## 3-A. 벡터 공간 8공리 (the 8 axioms of a vector space)

$V$가 **벡터 공간 (vector space)** 이려면, $\mathbf{u},\mathbf{v},\mathbf{w}\in V$, $c,d\in\mathbb{F}$에 대해:

**덧셈 (addition)**:
1. 결합법칙 (associativity): $\mathbf{u}+(\mathbf{v}+\mathbf{w})=(\mathbf{u}+\mathbf{v})+\mathbf{w}$
2. 교환법칙 (commutativity): $\mathbf{u}+\mathbf{v}=\mathbf{v}+\mathbf{u}$
3. 영벡터 존재 (zero vector): $\exists\,\mathbf{0}$ s.t. $\mathbf{v}+\mathbf{0}=\mathbf{v}$
4. 덧셈 역원 존재 (additive inverse): $\exists\,-\mathbf{v}$ s.t. $\mathbf{v}+(-\mathbf{v})=\mathbf{0}$

**스칼라곱 (scalar multiplication)**:
5. $c(d\mathbf{v})=(cd)\mathbf{v}$
6. $1\mathbf{v}=\mathbf{v}$ (항등원, identity)
7. $c(\mathbf{u}+\mathbf{v})=c\mathbf{u}+c\mathbf{v}$ (벡터 덧셈에 대한 분배)
8. $(c+d)\mathbf{v}=c\mathbf{v}+d\mathbf{v}$ (스칼라 덧셈에 대한 분배)

## 3-B. 부분공간 판정 (subspace test) — 두 조건만

- (i) $\mathbf{u}+\mathbf{w}\in V$ — 덧셈에 닫힘 (closed under addition)
- (ii) $c\mathbf{u}\in V$ — 스칼라곱에 닫힘 (closed under scalar multiplication)
- → **영벡터 $\mathbf{0}$ 포함은 자동** (스칼라곱 조건에서 $c=0$)

### ⭐ 실전 빠른 판별 패턴 (문제 3.1.9 같은 거 30초 컷)

| 집합 정의 방식 | 부분공간? | 이유 |
|:---|:---:|:---|
| **등식 $= 0$** ($b_1 = b_2$, $b_1+b_2+b_3=0$) | ✅ 거의 YES | 원점 포함 + 선형성 보존 |
| **등식 $=$ 상수** ($b_1 = 1$, $z = 5$) | ❌ NO | **원점 미포함** |
| **부등식** ($b_1 \leq b_2$, $x \geq 0$) | ❌ NO | 음수배로 부등호 뒤집힘 |
| **곱 $= 0$** ($b_1 b_2 b_3 = 0$) | ❌ NO | 덧셈 깨짐 (축들의 합집합) |
| **"모든 선형결합"** (span) | ✅ 항상 YES | 정의상 닫힘 |

> 🔑 **원점 체크 = 첫 관문**. $\mathbf{0}$이 안 속하면 바로 NO. 시간 절약.

**부분공간인 것 / 아닌 것 예시**:
- ✅ 원점 지나는 직선/평면
- ✅ $\{\mathbf{0}\}$, $\mathbb{R}^n$ 자체
- ❌ $\mathbb{R}^3$의 평면 $z=5$ (원점 미포함)
- ❌ 제1사분면 (스칼라곱 $c<0$에서 닫히지 않음)
- ❌ 양수 벡터의 집합 (덧셈 역원 없음)
- ❌ 정수 벡터 집합 ($0.5 \times$ 정수 = 정수 아님)
- ❌ 두 축의 합집합 ($(1,0)+(0,1)=(1,1)$ 어느 축도 아님)

**$\mathbb{R}^n$ 외의 벡터공간 예**:
- **행렬 공간 (matrix space)**: 모든 $n\times n$ 행렬 → 차원 $n^2$
- **함수 공간 (function space)**: 예) $y''=0$의 해공간, 기저 $\{1, x\}$, 차원 2

### $C(A)$, $N(A)$가 어디 사는지 — 절대 헷갈리지 말 것

$A \in \mathbb{R}^{m \times n}$일 때:

$$\boxed{C(A) \subseteq \mathbb{R}^m, \quad N(A) \subseteq \mathbb{R}^n}$$

- $C(A)$: 열벡터 길이 = $m$ → $\mathbb{R}^m$
- $N(A)$: $x \in Ax=0$에서 $x$ 성분 수 = $n$ → $\mathbb{R}^n$

## 3-C. 영공간 $N(A)$ 계산 — 기약행사다리꼴 (rref)

**$A\mathbf{x}=\mathbf{0}$의 모든 해**를 구하는 게 목표.

### rref의 정확한 4가지 조건 (기억할 것)

1. 각 피벗이 **1**
2. 피벗 **위아래 전부 0** (echelon form보다 강한 조건)
3. 피벗이 왼쪽→오른쪽, 위→아래 **계단형** 배치
4. **영행은 맨 아래**

> ⚠ Ch2의 $U$ (echelon form) 는 "피벗 아래만 0". rref는 "위아래 다 0" + "피벗 = 1".

### 핵심 불변량 ⭐

**소거해도 $N(A)$는 안 바뀐다**: $N(A) = N(R_0) = N(R)$.
- 증명 핵심: $R = EA$, $E$ 가역 → $Ax=0 \iff EAx=0$
- 그래서 편한 $R$에서 $N$을 구해도 됨

### 계산 절차

1. $A$를 소거해서 $R_0 = \mathrm{rref}(A)$로 만듦
2. $R_0$은 $r$개의 **피벗 열 (pivot columns)** 과 $n-r$개의 **자유 열 (free columns)** 을 가짐
3. 순열하면 $R = (I\ F)$ 형태

### 특수해 (special solutions) 구하는 법

- **자유변수 중 하나만 1**, 나머지 자유변수 0으로 놓고 $R\mathbf{x}=\mathbf{0}$ 풀기
- 이렇게 $n-r$개 나옴 → 이들이 **$N(A)$의 기저**
- 한꺼번에 표현: 특수해의 행렬은 $\binom{-F}{I}$의 열들

$$\boxed{(I\ F)\binom{-F}{I} = O}$$

> 💡 **왜 자유변수에 "1"을 넣나?** 아무 값이나 넣어도 되지만, 1을 넣으면 패턴이 $\binom{-F}{I}$로 깔끔해지고 기저 벡터들이 **자동으로 독립** (자유변수 자리에 단위벡터 꼴이라서). 계산 편의 + 독립 보장.

### Rank-Nullity 관계

$$\boxed{\text{rank}(A) + \text{nullity}(A) = n}$$

- $\text{rank}(A) = r$ = 피벗 수 = 독립 열 수
- $\text{nullity}(A) = n - r$ = 자유변수 수 = 특수해 수 = $\dim N(A)$
- **열 개수 = 피벗 열 + 자유 열 = $r + (n-r) = n$**. 항상 성립.

## 3-D. $A\mathbf{x}=\mathbf{b}$ 완전해 (complete solution)

**완전해 공식 — 항상 이 구조**:

$$\boxed{\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n,\qquad A\mathbf{x}_p=\mathbf{b},\ A\mathbf{x}_n=\mathbf{0}}$$

- **$\mathbf{x}_p$ (특수해, particular solution)**: 모든 자유변수를 **0으로** 놓고 $R_0\mathbf{x}=\mathbf{d}$ 풀기
- **$\mathbf{x}_n$ (영공간 부분, null-space part)**: 특수해 (special solutions)들의 결합, 즉 $N(A)$의 임의 벡터
- **풀림 조건 (solvability)**: $\mathbf{b} \in C(A)$ ⟺ $R_0$의 영행에 대응하는 $\mathbf{d}$ 성분이 모두 0

### ⭐ 6단계 표준 절차 (암기 후 기계적으로)

1. **Augmented matrix** $(A \mid \mathbf{b})$ 만들기
2. **rref로 소거** → $(R_0 \mid \mathbf{d})$
3. **해 존재 확인**: $R_0$의 영행 자리에 $\mathbf{d}$도 0인지 체크 (아니면 **해 없음**, 끝)
4. **$\mathbf{x}_p$ 구하기**: 자유변수 = 0 대입 → 피벗변수 값 계산
5. **$\mathbf{x}_n$ 구하기**: 특수해들 (3.2 방식) 찾기
6. **완전해 조합**: $\mathbf{x} = \mathbf{x}_p + c_1\mathbf{s}_1 + c_2\mathbf{s}_2 + \cdots$

> 💡 **왜 $\mathbf{x}_p$는 자유변수에 0을 넣나?** 아무 값이나 가능. 0이 계산 제일 깔끔.

### 해의 존재·유일성 — 두 축 분리

| 질문 | 무엇으로 판단? |
|:---|:---|
| 해가 **존재**하는가? | $\mathbf{b} \in C(A)$ 인가? |
| 해가 **유일**한가? | $N(A) = \{\mathbf{0}\}$ 인가? |
| 해 집합의 **차원**? | $\dim N(A) = n - r$ |

- **Full column rank** ($r = n$) → 해 있을 때 유일
- **Full row rank** ($r = m$) → 모든 $\mathbf{b}$에 해 존재

### 4가지 경우 판정표 (Ch3 일반판)

| 경우 | 형태 | $R_0$ 모양 | 해의 수 |
|:---|:---|:---|:---|
| $r=m=n$ | 정방·가역 (square, invertible) | $I$ | **유일 (unique)** |
| $r=m<n$ | 짧고 넓음 (short-wide, full row rank) | $(I\ F)$ | 항상 존재, **무한히 많음 ($\infty$)** |
| $r=n<m$ | 크고 좁음 (tall-thin, full column rank) | $\binom{I}{0}$ | **0개 또는 1개** |
| $r<m,\,r<n$ | 랭크 부족 (rank deficient) | $\binom{I\ F}{0\ 0}$ | **0개 또는 무한히 많음** |

## 3-E. 선형독립·생성·기저·차원 (independence · span · basis · dimension)

- **선형독립 (linearly independent, LI)**: $c_1\mathbf{v}_1+\cdots+c_n\mathbf{v}_n=\mathbf{0}$의 유일 해가 모든 $c_i=0$
  → 동치: $A\mathbf{x}=\mathbf{0}$의 유일 해가 $\mathbf{x}=\mathbf{0}$, 즉 $N(A)=\{\mathbf{0}\}$
- **생성 (span)**: 주어진 벡터들의 모든 선형결합의 집합
- **기저 (basis)**: 독립 + 생성 → 모든 벡터가 **유일하게** 표현됨
- **차원 (dimension)**: 기저의 크기 (어떤 기저를 뽑든 불변)

**중요 사실**:
- $\mathbb{R}^m$에서 $n > m$이면 $n$개 벡터는 **반드시 종속 (dependent)**
- 모든 가역 $n\times n$ 행렬의 열들은 $\mathbb{R}^n$의 기저

## 3-F. 4대 부분공간 & 기본 정리 ★ (Fundamental Theorem) ★

$A \in \mathbb{R}^{m\times n}$, 랭크 $r$일 때:

| 부분공간 (subspace) | 정의 | 사는 곳 | 차원 (dimension) |
|:---|:---|:---|:---:|
| **열공간 (column space) $C(A)$** | $\{A\mathbf{x}\}$ | $\mathbb{R}^m$ | $r$ |
| **행공간 (row space) $C(A^T)$** | $\{A^T\mathbf{y}\}$ | $\mathbb{R}^n$ | $r$ |
| **영공간 (null space) $N(A)$** | $A\mathbf{x}=\mathbf{0}$ | $\mathbb{R}^n$ | $n-r$ |
| **좌영공간 (left null space) $N(A^T)$** | $A^T\mathbf{y}=\mathbf{0}$ | $\mathbb{R}^m$ | $m-r$ |

**직교 여공간 쌍 (orthogonal complement pairs)** — 기본 정리의 핵심:

$$\underbrace{C(A^T) \perp N(A)}_{\mathbb{R}^n\text{에서 직교 여공간}}\qquad \underbrace{C(A) \perp N(A^T)}_{\mathbb{R}^m\text{에서 직교 여공간}}$$

$$r + (n-r) = n,\qquad r + (m-r) = m$$

### ⭐ 소거의 불변량 — 왜 기저 뽑는 법이 공간마다 다른가

**핵심 원리**: 소거는 **행 연산**만 함. 그래서:

| 공간 | 소거 후 보존? | 이유 |
|:---|:---:|:---|
| 행공간 $C(A^T)$ | ✅ 보존 | 행끼리 조작해도 span은 그대로 |
| 영공간 $N(A)$ | ✅ 보존 | $Ax=0 \iff EAx=0$ (가역) |
| 열공간 $C(A)$ | ❌ **바뀜** | 열의 실제 벡터가 변함 |
| 좌영공간 $N(A^T)$ | ❌ **바뀜** | $A^T$ 기준으로는 행 연산 아님 |

**구체 예시**:
$A = \begin{pmatrix}1&2&4\\2&4&8\end{pmatrix}$ 의 열 $(1,2)$ 방향 직선 $\ne$ $R_0 = \begin{pmatrix}1&2&4\\0&0&0\end{pmatrix}$ 의 열 $(1,0)$ 방향 = $x$축.

### 기저 뽑는 법 (how to extract a basis) — 기계적 규칙

- **$C(A)$** → **원래 $A$의 피벗 열** (⚠ $R_0$에서 위치만 확인, 벡터는 $A$에서!)
- **$N(A)$** → $R_0$에서 **특수해 (special solutions)**, $\binom{-F}{I}$의 열들
- **$C(A^T)$** → $R_0$의 **영이 아닌 행** (보존되니 편한 쪽에서)
- **$N(A^T)$** → $A^T$를 따로 소거해서 특수해

> 🔑 **한 줄 요약**: "보존되는 공간은 $R$에서, 변하는 공간은 $A$ (혹은 $A^T$) 에서 뽑는다."

### $A$와 $R_0$의 관계 요약

| 공유하는 것 | 안 공유하는 것 |
|:---|:---|
| 행공간 $C(A^T) = C(R_0^T)$ | 열공간 $C(A) \ne C(R_0)$ (차원만 같음) |
| 영공간 $N(A) = N(R_0)$ | 좌영공간 $N(A^T) \ne N(R_0^T)$ |
| rank $r$ | |

## 3-G. Ch3 자주 헷갈리는 것

| ❌ 헷갈림 | ✅ 바로잡기 |
|:---|:---|
| 제1사분면, $z=5$ 평면은 부분공간 | **아니다** (원점 미포함 / 스칼라곱에 안 닫힘) |
| 소거하면 열공간도 보존됨 | **아니다**. 행공간·영공간만 보존. 열공간은 **바뀜** |
| $C(A)$ 기저를 $R_0$의 피벗 열로 뽑기 | **$A$ 자체**의 피벗 열로 뽑아야 함 |
| $\mathbf{x}_p$와 특수해 (special sol.)가 같음 | 한글은 같지만 다른 개념! $\mathbf{x}_p$는 $A\mathbf{x}=\mathbf{b}$의 한 해, special은 $N(A)$의 기저 벡터 |
| $A\mathbf{x}=\mathbf{0}$의 유일 해는 항상 $\mathbf{0}$ | $r<n$이면 **비자명 해 (nontrivial)** 있음 |
| 기저가 하나뿐 | **무한히 많음**. 차원만 불변 |
| $n > m$인 $n$개 벡터도 독립일 수 있음 | **항상 종속** |
| $b$를 $2\times 3$ 행렬 같은 걸로 생각 | $\mathbf{b}$는 **항상 벡터 하나** ($m \times 1$). $\mathbf{x}$도 벡터 ($n \times 1$) |
| "rank 2면 평면이니 해 무한" (해 개수를 $C(A)$로 판단) | 해 개수는 **$N(A)$로** 판단. $C(A)$는 해 **존재** 판정용 |
| "왜 어떨 땐 행, 어떨 땐 열로 뽑지?" | 소거는 **행 연산**이라 행 관련 공간만 보존. 열 관련은 원래 $A$ 필요 |
| 독립 = "두 벡터 사이 스칼라배 아님" | 벡터 3개 이상이면 부족. **$c_1v_1+\cdots+c_nv_n=0$의 유일 해가 모두 0** 이 정확한 정의 |
| 영벡터 포함한 집합도 독립일 수 있음 | **무조건 종속** ($1\cdot\mathbf{0} + 0\cdot v = \mathbf{0}$) |
| 완전해에서 해 존재 판단을 원래 $A$로 | 소거 후 $(R_0 \mid \mathbf{d})$에서 영행 대응 $\mathbf{d}$ 성분이 0인지 보면 됨 |

## 3-H. 용어 미묘한 차이

- **피벗 (pivot)** vs **피벗 열 (pivot column)**: 피벗 = 수 (숫자), 피벗 열 = 열 (벡터)
- **$\mathbf{x}_p$ (particular solution)** vs **$\mathbf{s}_j$ (special solution)**:
  - $\mathbf{x}_p$: $A\mathbf{x}=\mathbf{b}$ 만족, 자유변수 전부 0
  - $\mathbf{s}_j$: $A\mathbf{x}=\mathbf{0}$의 기저, 자유변수 하나만 1
- **특이 (singular)** = 비가역 / **정칙 (nonsingular)** = 가역
- **동차 (homogeneous)** = 우변 $\mathbf{0}$ / **비동차 (nonhomogeneous)** = 우변 $\mathbf{0}$ 아님
- **부정 (underdetermined)** = $n>m$ (보통 해 $\infty$) / **과결정 (overdetermined)** = $m>n$ (보통 해 없음)
- **echelon form ($U$)** (Ch2) vs **reduced row echelon form ($R_0$)** (Ch3):
  - $U$: 피벗 **아래**만 0. 피벗이 1일 필요 없음
  - $R_0$: 피벗 **위아래** 전부 0. 피벗 **= 1** 필수
- **$C(A)$ vs $C(R_0)$**: 차원은 같지만 **공간 자체는 다름** (소거로 열이 변하니까)

## 3-I. Ch3 체크리스트

- [ ] 벡터 공간 8공리를 떠올릴 수 있는가?
- [ ] 부분공간 판정 두 조건?
- [ ] 부분공간 빠른 판별 패턴 (등식/상수/부등식/곱/span)?
- [ ] $R_0 = \mathrm{rref}(A)$를 직접 구할 수 있는가? (rref 4조건)
- [ ] 왜 자유변수에 "1"을 넣는지 설명할 수 있는가?
- [ ] 특수해 (special solutions)로 $N(A)$의 기저 뽑기?
- [ ] $\mathbf{x}_p + \mathbf{x}_n$ 형태로 완전해 쓰기 (6단계 절차)?
- [ ] $A\mathbf{x}=\mathbf{b}$ 해 존재 조건 ($R_0$의 영행에 $\mathbf{d}$도 0)?
- [ ] 4가지 경우 판정표?
- [ ] 4대 부분공간의 차원 $(r, r, n-r, m-r)$?
- [ ] 4대 부분공간의 기저 각각 뽑는 법? (왜 $C(A)$만 원래 $A$에서?)
- [ ] 소거의 불변량 (보존: $C(A^T), N(A)$ / 변함: $C(A), N(A^T)$)?
- [ ] 직교 여공간 2쌍?
- [ ] 기본 정리 (Fundamental Theorem of Linear Algebra)를 한 문장으로?

---
---

## 🎯 한 줄 마무리 (the one-liner)

> **선형대수는 결국 두 질문으로 귀결된다:**
> **"$\mathbf{b}$가 $C(A)$에 있는가?"** 와 **"$N(A)$는 얼마나 큰가?"**
> (Linear algebra boils down to: "Is $\mathbf{b}$ in $C(A)$?" and "How big is $N(A)$?")
> 나머지 — 소거법, $LU$, 전치, 직교성 — 는 전부 이 두 질문에 답하기 위한 도구.

---

## 부록 A. 영문 용어 색인 (English Glossary)

영어 시험 대비용 용어 사전 — 알파벳순. 괄호 안 숫자 = 해당 챕터 (e.g. [2] = Ch2에서 주로 등장).

| English | 한글 | 비고 / 기호 |
|:---|:---|:---|
| associative [1] | 결합법칙 | $(AB)C=A(BC)$ |
| augmented matrix [2,3] | 확대(첨가) 행렬 | $(A\,|\,\mathbf{b})$ |
| back-substitution [2] | 후진 대입 | $U\mathbf{x}=\mathbf{c}$ 풀기 |
| basis [3] | 기저 | 독립 + 생성 |
| column space [1,3] | 열공간 | $C(A)$ |
| commutative [1] | 교환법칙 | $AB=BA$ (일반적으로 거짓) |
| complete solution [3] | 완전해 | $\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n$ |
| determinant [2] | 행렬식 | $\det A$ |
| diagonal matrix [2] | 대각 행렬 | $D$ |
| dimension [3] | 차원 | 기저의 크기 |
| dot (inner) product [1] | 내적 | $\mathbf{v}\cdot\mathbf{w}=\mathbf{v}^T\mathbf{w}$ |
| echelon form [3] | 사다리꼴 | |
| elimination [2] | 소거법 | 가우스 소거 |
| elimination matrix [2] | 소거 행렬 | $E_{ij}$ |
| equivalent [2] | 동치 | $\iff$ |
| free variable [3] | 자유변수 | 피벗 없는 변수 |
| full column rank [3] | 열 완전 계수 | $r=n$ |
| full row rank [3] | 행 완전 계수 | $r=m$ |
| Fundamental Theorem of Linear Algebra [3] | 선형대수학의 기본 정리 | |
| homogeneous [3] | 동차 | $A\mathbf{x}=\mathbf{0}$ |
| identity matrix [1] | 항등(단위) 행렬 | $I$ |
| independent (linearly) [1,3] | 선형독립 | LI |
| inverse [2] | 역행렬 | $A^{-1}$ |
| invertible / nonsingular [2] | 가역 / 정칙 | $\det\ne 0$ |
| left null space [3] | 좌영공간 | $N(A^T)$ |
| length (norm) [1] | 길이 (노름) | $\|\mathbf{v}\|$ |
| linear combination [1] | 선형 결합 | $c\mathbf{v}+d\mathbf{w}$ |
| lower triangular [2] | 하삼각 | $L$ |
| null space [3] | 영공간 | $N(A)$ |
| orthogonal / perpendicular [1] | 직교 / 수직 | $\mathbf{v}\cdot\mathbf{w}=0$ |
| orthogonal complement [3] | 직교 여공간 | $N(A)\perp C(A^T)$ |
| outer product [1] | 외적 | $\mathbf{u}\mathbf{v}^T$ (rank 1) |
| particular solution [3] | 특수해 (해 하나) | $\mathbf{x}_p$ |
| partial pivoting [2] | 부분 피벗팅 | 안정성 향상 |
| permutation matrix [2] | 치환(순열) 행렬 | $P$; $P^{-1}=P^T$ |
| pivot [2] | 피벗 | $U$의 0이 아닌 대각 |
| pivot column [2,3] | 피벗 열 | 피벗을 포함하는 열 |
| positive definite [2] | 양의 정부호 | $\mathbf{x}^T K\mathbf{x}>0$ |
| positive semi-definite [2] | 양의 준정부호 | $\mathbf{x}^T K\mathbf{x}\ge 0$ |
| Pythagoras theorem [1] | 피타고라스 정리 | |
| rank [1,3] | 랭크 (계수) | $r$ |
| reduced row echelon form (rref) [3] | 기약행사다리꼴 | $R_0$ |
| row space [3] | 행공간 | $C(A^T)$ |
| Schwarz inequality [1] | 슈바르츠 부등식 | $|\mathbf{v}\cdot\mathbf{w}|\le\|\mathbf{v}\|\|\mathbf{w}\|$ |
| singular [2] | 특이 (비가역) | $\det = 0$ |
| solvable [3] | 풀림 가능 | $\mathbf{b}\in C(A)$ |
| span [1,3] | 생성 | 모든 선형 결합 |
| special solution [3] | 특수해 ($N(A)$의 기저) | $N(A)$의 기저 벡터 |
| square matrix [2] | 정방(정사각) 행렬 | $m=n$ |
| subspace [3] | 부분공간 | |
| symmetric [2] | 대칭 | $S^T=S$ |
| Taylor series [2] | 테일러 급수 | |
| transpose [2] | 전치 | $A^T$ |
| triangle inequality [1] | 삼각 부등식 | $\|\mathbf{v}+\mathbf{w}\|\le\|\mathbf{v}\|+\|\mathbf{w}\|$ |
| tridiagonal [2] | 삼중대각 | $K$ |
| trivial solution [3] | 자명해 | $\mathbf{x}=\mathbf{0}$ |
| unique solution [2] | 유일 해 | $\exists!$ |
| unit vector [1] | 단위벡터 | $\|\mathbf{v}\|=1$ |
| upper triangular [2] | 상삼각 | $U$ |
| vector space [3] | 벡터 공간 | 8공리 |
| zero vector [3] | 영벡터 | $\mathbf{0}$ |
