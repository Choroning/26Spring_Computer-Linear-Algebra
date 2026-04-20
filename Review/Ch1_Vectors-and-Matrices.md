# Chapter 1 복습 시트 — 벡터와 행렬 (Introduction to Vectors)

> Strang, *Introduction to Linear Algebra* 6th Ed., Chapter 1
> 다루는 절: 1.1 선형결합 · 1.2 내적·길이·각도 · 1.3 행렬과 열공간 · 1.4 행렬 곱셈 $AB$와 $CR$

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
| $R_0 = \mathrm{rref}(A)$ | 기약행사다리꼴 (reduced row echelon form) | — |

---

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

## 부록. Ch1 영문 용어 색인 (English Glossary)

| English | 한글 | 비고 / 기호 |
|:---|:---|:---|
| associative | 결합법칙 | $(AB)C=A(BC)$ |
| column space | 열공간 | $C(A)$ |
| commutative | 교환법칙 | $AB=BA$ (일반적으로 거짓) |
| dot (inner) product | 내적 | $\mathbf{v}\cdot\mathbf{w}=\mathbf{v}^T\mathbf{w}$ |
| identity matrix | 항등(단위) 행렬 | $I$ |
| independent (linearly) | 선형독립 | LI |
| length (norm) | 길이 (노름) | $\|\mathbf{v}\|$ |
| linear combination | 선형 결합 | $c\mathbf{v}+d\mathbf{w}$ |
| orthogonal / perpendicular | 직교 / 수직 | $\mathbf{v}\cdot\mathbf{w}=0$ |
| outer product | 외적 | $\mathbf{u}\mathbf{v}^T$ (rank 1) |
| Pythagoras theorem | 피타고라스 정리 | |
| rank | 랭크 (계수) | $r$ |
| Schwarz inequality | 슈바르츠 부등식 | $|\mathbf{v}\cdot\mathbf{w}|\le\|\mathbf{v}\|\|\mathbf{w}\|$ |
| span | 생성 | 모든 선형 결합 |
| triangle inequality | 삼각 부등식 | $\|\mathbf{v}+\mathbf{w}\|\le\|\mathbf{v}\|+\|\mathbf{w}\|$ |
| unit vector | 단위벡터 | $\|\mathbf{v}\|=1$ |
