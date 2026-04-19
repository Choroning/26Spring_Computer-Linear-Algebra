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

## 1-B. 내적·길이·각도 (dot product, length, angle)

$$\mathbf{v}\cdot\mathbf{w} = \sum v_i w_i = \mathbf{v}^T\mathbf{w},\qquad \|\mathbf{v}\| = \sqrt{\mathbf{v}\cdot\mathbf{v}}$$

$$\cos\theta = \frac{\mathbf{v}\cdot\mathbf{w}}{\|\mathbf{v}\|\,\|\mathbf{w}\|}$$

- **슈바르츠 부등식 (Schwarz inequality)**: $|\mathbf{v}\cdot\mathbf{w}| \le \|\mathbf{v}\|\,\|\mathbf{w}\|$ (등호는 $\mathbf{v}\parallel\mathbf{w}$일 때)
- **삼각 부등식 (triangle inequality)**: $\|\mathbf{v}+\mathbf{w}\| \le \|\mathbf{v}\| + \|\mathbf{w}\|$
- **직교 (orthogonal / perpendicular)** $\Leftrightarrow$ $\mathbf{v}\cdot\mathbf{w}=0$
  → **피타고라스 (Pythagoras)**: $\|\mathbf{v}\pm\mathbf{w}\|^2 = \|\mathbf{v}\|^2+\|\mathbf{w}\|^2$
- **단위벡터 (unit vector)**: $\|\mathbf{v}\|=1$. 정규화 (normalize) = $\mathbf{v}/\|\mathbf{v}\|$
- **$\mathbf{v}\cdot\mathbf{w}$의 부호**: $\theta<90°$면 $>0$, $\theta=90°$면 $=0$, $\theta>90°$면 $<0$

## 1-C. 행렬-벡터 곱 $A\mathbf{x}$의 두 관점

$A \in \mathbb{R}^{m\times n}$, $\mathbf{x} \in \mathbb{R}^n$일 때:

1. **행 관점 (row view, 내적)**: $A\mathbf{x}$의 $i$번째 성분 = $A$의 $i$행 · $\mathbf{x}$
2. **열 관점 (column view, 선형결합)**:
$$A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n$$

## 1-D. 열공간, 생성, 독립, 랭크 (기초)

- **열공간 (column space) $C(A)$**: 모든 $A\mathbf{x}$ = 열들의 모든 선형결합
- **생성 (span)**: 주어진 벡터 집합의 모든 선형결합의 집합
- **독립 (independent)**: $A\mathbf{x}=\mathbf{0}$의 유일 해가 $\mathbf{x}=\mathbf{0}$일 때
- **랭크 (rank) $r$**: 독립 열의 수 ★ **행 랭크 = 열 랭크** 항상 성립 ★
- **기저 (basis)**: $C(A)$의 기저 = $A$의 처음 $r$개의 독립 열

**예시 ($3\times 3$ 행렬)**:
- $C(A) = \mathbb{R}^3$ (독립 열 3개)
- $C(A)$ = 평면 (독립 열 2개)
- $C(A)$ = 직선 (독립 열 1개)
- $C(A) = \{\mathbf{0}\}$ ($A = O$)

## 1-E. 랭크 1 행렬 (rank-1 matrix)

- 모든 열벡터가 **같은 직선 위**에 놓임
- 모든 행도 같은 직선 위 (∵ 행 랭크 = 열 랭크)
- $A = \mathbf{u}\mathbf{v}^T$로 쓸 수 있음 (외적, outer product)

## 1-F. 행렬 곱 $AB$ — 세 가지 관점 (three views)

$A \in \mathbb{R}^{m\times n}$, $B \in \mathbb{R}^{n\times p}$일 때:

1. **내적 (inner product)**: $(AB)_{ij} = \text{row}_i(A)\cdot\text{col}_j(B)$
2. **열 결합 (column combination)**: $AB = [A\mathbf{b}_1 \mid A\mathbf{b}_2 \mid \cdots \mid A\mathbf{b}_p]$
3. **외적 합 (sum of outer products)**: $AB = \sum_{k=1}^{n} \mathbf{a}_k\mathbf{b}_k^*$ ← **랭크 1 행렬의 합 (sum of rank-1 matrices)**

**곱셈 비용 (cost)**: $mnp$번 (정방이면 $n^3$)

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

## 1-H. Ch1 자주 헷갈리는 것 (pitfalls)

| ❌ 헷갈림 | ✅ 바로잡기 |
|:---|:---|
| 행 랭크 ≠ 열 랭크 | **항상 같다** $= r$ |
| $\mathbf{v}, \mathbf{w} \in \mathbb{R}^3$의 결합이 $\mathbb{R}^3$을 채운다 | 최대 **2D 평면**까지. 3차원은 독립 벡터 3개 필요 |
| $AB = BA$ | **일반적으로 거짓**. 결합·분배만 성립 |
| 랭크 1 행렬은 열만 같은 방향 | 행도 같은 방향 (행 랭크 = 열 랭크) |

## 1-I. Ch1 체크리스트

- [ ] 내적·길이·각도 공식을 바로 쓸 수 있는가?
- [ ] 슈바르츠 / 삼각 부등식을 증명 없이 쓸 수 있는가?
- [ ] $A\mathbf{x}$를 행/열 두 관점으로 계산할 수 있는가?
- [ ] $AB$의 세 가지 관점 (내적·열결합·외적합)을 모두 쓸 수 있는가?
- [ ] $A = CR$의 $C, R$을 직접 뽑을 수 있는가?
- [ ] 랭크 1 행렬을 보고 $A = \mathbf{u}\mathbf{v}^T$ 형태로 쓸 수 있는가?

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

## 2-C. 소거 행렬 $E_{ij}$과 $A = LU$

- **소거 행렬 (elimination matrix)**: 단위행렬에서 $(i,j)$ 위치만 $-l_{ij}$
- **한꺼번에**: $E_{n2}\cdots E_{32} E_{31} E_{21} A = U$, 즉 $EA = U$
- **역변환**: $A = E^{-1} U = LU$

$$\boxed{A = LU}$$

- $L$ = 하삼각 (lower triangular), 대각은 1, **승수 $l_{ij}$가 그 위치에 그대로** 들어감 ★
- $U$ = 상삼각 (upper triangular), 대각은 피벗
- **성립 조건**: 행 교환 불필요 (모든 좌상단 $k\times k$ 부분행렬이 가역)

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

$$\sin h = h - \frac{h^3}{6} + \frac{h^5}{120} - \cdots$$

$$\cos h = 1 - \frac{h^2}{2} + \frac{h^4}{24} - \cdots$$

$$e^h = 1 + h + \frac{h^2}{2} + \frac{h^3}{6} + \frac{h^4}{24} + \cdots$$

**쓰임새**:
- 유한 차분 공식의 **절단 오차 차수 (truncation order)** 판정 ($O(h)$ vs $O(h^2)$): $y(x\pm h)$를 테일러 전개해서 더하거나 빼 보면 어떤 항이 남는지 보임
- $\sin$은 홀함수 (odd) → 홀수 차수만, $\cos$는 짝함수 (even) → 짝수 차수만
- **중심 차분이 $O(h^2)$인 이유**: $y(x+h)-y(x-h)$에서 짝수 차수 전부 소거, 홀수 차수만 남음 → 주 오차 $h^3$
- $e^h$는 해석 함수 (entire / analytic) — 모든 차수 존재

## 2-M. Ch2 자주 헷갈리는 것

| ❌ 헷갈림 | ✅ 바로잡기 |
|:---|:---|
| $(AB)^{-1} = A^{-1}B^{-1}$ | **틀림**. $(AB)^{-1} = B^{-1}A^{-1}$ — **역순** |
| $(AB)^T = A^T B^T$ | **틀림**. $(AB)^T = B^T A^T$ — **역순** |
| 피벗이 0이면 $A$는 항상 특이 | 행 교환하면 가역일 수 있음 ($PA = LU$) |
| $A = LU$ 하면 $L$ 계산이 복잡 | **승수 $l_{ij}$가 그 자리에 그대로** — 별도 계산 불필요 |
| $L$은 대각이 피벗 | $L$의 대각은 **모두 1**. 피벗은 $U$의 대각 |
| $A^{-1}$ 구해서 $\mathbf{x}=A^{-1}\mathbf{b}$가 효율적 | **낭비**. 소거로 바로 풀 것 |

## 2-N. Ch2 체크리스트

- [ ] 소거법과 후진 대입으로 $A\mathbf{x}=\mathbf{b}$를 손으로 풀 수 있는가?
- [ ] $A=LU$에서 $L, U$를 바로 뽑을 수 있는가? (승수가 어디로 가는지)
- [ ] $PA=LU$에서 $P$를 추적하는 방법?
- [ ] 가역성 동치 조건을 5개 이상 말할 수 있는가?
- [ ] $(AB)^{-1}$, $(AB)^T$의 역순 규칙?
- [ ] $2\times 2$ 역행렬 공식?
- [ ] 가우스-조르단으로 $A^{-1}$ 구하기?
- [ ] $K, T, B$ 행렬의 성질 구분?
- [ ] 테일러 급수로 유한 차분 오차 차수 도출?

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
- → **영벡터 $\mathbf{0}$ 포함은 자동**

**부분공간인 것 / 아닌 것 예시**:
- ✅ 원점 지나는 직선/평면
- ✅ $\{\mathbf{0}\}$, $\mathbb{R}^n$ 자체
- ❌ $\mathbb{R}^3$의 평면 $z=5$ (원점 미포함)
- ❌ 제1사분면 (스칼라곱 $c<0$에서 닫히지 않음)
- ❌ 양수 벡터의 집합 (덧셈 역원 없음)

**$\mathbb{R}^n$ 외의 벡터공간 예**:
- **행렬 공간 (matrix space)**: 모든 $n\times n$ 행렬 → 차원 $n^2$
- **함수 공간 (function space)**: 예) $y''=0$의 해공간, 기저 $\{1, x\}$, 차원 2

## 3-C. 영공간 $N(A)$ 계산 — 기약행사다리꼴 (rref)

**$A\mathbf{x}=\mathbf{0}$의 모든 해**를 구하는 게 목표.

1. $A$를 소거해서 $R_0 = \mathrm{rref}(A)$로 만듦
2. $R_0$은 $r$개의 **피벗 열 (pivot columns)** 과 $n-r$개의 **자유 열 (free columns)** 을 가짐
3. 순열하면 $R = (I\ F)$ 형태

**특수해 (special solutions) 구하는 법**:
- 자유변수 하나만 1, 나머지 0으로 놓고 $R\mathbf{x}=\mathbf{0}$ 풀기
- 이렇게 $n-r$개 나옴 → 이들이 **$N(A)$의 기저**
- 한꺼번에 표현: 특수해의 행렬은 $\binom{-F}{I}$의 열들

$$\boxed{(I\ F)\binom{-F}{I} = O}$$

## 3-D. $A\mathbf{x}=\mathbf{b}$ 완전해 (complete solution)

**완전해 공식 — 항상 이 구조**:

$$\boxed{\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n,\qquad A\mathbf{x}_p=\mathbf{b},\ A\mathbf{x}_n=\mathbf{0}}$$

- **$\mathbf{x}_p$ (특수해, particular solution)**: 모든 자유변수를 **0으로** 놓고 $R_0\mathbf{x}=\mathbf{d}$ 풀기
- **$\mathbf{x}_n$ (영공간 부분, null-space part)**: 특수해 (special solutions)들의 결합, 즉 $N(A)$의 임의 벡터
- **풀림 조건 (solvability)**: $\mathbf{b} \in C(A)$ ⟺ $R_0$의 영행에 대응하는 $\mathbf{d}$ 성분이 모두 0

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

**기저 뽑는 법 (how to extract a basis)**:
- $C(A)$의 기저 → **$A$의 피벗 열** (⚠ $R_0$의 피벗 열이 **아님!**)
- $C(A^T)$의 기저 → $R_0$의 영이 아닌 행 (or $R$의 모든 행)
- $N(A)$의 기저 → **특수해 (special solutions)**, $\binom{-F}{I}$의 열들
- $N(A^T)$의 기저 → $A^T\mathbf{y}=\mathbf{0}$의 특수해

**$A$와 $R_0$의 관계 요약**:
- $A$와 $R_0$는 **행공간·영공간 동일** (소거가 행공간·영공간을 보존)
- 열공간은 **다름** ($C(A) \ne C(R_0)$). 차원만 같음

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

## 3-H. 용어 미묘한 차이

- **피벗 (pivot)** vs **피벗 열 (pivot column)**: 피벗 = 수 (숫자), 피벗 열 = 열 (벡터)
- **$\mathbf{x}_p$ (particular solution)** vs **$\mathbf{s}_j$ (special solution)**:
  - $\mathbf{x}_p$: $A\mathbf{x}=\mathbf{b}$ 만족, 자유변수 전부 0
  - $\mathbf{s}_j$: $A\mathbf{x}=\mathbf{0}$의 기저, 자유변수 하나만 1
- **특이 (singular)** = 비가역 / **정칙 (nonsingular)** = 가역
- **동차 (homogeneous)** = 우변 $\mathbf{0}$ / **비동차 (nonhomogeneous)** = 우변 $\mathbf{0}$ 아님
- **부정 (underdetermined)** = $n>m$ (보통 해 $\infty$) / **과결정 (overdetermined)** = $m>n$ (보통 해 없음)

## 3-I. Ch3 체크리스트

- [ ] 벡터 공간 8공리를 떠올릴 수 있는가?
- [ ] 부분공간 판정 두 조건?
- [ ] $R_0 = \mathrm{rref}(A)$를 직접 구할 수 있는가?
- [ ] 특수해 (special solutions)로 $N(A)$의 기저 뽑기?
- [ ] $\mathbf{x}_p + \mathbf{x}_n$ 형태로 완전해 쓰기?
- [ ] 4가지 경우 판정표?
- [ ] 4대 부분공간의 차원 $(r, r, n-r, m-r)$?
- [ ] 4대 부분공간의 기저 각각 뽑는 법?
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
